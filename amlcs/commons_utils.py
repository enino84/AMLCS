import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as spa
import pandas as pd
import time
from netCDF4 import Dataset
from sklearn.linear_model import Ridge, Lasso

####################################################################################
####################################################################################
####################################################################################
def compute_modified_Cholesky_decomposition(DX, pre_info, thr=0.15):
    local_pre = pre_info[0];
    npr = pre_info[1];
    n,N = DX.shape;
    D_sqrt = np.zeros(n);
    non_zeros = n+npr; #predecessors + diagonal elements
    I = np.zeros(non_zeros);
    J = np.zeros(non_zeros);
    V = np.zeros(non_zeros);
    #if test==4: print('npr = '+str(npr)+' n^2= '+str(n**2)+' non_zeros = '+str(non_zeros));
    ind = 0;
    for i in range(0,n):
        pre_i = np.array(local_pre[i]).astype('int32');
        #print('i {0} pre_i {1}'.format(i,pre_i));
        xi = DX[i,:].T;
        pi = DX[pre_i].T;
        if pre_i.size>0:
           beta_i = compute_coef_SVD(pi, xi, thr);
           ##if test==4: print(str(i) + ' - ' + str(pre_i.size) + '('+str(ind)+','+str(ind+pre_i.size)+')');
           I[ind:ind+pre_i.size] = i;
           J[ind:ind+pre_i.size] = pre_i;
           V[ind:ind+pre_i.size] = -beta_i;
           std_i = np.std(xi-pi @ beta_i);
           D_sqrt[i] = 1/std_i if std_i > 1e-8 else 1/1e-8;
           ind += pre_i.size;
        else:
           #D_sqrt[i] = -1;
           std_i = np.std(xi);
           D_sqrt[i] = 1/std_i if std_i > 1e-8 else 1/1e-8;
    I[ind:ind+n] = np.arange(0,n);
    J[ind:ind+n] = np.arange(0,n);
    V[ind:ind+n] = np.ones(n);
    					
    #if thr==0.16: pd.DataFrame(np.concatenate((D_sqrt,I,J,V), axis=0)).to_csv('All_OBS_2.csv');
    
    #L factor
    L = spa.coo_matrix((V, (I, J)), shape=(n, n));
    
    In = np.arange(0,n);
    
    #D^{-1/2} factor
    Dinv = spa.coo_matrix((D_sqrt,(In, In)), shape=(n,n));
    
    Binv_sqrt = L.T @ Dinv;
    return Binv_sqrt;

####################################################################################
####################################################################################
####################################################################################
def get_random_code():
    return str(int(pd.to_datetime('now').timestamp()));
####################################################################################
####################################################################################
####################################################################################
def _robust_svd(A):
    """SVD with a fallback for the LAPACK 'SVD did not converge' failure.

    np.linalg.svd uses the divide-and-conquer driver (gesdd), which can fail
    to converge on ill-conditioned local design matrices (common when a wide
    radius packs many predecessors against few ensemble members). On failure
    we retry with scipy's slower but more robust gesvd driver. If that also
    fails, we return a zero spectrum so the caller degrades gracefully
    (the affected row keeps only its diagonal entry instead of crashing).
    """
    # sanitize non-finite entries that would otherwise break LAPACK
    if not np.all(np.isfinite(A)):
        A = np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0);
    try:
        return np.linalg.svd(A, full_matrices=False);
    except np.linalg.LinAlgError:
        try:
            import scipy.linalg as sla;
            return sla.svd(A, full_matrices=False, lapack_driver='gesvd');
        except Exception:
            N, n = A.shape;
            mn = min(N, n);
            return (np.zeros((N, mn)), np.zeros(mn), np.zeros((mn, n)));


def compute_coef_SVD(A,b,thr):
    N,n = A.shape;
    Ui,Si,Vi = _robust_svd(A);
    if Si.size == 0 or np.max(Si) == 0:
        return np.zeros(n);
    Smax = np.max(Si);
    Vi = Vi.T;
    beta = np.zeros(n);
    minn = min(N,n);
    for i in range(0,minn):
        #dxi = Vi[:,i]*((Ui[:,i].T @ b)/Si[i]);
        #print('* shape '+str(beta.shape));
        if Si[i]/Smax>thr:
           beta += Vi[:,i]*((Ui[:,i].T @ b)/Si[i]);
        else:
           #if i>20: print('* Hago break y me salgo en {0} de {1}'.format(i,minn))
           break;
    return beta;


####################################################################################
####################################################################################
####################################################################################
#  SHRINKAGE BASED ON INVERSE COVARIANCE (precision-space shrinkage)
#
#  These helpers implement a family of estimators of the background-error
#  *precision* matrix B^{-1} as a convex combination of two precision
#  estimators:
#
#        Binv_shrunk(alpha) = alpha * Binv_(1) + (1-alpha) * Binv_(2),
#
#  with alpha in [0,1] chosen by one of three principled criteria:
#
#     - "mse"   : Frobenius mean-squared-error weight. The two estimators
#                 are two modified-Cholesky precisions at two different
#                 radii of influence (r_narrow < r_wide). Both are sparse,
#                 banded, full-rank and SPD.
#     - "stein" : Stein-loss (information-geometric) weight. Combines the
#                 modified-Cholesky precision (radius r) with the
#                 truncated-SVD pseudo-inverse precision  W W^T.
#     - "da"    : data-assimilation-aware weight minimising the trace of the
#                 posterior analysis covariance. Same two estimators as Stein.
#
#  This is the observation/precision-space version consistent with:
#  Nino-Ruiz, E. D., Sandu, A., & Deng, X. (2018). An ensemble Kalman filter
#  implementation based on modified Cholesky decomposition for inverse
#  covariance matrix estimation. SIAM J. Sci. Comput., 40(2), A867-A886.
####################################################################################
####################################################################################
####################################################################################
def compute_modified_Cholesky_precision(DX, pre_info, thr=0.15):
    """Modified-Cholesky precision factor at a given predecessor radius.

    Returns the sparse lower/upper factor ``Binv_sqrt`` such that the
    estimated precision is ``Binv = Binv_sqrt @ Binv_sqrt.T``. This is just
    a thin alias around :func:`compute_modified_Cholesky_decomposition` so
    the shrinkage code reads in precision-space terms.
    """
    return compute_modified_Cholesky_decomposition(DX, pre_info, thr=thr);


####################################################################################
####################################################################################
####################################################################################
def _fro2_sparse(A):
    """Squared Frobenius norm of a (sparse or dense) matrix."""
    if spa.issparse(A):
        d = A.tocsr().data;
        return float((d**2).sum());
    return float((np.asarray(A)**2).sum());


####################################################################################
####################################################################################
####################################################################################
def _estimate_bl_decay(Binv, r):
    """Estimate (beta, C) for the Bickel-Levina banded-precision bound.

    The off-diagonal mass at lag l of the precision Binv is assumed to decay
    like C * l^{-beta}. We read it off by fitting log|mass(l)| vs log(l) for
    lags 1..r, with robust fallbacks when the fit is ill-determined.
    """
    A = Binv.toarray() if spa.issparse(Binv) else np.asarray(Binv);
    lags = [];
    mass = [];
    for ell in range(1, int(r)+1):
        diag_vals = np.abs(np.diagonal(A, offset=-ell));
        if diag_vals.size == 0:
            continue;
        m = float(diag_vals.mean());
        if m > 0:
            lags.append(ell);
            mass.append(m);
    if len(lags) >= 2:
        x = np.log(np.array(lags));
        y = np.log(np.array(mass));
        slope = np.polyfit(x, y, 1)[0];
        beta = float(np.clip(-slope, 0.5, 4.0));
        C = float(np.exp(y[0] + beta*x[0]));
    else:
        beta = 1.0;
        C = float(mass[0]) if mass else float(np.abs(np.diagonal(A)).mean());
    if not np.isfinite(C) or C <= 0:
        C = 1.0;
    return beta, C;


####################################################################################
####################################################################################
####################################################################################
def compute_alpha_mse(Binv1, Binv2, n, r1, r2):
    """Frobenius-MSE shrinkage weight between two modified-Cholesky targets.

    Binv1 is the narrow-radius (r1) target: more bias, less variance (stable).
    Binv2 is the wide-radius  (r2) target: less bias, more variance.

    The first-order optimum is alpha* = -rho/gamma with gamma=||Binv1-Binv2||_F^2
    and |rho| bounded by the Bickel-Levina banded-precision rate. Lacking the
    sign of rho we take the admissible weight closest to 1 (favouring the
    stable narrow target): alpha = clip_{[0,1]}(1 - B_max/gamma).
    """
    gamma = _fro2_sparse(Binv1 - Binv2);
    if gamma <= 0.0:
        return 1.0;
    beta, C = _estimate_bl_decay(Binv2, r2);
    B_max = (C**2) * n * (r2**(-beta)) * (r1**(-beta) + r2**(-beta));
    half_width = B_max / gamma;
    alpha = 1.0 - half_width;
    return float(min(1.0, max(0.0, alpha)));


####################################################################################
####################################################################################
####################################################################################
def compute_pseudo_inverse_factor(DX, rtol=0.25):
    """Truncated-SVD pseudo-inverse factor W with Binv_SVD = W W^T.

    Uses the thin SVD of DX. Returns W = sqrt(N-1) * U_k * diag(1/sigma_k)
    and the retained singular values sigma_k. Singular values with
    sigma/sigma_1 <= rtol are dropped to keep noise-dominated modes out.
    """
    n, N = DX.shape;
    U, s, _ = np.linalg.svd(DX, full_matrices=False);
    if s.size == 0 or s[0] == 0:
        return np.zeros((n, 0)), np.zeros(0);
    keep = s/s[0] > rtol;
    Uk = U[:, keep];
    sk = s[keep];
    W = np.sqrt(N-1.0) * Uk * (1.0/sk)[None, :];
    return W, sk;


####################################################################################
####################################################################################
####################################################################################
def compute_alpha_stein(Binv_mc, W, sk, N):
    """Stein-loss shrinkage weight between modified-Cholesky and SVD precisions.

    Closed-form approximation in the k-dimensional SVD subspace. Stein's loss
    is scale-invariant, so the modified-Cholesky target is used as-is.
    """
    k = W.shape[1];
    if k == 0:
        return 1.0;
    BW = Binv_mc @ W;
    if hasattr(BW, 'toarray'):
        BW = BW.toarray();
    BW = np.asarray(BW);
    WtBW = W.T @ BW;                                   # (k,k)
    UtBU = WtBW * np.outer(sk, sk) / (N-1.0);
    tr_Bmc_BsvdInv = float(np.trace(UtBU @ np.diag(sk**2) / (N-1.0)));
    sign, logdet_UtBU = np.linalg.slogdet(UtBU);
    if sign <= 0:
        return 1.0;
    logdet_Bsvd_sub = float(-2.0*np.log(sk).sum() + k*np.log(N-1.0));
    logdet_diff = logdet_UtBU - logdet_Bsvd_sub;
    inv_s2 = 1.0/(sk**2);
    tr_Bsvd2 = float(((N-1.0)**2) * (inv_s2**2).sum());
    tr_Bmc2 = _fro2_sparse(Binv_mc);
    tr_Bmc_Bsvd = float(np.einsum('ij,ij->', W, BW));
    denom = max(tr_Bmc2 + tr_Bsvd2 - 2.0*tr_Bmc_Bsvd, 1e-30);
    alpha = 0.5 + (tr_Bmc_BsvdInv - k - logdet_diff) / (2.0*denom);
    return float(min(1.0, max(0.0, alpha)));


####################################################################################
####################################################################################
####################################################################################
def compute_alpha_da(Binv_mc, W, HtRinvH, n, n_probes=12, n_iter=40):
    """DA-aware shrinkage weight by golden-section minimisation.

    Minimises J(alpha) = tr[(alpha*Binv_mc + (1-alpha)*W W^T + HtRinvH)^{-1}]
    on [0,1]. J is strictly convex; the trace of the inverse is estimated by a
    Hutchinson probe using the Sherman-Morrison-Woodbury identity, never
    forming an n x n inverse.
    """
    import scipy.sparse.linalg as sla;
    k = W.shape[1];
    rng = np.random.default_rng(0);
    Z = rng.choice([-1.0, 1.0], size=(n, n_probes));

    Bmc = Binv_mc.tocsc() if spa.issparse(Binv_mc) else np.asarray(Binv_mc);
    HtH = HtRinvH.tocsc() if spa.issparse(HtRinvH) else np.asarray(HtRinvH);

    def make_solver(A0):
        if spa.issparse(A0):
            lu = sla.splu(A0.tocsc());
            return lambda B: lu.solve(np.asarray(B));
        Achol = np.linalg.inv(np.asarray(A0));
        return lambda B: Achol @ np.asarray(B);

    def trA(alpha):
        A0 = alpha*Bmc + HtH;
        solve = make_solver(A0);
        if k > 0 and (1.0-alpha) > 0:
            A0invW = solve(W);
            cap = np.eye(k) + (1.0-alpha)*(W.T @ A0invW);
            cap_inv = np.linalg.inv(cap);
            A0invZ = solve(Z);
            corr = (1.0-alpha) * A0invW @ (cap_inv @ (W.T @ A0invZ));
            MinvZ = A0invZ - corr;
        else:
            MinvZ = solve(Z);
        return float(np.einsum('ij,ij->', Z, MinvZ) / n_probes);

    gr = (np.sqrt(5.0)-1.0)/2.0;
    a, b = 0.0, 1.0;
    c = b - gr*(b-a);
    d = a + gr*(b-a);
    fc, fd = trA(c), trA(d);
    for _ in range(n_iter):
        if fc < fd:
            b, d, fd = d, c, fc;
            c = b - gr*(b-a);
            fc = trA(c);
        else:
            a, c, fc = c, d, fd;
            d = a + gr*(b-a);
            fd = trA(d);
        if abs(b-a) < 1e-4:
            break;
    return float(min(1.0, max(0.0, 0.5*(a+b))));