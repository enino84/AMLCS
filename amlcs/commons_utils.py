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
           D_sqrt[i] = 1/np.std(xi-pi @ beta_i);
           ind += pre_i.size;
        else:
           #D_sqrt[i] = -1;
           D_sqrt[i] = 1/np.std(xi);
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
def compute_coef_SVD(A,b,thr):
    N,n = A.shape;
    Ui,Si,Vi = np.linalg.svd(A,full_matrices=False);
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

    
