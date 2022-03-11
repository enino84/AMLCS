import os
import numpy as np
import sys
import scipy.sparse as spa
import pandas as pd
import time
from netCDF4 import Dataset
from commons_utils import compute_modified_Cholesky_decomposition


##########################################################################################
##########################################################################################
##########################################################################################
# General class - sequential ensemble data assimilation
##########################################################################################
##########################################################################################
##########################################################################################
class ensemble_DA:

    XB = None; 
    XA = None;
    XB_map = None;
    XA_map = None;
    Nens = None;
    nm = None;
    gr = None;
    XB_s = [];
    Ys = None;
    infla = None;
    
    def __init__(self, nm, infla, Nens):
        self.Nens = Nens;
        self.nm = nm;
        self.infla = infla;
        self.test = 1;
        #self.load_background_ensemble();
        
    

    def load_background_ensemble(self):
        nm = self.nm;
        nm.load_ensemble();
        self.XB = nm.X.copy();
    
    
    def get_ensemble_variable(self, v):
        Nens = self.Nens;
        v_info = v[0];
        v_reso = v[1];
        var_index = v_info[0];
        var_level = v_info[1];
        lat, lon = v_reso[0], v_reso[1];
        n = lat*lon;
        XB_v = np.zeros((n,Nens));
        nm = self.nm;
        for e in range(0, Nens):
            if 'PSG' in nm.var_names[var_index]:
                xb_e_v = self.XB[e][var_index][:,:].reshape((n,));
            else:
                xb_e_v = self.XB[e][var_index][var_level,:,:].reshape((n,));
            XB_v[:,e] = xb_e_v;
        return XB_v;
    
    
    def get_ensemble_block(self, msk_cor):
        n_vars = len(msk_cor);
        XB_block = np.empty(shape=(0,self.Nens));
        if self.test == 1: print('* ENDJ - Working on block {0}'.format(msk_cor));
        for v in msk_cor:
            XB_v = self.get_ensemble_variable(v);
            XB_block = np.concatenate((XB_block, XB_v), axis=0);
        return XB_block;
            
    
    def prepare_background(self):
        pass;
         
    def prepare_analysis(self, ob, k, args = None):
        pass;
    
    def perform_assimilation(self, ob):
        pass;
    
    def map_vector_states(self):
        self.XA = [];
        Nens = self.Nens;
        for e in range(0, Nens):
            xa_e = self.nm.map_vector_state(self.XA_map, e);
            nc_f = f'{self.nm.ensemble_0}ens_{e}/ensemble_member.nc';
            self.nm.map_state_netcdf(xa_e, nc_f); 
            self.XA.append(xa_e);
    
    def perform_forecast(self):
        self.nm.forecast_ensemble(self.Nens);
        
    def check_time_store(self, k, list_k):
        if k in list_k:
           nm = self.nm;
           xb_k = nm.compute_snapshot_mean(self.XB);
           fn_xb = f'{nm.time_snapshots}xb{k}.nc';
           xa_k = nm.compute_snapshot_mean(self.XA);
           fn_xa = f'{nm.time_snapshots}xa{k}.nc';
           nm.map_state_netcdf(xb_k, fn_xb);
           nm.map_state_netcdf(xa_k, fn_xa);
           
    
    def covariance_inflation(self, XA):  
        xa = XA.mean(axis=1).reshape(-1,1);
        DX = XA-xa;
        XA = xa + self.infla * DX;  
        return XA;
        
    def clear_all(self):
        self.nm.clear_all_folders();

##########################################################################################
##########################################################################################
##########################################################################################
# Observation space version of:
# Nino-Ruiz, E. D., Sandu, A., & Deng, X. (2018). An ensemble Kalman filter implementation based on modified Cholesky decomposition for inverse covariance matrix estimation. SIAM Journal on Scientific Computing, 40(2), A867-A886.
# To be published
##########################################################################################
##########################################################################################
##########################################################################################
class EnKF_MC_obs(ensemble_DA):
            
    
    def prepare_background(self):
        mask_cor = self.nm.mask_cor;
        var_names = self.nm.var_names;
        Nens = self.Nens;
        self.XB_map = [];
        gr = self.nm.gs;
        for msk_cor, pre_info in zip(mask_cor, gr.lpr):
            XB_block = self.get_ensemble_block(msk_cor);
            xb_block = XB_block.mean(axis=1).reshape(-1,1);
            DX_block = XB_block-xb_block;
            Binv_sqrt_block = compute_modified_Cholesky_decomposition(DX_block, pre_info, thr=0.15);

            self.XB_map.append({'XB_b':XB_block, 'xb_b':xb_block, 'Binv_s_b':Binv_sqrt_block});
     
    def prepare_analysis(self, ob, k, args = None):
        self.Ys = [];
        Nens = self.Nens;
        mask_cor = self.nm.mask_cor;
        N_blocks = len(mask_cor);
        for block in range(0, N_blocks):
            Ys_block = ob.get_perturbed_observations(block, k, Nens);  
            self.Ys.append(Ys_block); 
    
    def perform_assimilation(self, ob):
        self.XA_map = [];
        for XB_info, Ys_block, R_info, H_block in zip(self.XB_map, self.Ys, ob.obs_R_sparse, ob.obs_H_sparse):
            XB_block = XB_info['XB_b'];
            R_block  = R_info['R'];
            Binv_sqrt_block = XB_info['Binv_s_b'];
            XA_block = self.perform_assimilation_block(XB_block, Binv_sqrt_block, H_block, R_block, Ys_block);
            XA_block = self.covariance_inflation(XA_block);
            self.XA_map.append(XA_block);
        self.map_vector_states(); #Update ensemble folders
    
    def perform_assimilation_block(self, XB, Binv_sqrt, H, R, Ys):
      
          Ds = Ys - H @ XB;
           
          H_spar = H.toarray();
          
          P = spa.linalg.spsolve_triangular(Binv_sqrt, H_spar.T, lower=False, check_finite=False);
          
          Inno = R + P.T @ P;

          Q_temp = P @ spa.linalg.spsolve(Inno, Ds);
          
          #Q_temp = spa.linalg.spsolve_triangular(Binv_sqrt,Q_temp,lower=False);
          
          DXa = spa.linalg.spsolve_triangular(Binv_sqrt.T, Q_temp, lower=True, check_finite=False);

          XA = XB + DXa;  
          
          return XA;
      
         
##########################################################################################
##########################################################################################
##########################################################################################
# Miyoshi, T., & Yamane, S. (2007). Local ensemble transform Kalman filtering with an AGCM at a T159/L48 resolution. Monthly Weather Review, 135(11), 3841-3861.
##########################################################################################
##########################################################################################
##########################################################################################
class LETKF(ensemble_DA):
            
    
    def prepare_background(self):
        mask_cor = self.nm.mask_cor;
        var_names = self.nm.var_names;
        Nens = self.Nens;
        self.XB_map = [];
        gr = self.nm.gs;
        for msk_cor, lbo_info in zip(mask_cor, gr.lbo):
            XB_block = self.get_ensemble_block(msk_cor);
            xb_block = XB_block.mean(axis=1).reshape(-1,1);
            DX_block = XB_block-xb_block;
            self.XB_map.append({'XB_b':XB_block, 'xb_b':xb_block, 'DX_b':DX_block, 'lbo_b':lbo_info});


    def prepare_analysis(self, ob, k, args = None):
        self.y = [];
        Nens = self.Nens;
        mask_cor = self.nm.mask_cor;
        N_blocks = len(mask_cor);
        for block in range(0, N_blocks):
            self.y.append(ob.y_obs[k][block]); 
        self.nm.gs.compute_local_boxobs(ob.obs_H_sparse);
    
    def perform_assimilation(self, ob):
        self.XA_map = [];
        for XB_info, y_block, R_info, H_block, lobs_block in zip(self.XB_map, self.y, ob.obs_R_sparse, ob.obs_H_sparse, self.nm.gs.lbo_obs):
            XB_block = XB_info['XB_b'];
            xb_block = XB_info['xb_b'];
            DX_block = XB_info['DX_b'];
            lbo_block = XB_info['lbo_b'];
            Ri_block  = R_info['Ri'];
            XA_block = self.perform_assimilation_block(XB_block, xb_block, DX_block, H_block, Ri_block, y_block, lbo_block, lobs_block);
            self.XA_map.append(XA_block);
        self.map_vector_states(); #Update ensemble folders
        
    
    def perform_assimilation_local_box(self, XB, xb, DX, H, Ri, y):
        
        n, Nens = XB.shape;
        d = y - H @ xb;
        Yb = H @ DX;
        
        #print(f'Ri {Ri.shape}');
        #print(f'Yb {Yb.shape}');
        #print(f'd {d.shape}');
        #print(f'H {H.shape}');
        #print(f'xb {xb.shape}');
        #print(f'y {y.shape}');
        
        Pa_Nens = (Nens-1)*np.eye(Nens) + Yb.T @ ( Ri @ Yb );
        Q_temp = Yb.T @ (Ri @ d);
        
        U, S, _ = np.linalg.svd(Pa_Nens, full_matrices=False);
        
        Pa_sqrt = U @ ( np.diag(np.sqrt(Nens/S)) @ U.T );
        Pa_invs = U @ ( np.diag(1/S) @ U.T );
        
        wa = Pa_invs @ Q_temp;
        
        xa = xb + DX @ wa;
        
        XA = xa + DX @ Pa_sqrt;
        
        return XA;
        
        
    
    def perform_assimilation_block(self, XB, xb, DX, H, Ri, y, lbo_info, lobs):

          
          n, Nens = XB.shape;
          
          lbo, nlbo = lbo_info; #local boxes information (indexes, number of components in all boxes)

          y_model = H.T @ y;
          
          
          Ri_space = H.T @ (Ri @ H);
          Ri_space = Ri_space.toarray();
          
          XA = np.zeros((n, Nens));
          
          
          for i in range(0, n):
              #local box for model component i
              lbo_i = np.array(lbo[i]).astype('int32');
              gp_i, = np.where(lbo_i==i);
              xb_i = xb[lbo_i];
              XB_i = XB[lbo_i];
              #local observation operator
              H_ind = np.array(lobs[i]).astype('int32'); #local observed components
              #print(f'H_ind {H_ind}');
              m_i = H_ind.size;
              if m_i>0:
                 n_i = xb_i.size;
                 I = np.arange(0, m_i);
                 J = H_ind;
                 H_i = spa.coo_matrix((np.ones(m_i),(I,J)), shape=(m_i, n_i));
                 DX_i = DX[lbo_i];
                 y_i  = y_model[lbo_i[H_ind]];
                 Ri_i = np.diag(Ri_space[lbo_i[H_ind], lbo_i[H_ind]]).reshape((m_i, m_i)); #local data error covariance matrix
                 XA_i = self.perform_assimilation_local_box(XB_i, xb_i, DX_i, H_i, Ri_i, y_i);
                 #perform_assimilation_local_box(self, XB, xb, DX, H, Ri, y)
              else:
                 XA_i = XB_i;
                 
              XA[i, :] = XA_i[gp_i, :]; 
          
          return XA;
          #XB, xb, DX, H, Ri, y):
          



##########################################################################################
##########################################################################################
##########################################################################################
# Ott, E., Hunt, B. R., Szunyogh, I., Zimin, A. V., Kostelich, E. J., Corazza, M., ... & Yorke, J. A. (2004). A local ensemble Kalman filter for atmospheric data assimilation. Tellus A: Dynamic Meteorology and Oceanography, 56(5), 415-428.
##########################################################################################
##########################################################################################
##########################################################################################
class LEnKF(ensemble_DA):
            
    
    def prepare_background(self):
        mask_cor = self.nm.mask_cor;
        var_names = self.nm.var_names;
        Nens = self.Nens;
        self.XB_map = [];
        gr = self.nm.gs;
        for msk_cor, lbo_info in zip(mask_cor, gr.lbo):
            XB_block = self.get_ensemble_block(msk_cor);
            xb_block = XB_block.mean(axis=1).reshape(-1,1);
            DX_block = XB_block-xb_block;
            self.XB_map.append({'XB_b':XB_block, 'DX_b':DX_block, 'lbo_b':lbo_info});


    def prepare_analysis(self, ob, k, args = None):
        self.Ys = [];
        Nens = self.Nens;
        mask_cor = self.nm.mask_cor;
        N_blocks = len(mask_cor);
        for block in range(0, N_blocks):
            Ys_block = ob.get_perturbed_observations(block, k, Nens);  
            self.Ys.append(Ys_block);
        self.nm.gs.compute_local_boxobs(ob.obs_H_sparse); 
    
    def perform_assimilation(self, ob):
        self.XA_map = [];
        for XB_info, Ys_block, R_info, H_block, lobs_block in zip(self.XB_map, self.Ys, ob.obs_R_sparse, ob.obs_H_sparse, self.nm.gs.lbo_obs):
            XB_block = XB_info['XB_b'];
            lbo_block = XB_info['lbo_b'];
            Ri_block  = R_info['Ri'];
            XA_block = self.perform_assimilation_block(XB_block, H_block, Ri_block, Ys_block, lbo_block, lobs_block);
            self.XA_map.append(XA_block);
        self.map_vector_states(); #Update ensemble folders
        
    
    def perform_assimilation_local_box(self, XB, H, R, Ys):
        
        Ds = Ys - H @ XB;
        
        Pb = np.cov(XB);
        
        XA = XB + Pb @ ( H.T @ np.linalg.solve(R + H @ (Pb @ H.T), Ds) );
        
        return XA;
        
        
    
    def perform_assimilation_block(self, XB, H, Ri, Ys, lbo_info, lobs):

          
          n, Nens = XB.shape;
          
          lbo, nlbo = lbo_info; #local boxes information (indexes, number of components in all boxes)

          Ys_model = H.T @ Ys;
          
          
          Ri_space = H.T @ (Ri @ H);
          Ri_space = Ri_space.toarray();
          
          XA = np.zeros((n, Nens));
          
          
          for i in range(0, n):
              #local box for model component i
              lbo_i = np.array(lbo[i]).astype('int32');
              gp_i, = np.where(lbo_i==i);
              XB_i = XB[lbo_i];
              #local observation operator
              H_ind = np.array(lobs[i]).astype('int32'); #local observed components
              #print(f'H_ind {H_ind}');
              m_i = H_ind.size;
              if m_i>0:
                 n_i,_ = XB_i.shape;
                 I = np.arange(0, m_i);
                 J = H_ind;
                 H_i = spa.coo_matrix((np.ones(m_i),(I,J)), shape=(m_i, n_i));
                 Ys_i  = Ys_model[lbo_i[H_ind]];
                 Ri_i = np.diag(Ri_space[lbo_i[H_ind], lbo_i[H_ind]]).reshape((m_i, m_i)); #local data error covariance matrix
                 XA_i = self.perform_assimilation_local_box(XB_i, H_i, Ri_i, Ys_i);
                 #perform_assimilation_local_box(self, XB, xb, DX, H, Ri, y)
              else:
                 XA_i = XB_i;
                 
              XA[i, :] = XA_i[gp_i, :]; 
          
          return XA;
          #XB, xb, DX, H, Ri, y):   
##########################################################################################
##########################################################################################
##########################################################################################
# General factory - sequential ensemble data assimilation
##########################################################################################
##########################################################################################
##########################################################################################
class sequential_method:
      
      method_name = None;
      
      def __init__(self, method_name):
          self.method_name = method_name;
      
      def get_instance(self, nm, infla, Nens):
          if self.method_name=='EnKF_MC_obs': return EnKF_MC_obs(nm, infla, Nens);
          if self.method_name=='LETKF': return LETKF(nm, infla, Nens);
          if self.method_name=='LEnKF': return LEnKF(nm, infla, Nens);
          



    
    
    
    
        
    
