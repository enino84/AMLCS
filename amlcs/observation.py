import numpy as np
import scipy.sparse as spa
from netCDF4 import Dataset


####################################################################################
####################################################################################
####################################################################################
class observation:
      H_vect = None; #Observation operator as vector
      H_spar = None; #Observation operator as matrix sparse
      H_grid = None; #Observation operator as numerical grid
      R_invs = None;
      R_spar = None;
      R_vect = None;
      R_sqrt = None;
      
      obs_network = None;
      obs_H = None;
      obs_n = None;
      obs_m = None;
      obs_H_sparse = None;
      obs_R_sparse = None;
      
      m = -1; #Number of observations
      m_var = -1; #Number of observations for each variable.
      err_obs = None;
      y = None;
      y_obs = None;
      
      #s denotes the sparsity pattern per layer
      #nm - numerical model
      #gr - grid
      
      obs_var = None;
      
      def __init__(self,err_obs, obs_indexes):
          self.err_obs = err_obs;
          self.obs_var = obs_indexes;
          
          
 
      def get_stations_variables(self, lat, lon, p, s):
          S = np.zeros((lat, lon));
          for i in range(0,lat,s):
              for j in range(0,lon,s):
                  S[i,j] = 1;
          S_res = S.reshape(lon*lat, order='C');
          pos,=np.where(S_res>0);
          return list(pos+p);
          
          
      def build_observational_network(self, gr, nm, s):
          mesh = gr.mesh;
          mask = nm.mask_cor;
          self.obs_H = [];
          self.obs_n = [];
          self.obs_m = [];
          for mesh_, mask_ in zip(mesh, mask): #moving across meshes
              n_mesh = 0;
              m_mesh = 0;
              H_vect = [];
              o_data = [];
              p = 0;
              for v, n in zip(mesh_, mask_): #moving across variables
                  #exit();
                  mes_info = v[0];
                  res_info = v[1];
                  lat, lon = res_info[0], res_info[1];
                  #print(f'n[0] reads {n[0]}');
                  var_index = n[0][0]; #n[0] = [var_index, level]
                  if self.obs_var[var_index]: #If observations are available in this layer and this variable
                     H_l = self.get_stations_variables(lat, lon, n_mesh, s);
                     m_mesh = len(H_l);
                     H_vect.extend(H_l);
                     o_data.append([n[0], m_mesh]);
                  n_mesh+= lat * lon; #we skip indexes from being observed at this layer-var
              self.obs_H.append(np.array(H_vect));
              self.obs_n.append(n_mesh);
              self.obs_m.append(o_data);
          #print('* self.obs_H {0}'.format(self.obs_H));
          #print('* self.obs_n {0}'.format(self.obs_n));
          #print('* self.obs_m {0}'.format(self.obs_m));
          
          self.build_observational_operator();
          self.build_data_error_covariance(gr, nm);
      
     
      def build_observational_operator(self):
          obs_H = self.obs_H;
          obs_n = self.obs_n;
          n_sub = len(obs_H);
          self.obs_H_sparse = [];
          for s in range(0, n_sub):
              n = obs_n[s];
              J = obs_H[s]; #indexes wherein observations are located
              V = np.ones_like(J);
              m = J.size;
              I = np.arange(0, m);
              H_sparse = spa.coo_matrix((V,(I,J)), shape=(m, n));
              self.obs_H_sparse.append(H_sparse);
      
      def build_data_error_covariance(self, gr, nm):
          obs_m = self.obs_m;
          self.obs_R_sparse = [];
          for om in obs_m:
              #print('* om es {0}'.format(om)); #[[var, level], number of observations]
              m_om = 0;
              Ig = [];
              V_R_g = [];
              V_Ri_g = [];
              V_Rs_g = [];
              for o in om:
                  var_info = o[0]
                  variable = var_info[0];
                  err_obm = self.err_obs[variable];
                  m = o[1];
                  I = np.arange(m_om, m_om+m);
                  Ig.extend(list(I));
                  V_R_g.extend(list((err_obm**2)*np.ones_like(I)))
                  V_Ri_g.extend(list((1/(err_obm**2))*np.ones_like(I)))
                  V_Rs_g.extend(list(err_obm*np.ones_like(I)))
                  m_om+=m;
               
              Ig = np.array(Ig);
              V_R_g = np.array(V_R_g); 
              V_Ri_g = np.array(V_Ri_g);  
              V_Rs_g = np.array(V_Rs_g);  
              
              #print([Ig.shape, V_R_g.shape]);   

              
              R_sparse = spa.coo_matrix((V_R_g,(Ig,Ig)), shape=(m_om, m_om));
              Rinv_sparse = spa.coo_matrix((V_Ri_g,(Ig,Ig)), shape=(m_om, m_om));
              Rsqr_sparse = spa.coo_matrix((V_Rs_g,(Ig,Ig)), shape=(m_om, m_om));
              
              self.obs_R_sparse.append({'R':R_sparse, 'Ri':Rinv_sparse, 'Rs':Rsqr_sparse, 'm':m_om});
           
 
          
          
      
      def build_synthetic_observations(self, nm, rs, M):
          reference_abs = nm.path+'reference/';
          reference_path = reference_abs+'snapshots/';
          mask_cor = nm.mask_cor;
          var_names = nm.var_names;
          self.y_obs = [];
          for s in range(0, M):
              xs = rs.x_ref[s]; 
              y_ma = [];
              for ma, R_data, H_sparse in zip(mask_cor, self.obs_R_sparse, self.obs_H_sparse):
                  x_data = [];
                  for m in ma:
                      var_info = m[0];
                      variable = var_info[0];
                      level = var_info[1];
                      var_name = var_names[variable];
                      #print([var_name, m[1]]);
                      if 'TRG' in var_name:
                         x_ma = xs[variable][level,:,:].reshape((-1,1));
                      elif 'PSG' in var_name:
                         x_ma = xs[variable][:,:].reshape((-1,1));
                      else:
                         x_ma = xs[variable][level,:,:].reshape((-1,1));
                      x_data.extend(list(x_ma));
                      #print(x_ma);
                      
                  x_data = np.array(x_data);
                  R_sqrt = R_data['Rs'];
                  m_obs = R_sqrt.size;
                  y = H_sparse @ x_data + R_sqrt @ np.random.randn(m_obs,1);
                  y_ma.append(y);
              
              self.y_obs.append(y_ma); 
          
          print('* ENDJ - Synthetic observations have been created');
          
  
      def get_perturbed_observations(self, block, k, N):
          np.random.seed(seed=10+k); #To replicate observations
          y_k = self.y_obs[k];
          y_block = y_k[block];
          R_sqrt_block = self.obs_R_sparse[block]['Rs'];
          H_block = self.obs_H_sparse[block];
          m_block = self.obs_R_sparse[block]['m'];
          Ys_block = y_block + R_sqrt_block @ np.random.randn(m_block, N);
          return Ys_block;
      
      def get_R(self, block):
          return self.obs_R_sparse[block]['R'];

      def get_H(self, block):
          return self.obs_H_sparse[block];

      def get_synthetic_noise(self,N):
          return self.R_sqrt @ np.random.randn(self.m,N);
