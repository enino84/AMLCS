import os
import numpy as np
import pandas as pd

class error_metric:
      
      nm = None;
      M = None;
      error_steps = None;
      path = None;
      kind = None;
      
      def __init__(self, nm, kind, M):
          self.nm = nm;
          self.path = self.nm.path+'results/';
          self.M = M;
          self.kind = kind;
          self.create_error_steps(M);
       
      def create_error_steps(self, M):
          self.error_steps = [];
          nvars = len(self.nm.var_names);
          for v in range(0, nvars):
              if 'PSG' in self.nm.var_names[v]:
                 self.error_steps.append(np.zeros((1,M)));
              else:
                 self.error_steps.append(np.zeros((8,M)))
          os.system('mkdir '+self.path);
          
      def compute_error_step(self, k, X, xr):
          
          xm = self.nm.compute_snapshot_mean(X);
          var_names = self.nm.var_names;
          var_resol = self.nm.var_resol;
          n_vars = len(var_names);
          for v in range(0, n_vars):
              var_data = self.error_steps[v];
              #print(f'* var_data {var_data}');
              if 'PSG' in var_names[v]:
                  err_lev = np.linalg.norm(xm[v][:,:].reshape(-1,)-xr[v][:,:].reshape(-1,1));
                  var_data[0,k] = err_lev;
              else:
                  
                  lat, lon = var_resol[v][0];
                  lev = var_resol[v][1];
                  #print('* lev reads {0}'.format(var_resol));
                  err_lev = np.zeros((lev,));
                  for l in range(0, lev):
                      err_l = np.linalg.norm(xm[v][l,:,:].reshape(-1,1)-xr[v][l,:,:].reshape(-1,1));
                      err_lev[l] = err_l;
                  var_data[:,k] = err_lev;
              self.error_steps[v] = var_data;
              
      def get_rmse_per_variable_level(self, var, lev):
          var_data = self.error_steps[var];
          if 'PSG' in self.nm.var_names[var]:
              rmse_var = np.sqrt(np.square(var_data[0, :]).mean());
          else:
              rmse_var = np.sqrt(np.square(var_data[lev, :].mean()));
          return rmse_var;
          
          
      def get_error_per_variable_level(self, var, lev, k):
          var_data = self.error_steps[var];
          if 'PSG' in self.nm.var_names[var]:
              rmse_var = var_data[0, k];
          else:
              rmse_var = var_data[lev, k];
          return rmse_var;
      

      def store_all_results(self):
          nvars = len(self.nm.var_names);
          for v in range(0, nvars):
              err_v = pd.DataFrame(self.error_steps[v].T);
              err_v.to_csv(self.path+self.nm.var_names[v]+'_'+self.kind+'.csv');
             
                  
          
              
          
      
      
