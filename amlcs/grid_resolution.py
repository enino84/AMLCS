import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as spa
from sklearn.utils.extmath import randomized_svd
import pandas as pd
import time
from sklearn.linear_model import Ridge, Lasso
from netCDF4 import Dataset
np.random.seed(seed=10);

####################################################################################
####################################################################################
####################################################################################
class grid_resolution:
      name = None;
      lon = None;
      lat = None;
      lev = None;
      sub_doms = None;
      local_box = [];
      local_pre = [];
      local_obs = [];
      local_pre_obs = [];
      npr = -1;
      mpr = -1;
      L = None;
      lbo = None;
      lpr = None;
      
      def __init__(self,name):
          self.name = name;
          self.set_resolution(name);
          
      def compute_local_box(self, i, j, me, v1, r):
          loc_ind = [];
          mes_info = me[v1][0];
          res_info = me[v1][1];
          lat, lon = res_info[0], res_info[1];
          for i1 in range(i-r, i+r+1): 
              if i1>=0 and i1<lat:
                 for j1 in range(j-r, j+r+1):
                     jlon = j1;
                     if j1>=lon: jlon = j1-lon;
                     loc_ind.append(mes_info[i1, jlon]);
          lbo = np.array(loc_ind);
          return lbo
          
      
      def compute_local_boxobs(self, H):
          self.lbo_obs = [];
          for lbo_block, Hbo in zip(self.lbo, H):
              lbo_o = [];
              m, n = Hbo.shape;
              isobs = Hbo.T @ np.ones((m,)); #1 there is an observation 0 there is no observation
              lbo, nlblo = lbo_block;
              for lb_i in lbo:
                  lb = np.array(lb_i);
                  #print('lb {0}'.format(lb));
                  lb = lb.astype('int32');
                  ind_obs, = np.where(isobs[lb_i]>0);
                  ind_obs = ind_obs.astype('int32'); #local indexes
                  #print(f'* ind_obs = {ind_obs} isobs = {isobs}');
                  obs_i = ind_obs; #lb[ind_obs];
                  #elias
                  lbo_o.append(obs_i);
              self.lbo_obs.append(lbo_o);
                  
      
      def compute_local_boxes(self, me, r):
          n_vars = len(me);
          npre = 0;
          nlbo = 0;
          lbo = [];
          pre = [];
          for v in range(0, n_vars):
              mes_info = me[v][0];
              res_info = me[v][1];
              latv, lonv = res_info[0], res_info[1];
              for i in range(0, latv):
                  for j in range(0, lonv):
                      lab = mes_info[i, j];
                      prev = [];
                      lbov = [];
                      for v1 in range(0, n_vars):
                          lbo1 = self.compute_local_box(i, j, me, v1, r);
                          pre1 = lbo1[np.where(lbo1<lab)];
                          prev.extend(pre1);
                          lbov.extend(lbo1);
                          npre+=pre1.size;
                          nlbo+=lbo1.size;
                      pre.append(prev);
                      lbo.append(lbov);
          return (lbo, nlbo), (pre, npre);
              
              
      
      def compute_sub_domains(self, r):
          self.lbo = [];
          self.lpr = [];
          for me in self.mesh:
              lbo_info, pre_info = self.compute_local_boxes(me, r);
              self.lbo.append(lbo_info);
              self.lpr.append(pre_info);
              
                                 
 
      def get_labels(self, p, lat, lon):
          n_lonlat = lon*lat;
          M = np.arange(0,n_lonlat).reshape((lat, lon), order='C');
          return M + p;
      
      def get_subset_labels(self, msk_cor):
          n_vars = len(msk_cor);
          L = [];
          #print(msk_cor);
          p = 0;
          for vs in msk_cor:
              var_info = vs[0];
              res_info = vs[1];
              #print('vs {0} var_info {1} res_info {2}'.format(vs, var_info, res_info));
              lat, lon = res_info[0], res_info[1];
              L.append((self.get_labels(p, lat, lon),(lat, lon)));
              p+=lon*lat;
          return L;    



                 
      def create_mesh(self, nm):
          mask_cor = nm.mask_cor;
          self.mesh = [];
          for msk_cor in mask_cor:
              vs_L = self.get_subset_labels(msk_cor);
              #print('* vs {0} and {1}'.format(vs, vs_L));
              self.mesh.append(vs_L);
          print('* ENDJ - Numerical mesh ready');
      
      
      def set_resolution(self,res_name):
          self.lev = 8;
          if res_name=='t21': 
             self.lat = 32; 
             self.lon = 64;
          elif res_name=='t30': 
             self.lat = 48; 
             self.lon = 96;
          elif res_name=='t47': 
             self.lat = 72; 
             self.lon = 144;
          elif res_name=='t63': 
             self.lat = 96; 
             self.lon = 192;
          elif res_name=='t106': 
             self.lat = 160; 
             self.lon = 320;
          else: 
             print('* Invalid resolution '+self.res); 
             exit();
      
      def get_resolution(self,res_name):
          if res_name=='t21': 
             lat = 32; 
             lon = 64;
          elif res_name=='t30': 
             lat = 48; 
             lon = 96;
          elif res_name=='t47': 
             lat = 72; 
             lon = 144;
          elif res_name=='t63': 
             lat = 96; 
             lon = 192;
          elif res_name=='t106': 
             lat = 160; 
             lon = 320;
          else: 
             print('* Invalid resolution '+self.res); 
             exit();
          return lat, lon;
              
          
          
          #ls = np.array(local_obs);
          #pd.DataFrame(ls).to_csv('local_obs.csv');

          #ls = np.array(local_pre_obs);
          #pd.DataFrame(ls).to_csv('local_pre_obs.csv');
          
          #exit();
          #return local_box,local_pre,local_obs,local_pre_obs,mpr,npr,p,n,H_relp;    

