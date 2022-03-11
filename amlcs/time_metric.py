import os
import numpy as np
import sys
import pandas as pd
import time
from netCDF4 import Dataset
np.random.seed(seed=10);

class time_metric:
      
      ini_time = None;
      end_time = None;
      time_col = None;
      
      def __init__(self, nm, kind):
          self.nm = nm;
          self.path = self.nm.path+'results/';
          self.kind = kind;
          self.time_col = [];
       
      def start_time(self):
          self.ini_time = time.time();

      def check_time(self):
          self.end_time = time.time();
          self.time_col.append(self.end_time - self.ini_time);
          self.ini_time = None;
          self.end_time = None;

      def store_all_results(self):
          time_col = np.array(self.time_col);
          pd_time = pd.DataFrame(time_col.reshape(-1,1));
          pd_time.columns = ['time'];
          pd_time.to_csv(f'{self.path}time_{self.kind}.csv');
             
                  
          
              
          
      
      
