import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as spa
import pandas as pd
import time
from netCDF4 import Dataset
np.random.seed(seed=10);

class reference_solution:
      nm = None;
      x_ref = None;
      M = None;

      def __init__(self, nm, M):
          self.nm = nm;
          self.M = M;
          self.load_reference_trajectory();
      
      def get_reference_state(self, k): 
          return self.x_ref[k];
          
      def load_reference_trajectory(self):
      
          M = self.M;
          nm = self.nm;
          self.x_ref = [];
          
          for s in range(0, M):
              nc_f = f'{nm.snapshots}reference_solution_{s}.nc';    
              xs = nm.load_netcdf_file(nc_f); 
              self.x_ref.append(xs);

