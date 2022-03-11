####################################################################################
####################################################################################
####################################################################################

import os
import numpy as np
import sys
import scipy.sparse as spa
import pandas as pd
import time
from netCDF4 import Dataset
np.random.seed(seed=10);

from numerical_model import numerical_model
from grid_resolution import grid_resolution
            

####################################################################################
####################################################################################
####################################################################################
def main():
    
    
    input_file = sys.argv[1];
    print('* ENDJ - Input file {0}'.format(input_file));
    df_par = pd.read_csv(input_file);
    for c in df_par.columns:
        print(' - {0} = {1}'.format(c, df_par[c].iloc[0]));

        
    #1. Prepaging
    Nens      = df_par['Nens'].iloc[0];
    M         = df_par['M'].iloc[0];
    res_name  = df_par['res_name'].iloc[0];
    per       = df_par['per'].iloc[0];
    obs_steps = df_par['obs_steps'].iloc[0];
    ini_steps = df_par['ini_steps'].iloc[0];
    ini_times = df_par['ini_times'].iloc[0];
    syn_tests = df_par['syn_tests'].iloc[0];
    data_prep = df_par['folder_prep'].iloc[0];
    code_prep = df_par['code'].iloc[0];
    par = df_par['par'].iloc[0].astype(bool);
    
    args = [0, obs_steps, 1];
    ini0 = [ini_steps, 0, 1];
    ini0_no_restart = [ini_times, 0, 0];

    if np.isnan(code_prep):
       path = data_prep+'/'+res_name+'_'+str(Nens)+'_'+str(per)+'_'+str(M)+'/';
       code_path = res_name+'_'+str(Nens)+'_'+str(per)+'_'+str(M);
    else:
       path = data_prep+'/'+df_par['code'].iloc[0]+'/';
       code_path = df_par['code'].iloc[0];
    
    os.system('mkdir '+data_prep+'/');
    os.system('mkdir '+path);
    
    gr = grid_resolution(res_name);

    nm = numerical_model(path, gr, Nens, par=par);
   
    #Creating the initial condition
    if syn_tests:
       nm.create_initial_condition(ini0_no_restart);   
       nm.create_reference_snapshots(ini0, args, M);
       nm.create_initial_ensemble(ini0, args, per, Nens);
       nm.create_free_run(ini0, args, M);
       nm.collect_ensemble_members();
    else:
       pass; #real data
    
    df_par['path'] = path;
    df_par['code_path'] = code_path;
    df_par.to_csv(f'{path}/config.csv', index=False);
    exit();
    
if __name__ == "__main__":
    main();
