####################################################################################
####################################################################################
####################################################################################

import os
import sys
import numpy as np
import scipy.sparse as spa
import pandas as pd
import time
from netCDF4 import Dataset

from numerical_model import numerical_model
from grid_resolution import grid_resolution
from observation import observation
from reference_solution import reference_solution
from sequential_methods import EnKF_MC_obs, sequential_method
from error_metric import error_metric
from time_metric import time_metric
            

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
    r            = int(df_par['r'].iloc[0]);
    s            = int(df_par['s'].iloc[0]);
    method       = df_par['method'].iloc[0].strip();	
    exp_settings = df_par['exp_settings'].iloc[0];
    infla        = float(df_par['infla'].iloc[0]);  
    data_obs     = df_par['err_obs'].iloc[0].strip().split(',');
    err_obs      = [float(v) for v in data_obs];   
    plac_obs     = df_par['obs_plc'].iloc[0].strip().split(',');
    obs_plc      = [int(v) for v in plac_obs];
    l_snap       = df_par['list_snapshots'].iloc[0].strip().split(',');
    option_mask  = int(df_par['option_mask'].iloc[0]);
    
    list_k = [int(v) for v in l_snap];

    print(f'plac_obs = {obs_plc}');
    print(f'err_obs  = {err_obs}');
    
    df_con = pd.read_csv(f'{exp_settings}/config.csv');
    print(df_con);
    
    Nens      = df_con['Nens'].iloc[0];
    M         = 1;#df_con['M'].iloc[0];
    res_name  = df_con['res_name'].iloc[0];
    per       = df_con['per'].iloc[0];
    obs_steps = df_con['obs_steps'].iloc[0];
    ini_steps = df_con['ini_steps'].iloc[0];
    ini_times = df_con['ini_times'].iloc[0];
    syn_tests = df_con['syn_tests'].iloc[0];
    data_prep = df_con['folder_prep'].iloc[0];
    code_prep = df_con['code'].iloc[0];
    code_path = df_con['code_path'].iloc[0];
    par       = df_con['par'].iloc[0].astype(bool);
    
    args = [0, obs_steps, 1];
    ini0 = [ini_steps, 0, 1];
    ini0_no_restart = [ini_times, 0, 0];
    #print(df_par['code'])
     
    method_path = code_path+'_'+method+'_'+str(r)+'_'+str(s)+'_'+str(int(100*infla));
    
    path = '../test/';

    os.system(f'mkdir {path}');
    os.system(f'mkdir {path}{method_path}');

    path_method = f'{path}{method_path}/';
    
    print('* The method reads {0}'.format(method));
    
    #exit();
    
    #1.2 Grid resolution
    gs = grid_resolution(res_name);
    data_obs = df_par['err_obs'].iloc[0].strip().split(',');
    err_obs = [float(v) for v in data_obs];#[0.1, 0.1, 0.1, 1e-4, 0.01, 0.1, 0.1, 0.1, 1e-4, 0.01]; #[u,v,T,H,p]
    #obs_plc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
    #print (err_obs)
    
    plac_obs = df_par['obs_plc'].iloc[0].strip().split(',');
    obs_plc = [int(v) for v in plac_obs];
    #print (obs_plc)
    ob = observation(err_obs, obs_plc);

    
    #1.1 Numerical model
    nm = numerical_model(path_method, gs, Nens, par=par);
    nm.define_relations(option=100);
    nm.load_settings(exp_settings, args, no_test=False);
    #print(nm.mask_cor);
     
    
    #1.3 Mesh grid is created
    gs.create_mesh(nm);
    gs.compute_sub_domains(r);
    rs = reference_solution(nm, M);
    
    
    #Setting the sequential data assimilation method
    seq_da = sequential_method(method).get_instance(nm, infla, Nens);
    
    ob.build_observational_network(gs, nm, s=s);
    ob.build_synthetic_observations(nm, rs, M); 
    
    #Metric
    em_bck = error_metric(nm, 'bck', M);
    em_ana = error_metric(nm, 'ana', M);
    tm_bck =  time_metric(nm, 'bck');
    tm_ana =  time_metric(nm, 'ana');
    
    seq_da.load_background_ensemble();
    tm_bck.start_time();
    seq_da.prepare_background();
    tm_bck.check_time();
    
    k = 0;
    
    tm_ana.start_time();
    seq_da.prepare_analysis(ob, k);
    seq_da.perform_assimilation(ob);
    tm_ana.check_time();
    
    em_bck.compute_error_step(k, seq_da.XB, rs.x_ref[k]);
    em_ana.compute_error_step(k, seq_da.XA, rs.x_ref[k]);
    
    em_bck.store_all_results();
    em_ana.store_all_results();
    tm_bck.store_all_results();
    tm_ana.store_all_results();
        
    seq_da.check_time_store(k, list_k);
        
    
    err_bac = em_bck.get_error_per_variable_level(0, 0, 0);
    err_ana = em_ana.get_error_per_variable_level(0, 0, 0);
    
    print('* ENDJ - background error = {0} analysis error = {1}'.format(err_bac, err_ana));
    
    
    err_bac = em_bck.get_error_per_variable_level(1, 0, 0);
    err_ana = em_ana.get_error_per_variable_level(1, 0, 0);
    
    
    print('* ENDJ - background error = {0} analysis error = {1}'.format(err_bac, err_ana));
    
    
    seq_da.clear_all();
    
    exit();
    
if __name__ == "__main__":
    main();
