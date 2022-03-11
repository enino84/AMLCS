import os
import numpy as np
from netCDF4 import Dataset
from grid_resolution import grid_resolution
import multiprocessing as mp


class numerical_model:
      name = 'speedy';
      
      path = None;
      res = None;
      source_model = None;
      source_local = None;
      ensemble_0 = None;
      snapshots = None;
      model_local = None;
      free_run = None;
      initial_condition = None;
      Nens = None;

      var_names = ['UG0','VG0','TG0','TRG0','PSG0','UG1','VG1','TG1','TRG1','PSG1'];

      x_ic = None;
      X = None;
      
      def __init__(self, path, gs, Nens, par=True):
          self.path = path;
          self.res = gs.name;
          self.source_model = '../models/speedy/'+self.res+'/';
          self.source_local = self.path+'source_local/'
          self.ensemble_0 = self.path+'ensemble_0/'
          self.snapshots = self.path+'snapshots/'
          self.model_local = self.path+'model_local/';
          self.free_run = self.path+'free_run/'
          self.initial_condition = self.path+'initial_condition/'
          self.Nens = Nens;
          self.test = 1;
          self.par = par;
          self.gs = gs;
          self.create_all_folders();
          
      def clear_all_folders(self):
          os.system(f'rm -rf {self.ensemble_0}');
          os.system(f'rm -rf {self.source_local}');
          os.system(f'rm -rf {self.model_local}');
          
      def create_all_folders(self):

          #
          os.system(f'mkdir {self.snapshots} ;mkdir {self.ensemble_0} ;mkdir {self.source_local} ;mkdir {self.free_run} ;mkdir {self.initial_condition} ;mkdir {self.model_local}');
          
          #
          Nens = self.Nens;
          for e in range(0, Nens):
              os.system(f'mkdir {self.ensemble_0}ens_{e}');
          
          #
          os.system(f'cp {self.source_model}* {self.source_local}');
          
      
      def set_resol_variables(self):
          self.var_resol = [];
          nvar = len(self.var_names);
          for v in range(0, nvar):
              if 'PSG' in self.var_names[v]:
                  self.var_resol.append([self.gs.get_resolution(self.res), 1]);
              else:
                  self.var_resol.append([self.gs.get_resolution(self.res), 8]);

      
      def set_time_integration(self,args):
          nmonths = args[0];
          days = args[1];
          restart = args[2];
          self.create_cls_instep_file(nmonths, days, restart);
          self.create_cls_indyns_file();
          os.system(f'mv cls_instep.h {self.source_local}cls_instep.h; mv cls_indyns.h {self.source_local}cls_indyns.h; cd {self.source_local}/ ; sh compile.sh>out.txt;');

      def load_netcdf_file(self,nc_f):
           var_names = self.var_names;
           n_vars = len(var_names);
           x = [];
           #print('dimension {0}'.format(x.shape));
           nc_fid = Dataset(nc_f, 'r');
           for v in range(0,n_vars):
               var_name = var_names[v];
               #print(var_name);
               var_full = nc_fid.variables[var_name][:];
               if 'TRG' in var_names[v]:
                   x.append(var_full[0,:,:,:]);
               elif 'PSG' in var_names[v]:
                   x.append(var_full[:,:]); #ask for 0 or -1
               else:
                   x.append(var_full[:,:,:]);
                #elias
                   #print('var_name is {0} and shape reads {1}'.format( var_name,var_data.shape));
                   #var_data = var_data.T; #The model stores the info lat x lon
           #print(x);
           nc_fid.close();
           return x;
           
      
      
      def create_initial_condition(self, ini0):
          #nmonths = args[0];
          #days = args[1];
          #restart = args[2];
          
          self.set_time_integration(ini0);
          
          os.system(f'cp {self.source_local}* {self.model_local}');          
          
          os.system(f'cd {self.model_local}; sh remove_grad_ctl.sh; ./imp.exe>out.txt');
          nc_f = self.model_local+'/ensemble_member.nc';
          gr_f = self.model_local+'/fort.10';
          
          os.system(f'cp {nc_f} {self.initial_condition}initial_condition.nc');
          os.system(f'cp {gr_f} {self.initial_condition}fort.3');
      
      
      def load_initial_condition(self):
          nc_f = self.initial_condition+'initial_condition.nc';
          self.x_ic = self.load_netcdf_file(nc_f);
                  
                  
      def create_perturbed_ensemble(self, per, Nens):
          np.random.seed(seed=10); #To replicate the initial ensemble
          self.load_initial_condition();
          n_vars = len(self.var_names);
          self.X = [];
          for e in range(0, Nens):
              Xe = [];
              for v in range(0, n_vars):
                  if self.test == 1: print('* Creating xb^[{0}] var {1}'.format(e, self.var_names[v]));
                  dim_v = self.x_ic[v].shape
                  ne = np.prod(dim_v);
                  x_e_v = self.x_ic[v].reshape((ne,1)) + per*np.random.randn(ne,1)*self.x_ic[v].reshape((ne,1));
                  Xe.append(x_e_v.reshape(dim_v));
                  #print(Xe[v].shape)
              fn_e = self.ensemble_0+f'ens_{e}/'+f'ensemble_member.nc';
              self.map_state_netcdf(Xe, fn_e); #Creating perturbed member in ensemble_0 local folder
              self.X.append(Xe);
              if self.test==1: print('* ENDJ - Stored the {0}-th ensemble member in = {1}'.format(e, fn_e));

      
      def map_state_netcdf(self, xs, fn):
      
          ds = Dataset(fn, 'w', format='NETCDF4');
          ntr = ds.createDimension('ntr',1);
          lev = ds.createDimension('lev', self.gs.lev);
          lat = ds.createDimension('lat', self.gs.lat);
          lon = ds.createDimension('lon', self.gs.lon);
          
          UG0 = ds.createVariable('UG0', np.float64, ('lev','lat', 'lon'));
          VG0 = ds.createVariable('VG0', np.float64, ('lev','lat', 'lon'));
          TG0 = ds.createVariable('TG0', np.float64, ('lev','lat', 'lon'));
          TRG0 = ds.createVariable('TRG0', np.float64, ('ntr','lev','lat', 'lon'));
          PSG0 = ds.createVariable('PSG0', np.float64, ('lat', 'lon',))
          UG1 = ds.createVariable('UG1', np.float64, ('lev','lat', 'lon'));
          VG1 = ds.createVariable('VG1', np.float64, ('lev','lat', 'lon'));
          TG1 = ds.createVariable('TG1', np.float64, ('lev','lat', 'lon'));
          TRG1 = ds.createVariable('TRG1', np.float64, ('ntr','lev', 'lat', 'lon'))
          PSG1 = ds.createVariable('PSG1', np.float64, ('lat', 'lon',))
          
          #print([UG0.shape, xs[0].shape]);
          
          UG0[:,:,:] = xs[0][:,:,:];
          VG0[:,:,:] = xs[1][:,:,:];
          TG0[:,:,:] = xs[2][:,:,:];
          TRG0[0,:,:,:] = xs[3][:,:,:];
          PSG0[:,:] = xs[4][:,:];
          UG1[:,:,:] = xs[5][:,:,:];
          VG1[:,:,:] = xs[6][:,:,:];
          TG1[:,:,:] = xs[7][:,:,:];
          TRG1[0,:,:,:] = xs[8][:,:,:];
          PSG1[:,:] = xs[9][:,:];
          
          ds.close();

###############################################################################
###############################################################################
#Routines regarding ensembles
###############################################################################
###############################################################################
          
      def copy_reference_restart(self, Nens):
          gr_f = self.initial_condition+'/fort.3';
          for e in range(0, Nens):
              os.system(f'cp {gr_f} {self.ensemble_0}ens_{e}/fort.3');

      
      def update_model_ensemble_folders(self, Nens):
          for e in range(0, Nens):
              ensemble_path = self.ensemble_0+'ens_'+str(e)+'/';
              os.system(f'cp {self.source_local}/* {ensemble_path}');


      def update_members_ensemble_folders(self, Nens):
          for e in range(0, Nens):
              if self.test==1: print('* ENDJ - Updating ensemble member {0}'.format(e));
              ensemble_path = self.ensemble_0+'ens_'+str(e)+'/';
              ornc_path = self.ensemble_0+'ensemble_member_'+str(e)+'.nc';
              orgr_path = self.ensemble_0+'fort_'+str(e)+'.3';
              os.system(f'mv {ornc_path} {ensemble_path}ensemble_member.nc');
              os.system(f'mv {orgr_path} {ensemble_path}fort.3');


      def perform_forecast(self, e):
          if self.test==1: print('* ENDJ - Performing forecast ensemble member {0}'.format(e));
          ensemble_path = f'{self.ensemble_0}ens_{e}/';
          os.system(f'cd {ensemble_path}; sh remove_grad_ctl.sh; ./imp.exe>out.txt; mv fort.10 fort.3');
          return ensemble_path;
      
      def forecast_ensemble(self, Nens):
          for e in range(0, Nens):
              #print('* Forecast the ensemble member {0}'.format(e));
              self.perform_forecast(e);
              
      def forecast_ensemble_parallel(self, Nens):
          pool = mp.Pool(mp.cpu_count())
          results = pool.map(self.perform_forecast, [e for e in range(0, Nens)]);
          pool.close();
          #print(results);
      
      def collect_ensemble_members(self):
          Nens = self.Nens;
          for e in range(0, Nens):
              os.system(f'mv {self.ensemble_0}/ens_{e}/ensemble_member.nc {self.ensemble_0}ensemble_member_{e}.nc');
              os.system(f'mv {self.ensemble_0}/ens_{e}/fort.3 {self.ensemble_0}fort_{e}.3');
              os.system(f'rm -rf {self.ensemble_0}/ens_{e}');
          if self.test==1: print('* ENDJ - All ensemble members have been collected');

      def load_ensemble(self):
          self.X = [];
          Nens = self.Nens;
          for e in range(0, Nens):
              ensemble_path = f'{self.ensemble_0}ens_{e}/ensemble_member.nc';
              xe = self.load_netcdf_file(ensemble_path);
              self.X.append(xe);

      
      def create_initial_ensemble(self, ini0, args, per, Nens):
          

          #self.update_members_ensemble_folders(Nens);
          self.set_time_integration(ini0);
          self.update_model_ensemble_folders(Nens);
          
          self.create_perturbed_ensemble(per, Nens);
          self.copy_reference_restart(Nens);
          
          if self.par:
             self.forecast_ensemble_parallel(Nens);
          else:
             self.forecast_ensemble(Nens);

          self.set_time_integration(args);
          self.update_model_ensemble_folders(Nens);
          
          
          if self.test == 1: print('* ENDJ - The initial ensemble has been created Nens = {0}'.format(Nens));
          
###############################################################################
###############################################################################
#Routines regarding reference and background trajectories
###############################################################################
###############################################################################
     
      def create_reference_snapshots(self, ini0, args, M):
          nc_f = self.initial_condition+'initial_condition.nc';
          gr_f = self.initial_condition+'fort.3';
          self.create_snapshots(nc_f, gr_f, self.snapshots, 'reference_solution', ini0, args, M);
     
      def create_free_run(self, ini0, args, M):
          self.load_ensemble();    
          xb = self.compute_snapshot_mean(self.X);
          nc_f = self.model_local+'ensemble_member.nc';
          self.map_state_netcdf(xb, nc_f);
          gr_f = self.initial_condition+'fort.3';
          os.system(f'cp {gr_f} {self.model_local}fort.10');
          self.create_snapshots(nc_f, gr_f, self.free_run, 'free_run', ini0, args, M, dyn_cons = False);
      
      def create_snapshots(self, nc_f, gr_f, folder_dest, name_conv, ini0, args, M, dyn_cons = True):
          
          #To months
          #args = [3, 6, 1]; #[months, hours, restart = 1(yes, we read the nc]
          #OJOOO
          
          if dyn_cons:
             self.set_time_integration(ini0);
             os.system(f'cp {nc_f} {folder_dest}{name_conv}.nc');
             os.system(f'cp {gr_f} {folder_dest}fort.3');
             os.system(f'cp {self.source_local}/* {self.model_local}');
             os.system(f'cp {folder_dest}{name_conv}.nc {self.model_local}ensemble_member.nc');
             os.system(f'cp {folder_dest}fort.3 {self.model_local}fort.3');
             os.system(f'cd {self.model_local}; sh remove_grad_ctl.sh; ./imp.exe>out.txt');
          
          #Copy the propagated model
          os.system(f'mv {self.model_local}ensemble_member.nc {folder_dest}{name_conv}_0.nc');
          os.system(f'mv {self.model_local}fort.10 {folder_dest}fort_0.3');
          #Change time integration to days
          self.set_time_integration(args);
          #Update the model folder
          os.system(f'cp {self.source_local}/* {self.model_local}');
           
          #Let's take the snapshots
          for s in range(0, M):
              
              print('* Working on snapshot {0}'.format(s));
              
              #Copy the reference solution to the Speedy model
              os.system(f'cp {folder_dest}{name_conv}_{s}.nc {self.model_local}ensemble_member.nc');
              os.system(f'cp {folder_dest}fort_{s}.3 {self.model_local}fort.3');
          
              #Enter in the reference model and run the model
              os.system(f'cd {self.model_local}; sh remove_grad_ctl.sh; ./imp.exe>out.txt');
          
              #Copy the propagated model
              os.system(f'mv {self.model_local}ensemble_member.nc {folder_dest}{name_conv}_{s+1}.nc');
              os.system(f'mv {self.model_local}fort.10 {folder_dest}fort_{s+1}.3');
          

          if self.test == 1: print(f'* ENDJ - Finishing creating the {name_conv} trajectory for M = {M}');
 

      def get_empty_state(self):
          n_vars = len(self.var_names);
          X = [];
          for v in range(0, n_vars):
              var_resol = self.var_resol[v];
              lat, lon = var_resol[0];
              lev = var_resol[1];
              if lev>1:
                 X.append(np.zeros((lev, lat, lon)));
              else:
                 X.append(np.zeros((lat, lon)));
          return X;
                 
      def map_vector_state(self, X_all, e):
          X = self.get_empty_state();
          for msk_cor, X_block in zip(self.mask_cor, X_all):
              ini = 0;
              for var in msk_cor:
                  var_info = var[0];
                  var_reso = var[1];
                  var_index = var_info[0];
                  var_level = var_info[1];
                  lat, lon = var_reso[0], var_reso[1];
                  X_v = X[var_index];
                  n = lat * lon;
                  fin = ini+n;
                  #print(f'ini = {ini} fin = {fin} n = {n} var_level = {var_level} e = {e}');
                  if 'PSG' in self.var_names[var_index]:
                     X_v = X_block[ini:fin, e].reshape((lat, lon));
                  else:
                     X_v[var_level,:,:] = X_block[ini:fin, e].reshape((lat, lon));
                  X[var_index] = X_v;
                  ini+=n;
          return X;  
          
      def compute_snapshot_mean(self, X):
        var_names = self.var_names;
        n_vars = len(var_names);
        samples = len(X);
        xm = [];
        for v in range(0, n_vars):
            x = np.zeros(X[0][v].shape);
            for e in range(0, samples):
                x+=X[e][v];
            x/=samples;
            xm.append(x);
        return xm;
        
###############################################################################
###############################################################################
#Routines regarding model parameters
###############################################################################
###############################################################################      
        
      def create_cls_instep_file(self, nmonths, days, restart):
          f = open('cls_instep.h','w')
          f.write('      NMONTS = '+str(nmonths)+'\n');
          f.write('      NDAYSL = '+str(days)+'\n');
          f.write('      HOURS = 0'+'\n');
          f.write('      ISTART = '+str(restart)+'\n');
          f.write('      NSTPPR = 6'+'\n');
          f.write('      NSTOUT = -1'+'\n');
          f.write('      IDOUT  = 0'+'\n');
          f.write('      NMONRS = -1'+'\n');
          f.write('      ISEASC = 1'+'\n');
          f.write('      IYEAR0 = 1979'+'\n');
          f.write('      IMONT0 = 1'+'\n');
          f.write('      NSTRAD = 3'+'\n');
          f.write('      NSTRDF = 0'+'\n');
          f.write('      INDRDF = 1'+'\n');
          f.write('      ICLAND = 1'+'\n');
          f.write('      ICSEA  = 0'+'\n');
          f.write('      ICICE  = 1'+'\n');
          f.write('      ISSTAN = 1'+'\n');
          f.write('      ISSTY0 = 1870'+'\n');
          f.write('      ISST0  = (IYEAR0-ISSTY0)*12+IMONT0'+'\n');
          f.write('      LPPRES = .true.'+'\n');
          f.write('      LCO2 = .false.'+'\n');
          
          if (self.res=='t21') or (self.res=='t30'):
             f.write('      NSTEPS = 36'+'\n');
             f.write('      NSTDIA = 36*5'+'\n');
          elif self.res=='t47':
             f.write('      NSTEPS = 72'+'\n');
             f.write('      NSTDIA = 72*5'+'\n');          
          elif self.res=='t63':
             f.write('      NSTEPS = 96'+'\n');
             f.write('      NSTDIA = 96*5'+'\n');          
          elif self.res=='t106':
             f.write('      NSTEPS = 106'+'\n');
             f.write('      NSTDIA = 106*5'+'\n');          
          else:
             print('* Invalid resolution '+self.res);
             exit();
          
          f.write('      BLOCKHOURS = 24./FLOAT(NSTEPS)'+'\n');
 
          f.close();
          print('* cls_instep.h has been created');


      def create_cls_indyns_file(self):
          f = open('cls_indyns.h','w')
          f.write('      GAMMA  = 6.'+'\n');
          f.write('      HSCALE = 7.5'+'\n');
          f.write('      HSHUM  = 2.5'+'\n');
          f.write('      REFRH1 = 0.7'+'\n');
          f.write('      THDS   = 12.'+'\n');
          f.write('      TDRS   = 24.*30.'+'\n');
          
          if (self.res=='t21') or (self.res=='t30'):
             f.write('      THD    = 2.4'+'\n');
             f.write('      THDD   = 2.4'+'\n');
          elif self.res=='t47':
             f.write('      THD    = 0.5'+'\n');
             f.write('      THDD   = 0.5'+'\n');      
          elif self.res=='t63':
             f.write('      THD    = 0.5'+'\n');
             f.write('      THDD   = 0.5'+'\n');         
          elif self.res=='t106':
             f.write('      THD    = 0.1'+'\n');
             f.write('      THDD   = 0.1'+'\n');                
          else:
             print('* Invalid resolution '+self.res);
             exit();
 
          f.close();
          print('* cls_indyns.h has been created');
          
###############################################################################
###############################################################################
#Routines regarding variable relationships
###############################################################################
###############################################################################
      def create_rel_per_level(self):
          self.mask_cor = [];
          nvar = len(self.var_names);
          for l in range(0, 8):
              var_l = [];
              for v in range(0, nvar):
                  ml = self.mask_lev[v];
                  if ml[l]:
                     var_l.append(([v,l],(self.gs.get_resolution(self.res))));
              self.mask_cor.append(var_l);
      
      def create_per_variable(self):
          self.mask_cor = [];
          nvar = len(self.var_names);
          for l in range(0, 8):
              for v in range(0, nvar):
                  ml = self.mask_lev[v];
                  if ml[l]:
                     self.mask_cor.append([([v,l],(self.gs.get_resolution(self.res)))]);
                     
      def create_rel_per_variable(self):
          self.mask_cor = [];
          nvar = int(len(self.var_names)/2);
          print(f'* nvar = {nvar}');
          for l in range(0, 8):
              for v in range(0, nvar):
                  ml = self.mask_lev[v];
                  if ml[l]:
                     self.mask_cor.append([([v,l],(self.gs.get_resolution(self.res))), ([v+nvar,l],(self.gs.get_resolution(self.res)))]);
                     
      def create_dummy(self):
          self.mask_cor = [[([0,0],(self.gs.get_resolution(self.res)))], [([1,0],(self.gs.get_resolution(self.res)))]];
                     
      def define_relations(self, option=1):
          var_all = [1,1,1,1,1,1,1,1];
          var_one = [1,0,0,0,0,0,0,0];
          var_two = [0,0,1,1,1,1,1,1];
          self.mask_lev = [var_all,var_all,var_all,var_two,var_one,var_all,var_all,var_all,var_two,var_one];
          if np.size(option)>1: self.mask_cor = option;
          if option==1: self.create_rel_per_level();
          if option==2: self.create_per_variable();
          if option==3: self.create_rel_per_variable();
          if option==100: self.create_dummy();
      
      def load_settings(self, path, args, no_test=True):
          os.system(f'cp -r {path}snapshots/* {self.snapshots}');
          os.system(f'cp -r {path}source_local/* {self.source_local}')
          os.system(f'cp -r {path}free_run/* {self.free_run}');
          os.system(f'cp -r {path}initial_condition/* {self.initial_condition}');
          os.system(f'cp -r {path}model_local/* {self.model_local}');
          
          self.time_snapshots = f'{self.path}time_snapshots/';
          os.system(f'mkdir {self.time_snapshots}');
          
          if no_test: self.set_time_integration(args);
          #
          Nens = self.Nens;
          for e in range(0, Nens):
              os.system(f'cp {self.source_local}/* {self.ensemble_0}ens_{e}/');
              os.system(f'cp {path}/ensemble_0/ensemble_member_{e}.nc {self.ensemble_0}ens_{e}/ensemble_member.nc');
              os.system(f'cp {path}/ensemble_0/fort_{e}.3 {self.ensemble_0}ens_{e}/fort.3');
          
          self.set_resol_variables();    
          
          #
      def set_resol_variables(self):
          self.var_resol = [];
          nvar = len(self.var_names);
          for v in range(0, nvar):
              if 'PSG' in self.var_names[v]:
                  self.var_resol.append([self.gs.get_resolution(self.res), 1]);
              else:
                  self.var_resol.append([self.gs.get_resolution(self.res), 8]);
          
         
