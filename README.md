# AMLCS
<img src="https://aml-cs.github.io/images/logo.jpg">

The `AMLCS Package` is a `FORTRAN`+`Python` based package which holds the following Data Assimilation methods for the SPEEDY model: 

1. Szunyogh, I., Kostelich, E. J., Gyarmati, G., Kalnay, E., Hunt, B. R., Ott, E., ... & Yorke, J. A. (2008). A local ensemble transform Kalman filter data assimilation system for the NCEP global model. Tellus A: Dynamic Meteorology and Oceanography, 60(1), 113-130.
2. Anderson, J. L. (2012). Localization and sampling error correction in ensemble Kalman filter data assimilation. Monthly Weather Review, 140(7), 2359-2371.
3. Nino-Ruiz, E. D., Sandu, A., & Deng, X. (2019). A parallel implementation of the ensemble Kalman filter based on modified Cholesky decomposition. Journal of Computational Science, 36, 100654.

If you use this repo, please cite the following paper:

1. Nino-Ruiz, E. D., Guzman-Reyes, L. G., & Beltran-Arrieta, R. (2020). An adjoint-free four-dimensional variational data assimilation method via a modified Cholesky decomposition and an iterative Woodbury matrix formula. Nonlinear Dynamics, 99(3), 2441-2457.
2. Nino-Ruiz, E. D., Sandu, A., & Deng, X. (2018). An ensemble Kalman filter implementation based on modified Cholesky decomposition for inverse covariance matrix estimation. SIAM Journal on Scientific Computing, 40(2), A867-A886.

# Code Documentation

## Description
The following code defines several variables and objects necessary for a data assimilation process:

- `gs`: an instance of `grid_resolution` class created with `res_name` parameter.
- `data_obs`: a string of observed data separated by commas extracted from the `err_obs` column of `df_par` dataframe.
- `err_obs`: a list of floating-point numbers created by converting the values in `data_obs`.
- `plac_obs`: a string of observation positions separated by commas extracted from the `obs_plc` column of `df_par` dataframe.
- `obs_plc`: a list of integers created by converting the values in `plac_obs`.
- `ob`: an instance of `observation` class created with `err_obs` and `obs_plc` as parameters.
- `nm`: an instance of `numerical_model` class created with `path_method`, `gs`, `Nens`, and `par` parameters.
- `rs`: an instance of `reference_solution` class created with `nm` and `M` parameters.
- `seq_da`: an instance of a class that inherits from `sequential_method` created with `nm`, `infla`, `Nens`, and `method` parameters.
- `em_bck`: an instance of `error_metric` class created with `nm`, `"bck"`, and `M` as parameters.
- `em_ana`: an instance of `error_metric` class created with `nm`, `"ana"`, and `M` as parameters.
- `tm_bck`: an instance of `time_metric` class created with `nm`, `"bck"`, and `M` as parameters.
- `tm_ana`: an instance of `time_metric` class created with `nm`, `"ana"`, and `M` as parameters.

## Parameters
- `res_name`: a string that defines the grid resolution.
- `df_par`: a dataframe that contains observational and model parameters.
- `path_method`: a string that defines the path to the numerical model.
- `Nens`: an integer that defines the number of ensemble members.
- `par`: a dictionary that contains model parameters.
- `option_mask`: a string that defines the relation between state variables and model parameters.
- `exp_settings`: a dictionary that contains experimental settings.
- `args`: a dictionary that contains additional arguments.
- `r`: a float that defines the sub-domain computation radius.
- `M`: an integer that defines the number of assimilation cycles.
- `infla`: a float that defines the inflation factor.
- `method`: a string that defines the assimilation method.
- `s`: a float that defines the observational noise standard deviation.
- `list_k`: a list that contains the assimilation cycle numbers to store.

## Example Usage
```python
gs = grid_resolution(res_name);
data_obs = df_par['err_obs'].iloc[0].strip().split(',');
err_obs = [float(v) for v in data_obs];
plac_obs = df_par['obs_plc'].iloc[0].strip().split(',');
obs_plc = [int(v) for v in plac_obs];
ob = observation(err_obs, obs_plc);

nm = numerical_model(path_method, gs, Nens, par=par);
nm.define_relations(option=option_mask);
nm.load_settings(exp_settings, args);
     
gs.create_mesh(nm);
gs.compute_sub_domains(r);
rs = reference_solution(nm, M);
    
seq_da = sequential_method(method).get_instance(nm, infla, Nens);
    
ob.build_observational_network(gs, nm, s=s);
ob.build_synthetic_observations(nm, rs, M); 
    
em_bck = error_metric(nm, 'bck', M);
em_ana = error_metric(nm, 'ana', M);
tm_bck =  time_metric(nm, 'bck');
tm_ana =  time_metric(nm, 'ana');
    
for k in range(0, M):
    seq_da.load_background_ensemble();
        
    tm_bck.start_time();
    seq_da.prepare_background();
    tm_bck.check_time();
        
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
        
    seq_da.perform_forecast();
```
    
# Adding a new model
 
## Incorporating Other DA Algorithms and Dynamical Models in AMLCS-DA Package

The AMLCS-DA package provides a flexible framework for incorporating other data assimilation (DA) algorithms and dynamical models.

## Adding a New DA Algorithm

To add a new DA algorithm, you need to create a new class that inherits from the `SequentialAlgorithm` abstract class. This class should implement the following methods:

- `prepare_background()`: prepares the background state for assimilation
- `prepare_analysis(observation, time_step)`: prepares the analysis state for assimilation based on the current observation and time step
- `perform_assimilation(observation)`: performs the assimilation step using the current observation
- `perform_forecast()`: performs the forecast step

## Adding a New Dynamical Model

To add a new dynamical model, you need to create a new class that inherits from the `NumericalModel` abstract class. This class should implement the following methods:

- `step_forward(state, time_step)`: steps forward the model state by the given time step
- `define_relations(option)`: defines the relations between the state variables of the model
- `load_settings(settings)`: loads the settings of the model

## Example

Here is an example of how to add a new DA algorithm and dynamical model in the AMLCS-DA package:

```python
from AMLCS.da.sequential_algorithm import SequentialAlgorithm
from AMLCS.da.numerical_model import NumericalModel

class MyAlgorithm(SequentialAlgorithm):
    # Implement the required methods here

class MyModel(NumericalModel):
    # Implement the required methods here

# Instantiate the new DA algorithm and dynamical model
algorithm = MyAlgorithm()
model = MyModel()

# Use the new algorithm and model in the AMLCS-DA package
da_system = AMLCS_DA(algorithm, model)
```

# Adding the Lorenz 96 model to AMLCS-DA:

To add a new dynamical model to the AMLCS-DA package, you need to define the model equations and implement them as a class in the `dynamical_model.py` file. This class should inherit the abstract class `DynamicalModel` and implement the following methods:

- `__init__(self, n)`: initialize the model parameters and state with the size `n`
- `integrate(self, t0, tf, dt)`: integrate the model equations from `t0` to `tf` with time step `dt`
- `get_state(self)`: return the current state of the model
- `set_state(self, state)`: set the state of the model to `state`

Once you have defined your new dynamical model class, you can use it in the AMLCS-DA package by specifying it as an argument to the `numerical_model` constructor. For example, if you have defined a Lorenz96 model in a class called `Lorenz96Model`, you can use it as follows:

```python
from amlcsda.numerical_model import numerical_model
from my_dynamical_model import Lorenz96Model

# Define the model parameters
path_method = 'EnSRF'
gs = 0.05
Nens = 20
par = {'F': 8.0, 'h': 1.0}

# Create an instance of the Lorenz96 model
model = Lorenz96Model(n=40)

# Create an instance of the numerical model using the Lorenz96 model
nm = numerical_model(path_method, gs, Nens, par=par, dynamical_model=model)

# Run the model
nm.run_model()
```

In this example, we create an instance of the Lorenz96Model class with n=40 and pass it as an argument to the numerical_model constructor. We then run the model using the run_model method of the numerical_model object.

```python
class Lorenz96Model(DynamicalModel):
    def __init__(self, n):
        self.n = n
        self.F = 8.0
        self.h = 1.0
        self.x = np.zeros(self.n)
        
    def integrate(self, t0, tf, dt):
        # Integration code here
        
    def get_state(self):
        return self.x
    
    def set_state(self, state):
        self.x = state
```

# Class Diagram

<img src="https://raw.githubusercontent.com/enino84/AMLCS/main/Class_Diagram_AMLCS.png">
