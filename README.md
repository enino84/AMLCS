# AMLCS
<img src="https://aml-cs.github.io/images/logo.jpg">

The `AMLCS Package` is a `FORTRAN`+`Python` based package which holds the following Data Assimilation methods for the SPEEDY model: 

1. Szunyogh, I., Kostelich, E. J., Gyarmati, G., Kalnay, E., Hunt, B. R., Ott, E., ... & Yorke, J. A. (2008). A local ensemble transform Kalman filter data assimilation system for the NCEP global model. Tellus A: Dynamic Meteorology and Oceanography, 60(1), 113-130.
2. Anderson, J. L. (2012). Localization and sampling error correction in ensemble Kalman filter data assimilation. Monthly Weather Review, 140(7), 2359-2371.
3. Nino-Ruiz, E. D., Sandu, A., & Deng, X. (2019). A parallel implementation of the ensemble Kalman filter based on modified Cholesky decomposition. Journal of Computational Science, 36, 100654.

If you use this repo, please cite the following papers:

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

The class diagram represents the organization and interconnection of the classes within the AMLCS package. AMLCS is a Python toolbox for testing sequential Data Assimilation (DA) methods using the SPEEDY model, which is a well-known Atmospheric General Circulation model in the DA community.

The classes in the diagram are divided into five categories: grid, numerical model, DA methods, input/output, and plot/covariance. The grid class contains a single Grid class that defines the grid resolution used in the numerical model. The numerical model category includes two classes: NumericalModel and ICTools. NumericalModel is the core class that defines the atmospheric model, while ICTools contains methods for creating initial conditions.

The DA methods category includes four classes: DA, LETKF, LEKF, and ETKF. These classes implement well-known ensemble-based DA methods, such as the Local Ensemble Transform Kalman Filter (LETKF) and the Local Ensemble Kalman Filter (LEKF). The input/output category includes a single IO class that contains methods for reading and writing data. The plot/covariance category includes two classes: PlotTools and Covariance. PlotTools contains methods for visualizing data, while Covariance contains methods for computing the background error covariance matrix.

Finally, the AMLCS package includes three scripts: AMLCS_PRE, AMLCS_MAIN, and AMLCS_POST, which are responsible for performing the preprocessing, assimilation, and post-processing steps, respectively. These scripts utilize the classes defined in the package to perform the various DA operations.

<img src="https://raw.githubusercontent.com/enino84/AMLCS/main/Class_Diagram_AMLCS.png">

## Class definitions

The AMLCS-DA package includes the following classes:

- `Grid`: Defines the spatial discretization of the domain using a specific grid resolution. It contains information about the domain size, grid spacing, and the number of grid points.
- `NumericalModel`: Implements the forecast step of the DA cycle using the AT-GCM SPEEDY by default. It is responsible for integrating the model forward in time, creating initial conditions for the ensemble members, and generating the reference trajectory. This class also provides methods for collecting ensemble members, creating free runs, and performing the forecast step of the DA cycle.
- `DA`: An abstract class that defines the basic structure of a DA method. It contains methods for updating the background state using the observation and the forecast error, and for computing the analysis error covariance matrix.
- `LETKF`: A class that implements the Local Ensemble Transform Kalman Filter method, which uses the local transformation of the ensemble members to update the background state.
- `LEKF`: A class that implements the Local Ensemble Kalman Filter method, which updates the background state using the ensemble mean and the perturbations around it.
- `IO`: A class that provides methods for reading and writing input and output files. It includes methods for reading the configuration file, the observations file, and the output files.
- `PlotTools`: A class that provides methods for generating plots of the results. It includes methods for plotting the time evolution of the analysis and forecast error, the spatial distribution of the analysis error, and the cross-section of the analysis error at different pressure levels.
- `Covariance`: A class that defines the background error covariance matrix structure. It includes methods for computing the correlation matrix based on different assumptions about the error structure.
- `AMLCS_PRE`: A class that defines the preprocessing step of the DA cycle. It reads the configuration file and sets up the initial conditions for the forecast and assimilation steps.
- `AMLCS_MAIN`: A class that runs the main loop of the DA cycle. It reads the observations file and performs the forecast and assimilation steps for each time step.
- `AMLCS_POST`: A class that postprocesses the output files. It generates plots of the results and computes the statistics of the analysis error.

The relationships between the classes are as follows:

- `Grid` is used by `NumericalModel` to define the spatial discretization of the domain.
- `NumericalModel` is used by all the DA methods (`LETKF`, `LEKF`) to perform the forecast step and generate the initial conditions for the ensemble members.
- `LETKF` and `LEKF` inherit from the abstract class `DA`, which defines the basic structure of a DA method.
- `IO` is used by `AMLCS_PRE` and `AMLCS_MAIN` to read input files and write output files.
- `PlotTools` is used by `AMLCS_MAIN` and `AMLCS_POST` to generate plots of the results.
- `Covariance` is used by `LETKF` and `LEKF` to compute the background error covariance matrix.
- `AMLCS_PRE` sets up the initial conditions for the forecast and assimilation steps, and `AMLCS_MAIN` performs the forecast and assimilation steps for each time step. `AMLCS_POST` postprocesses the output files generated by `AMLCS_MAIN`.


# Sequential Diagram

This sequential diagram shows the steps of the AMLCS (Adaptive Multilevel Monte Carlo Localization and Assimilation System) algorithm. The AMLCS algorithm is a data assimilation method used in numerical weather prediction and other fields to combine simulation models and real-world observational data in order to produce more accurate predictions.

The diagram begins with the user running the AMLCS program. The program first instantiates the numerical model (B), which is used to represent the simulated system being studied. The program then instantiates the grid resolution (C) to define the spatial resolution of the model. The grid resolution is then used to create a mesh and compute subdomains for the numerical model.

The program then defines relations between the numerical model and the grid resolution. The program loads the settings for the experiment and instantiates the reference solution (E), which represents the "true" state of the system being simulated. The program also instantiates the sequential method (F), which is used to perform the data assimilation.

Next, the user provides observational data to the program. The observation class (D) is used to build the observational network and generate synthetic observations.

The program then enters a loop for each time step. At each time step, the sequential method loads the background ensemble, prepares the background, and checks the background time. The program then starts the analysis time and prepares the analysis using the observational data. The analysis is then performed by assimilating the observational data into the model using the observation class. The program then checks the analysis time and computes the error step using the error metric class (G). The error results are then stored, along with the time results, using the time metric class (H). The program then checks the time store and performs a forecast for the next time step.

The loop continues until all time steps have been processed. Finally, the program executes the complete AMLCS algorithm and returns the results to the user.

<img src="https://raw.githubusercontent.com/enino84/AMLCS/main/Seq_Diagram_AMLCS.png">
