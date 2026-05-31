python3 amlcs_pre.py amlcs_pre_t21.csv
python3 amlcs_da.py amlcs_da_t21_0.csv   # EnKF_MC_obs
python3 amlcs_da.py amlcs_da_t21_1.csv   # LETKF
python3 amlcs_da.py amlcs_da_t21_2.csv   # LEnKF
python3 amlcs_da.py amlcs_da_t21_3.csv   # EnKF_Sh_Binv_MSE   (shrinkage inverse-cov, MSE)
python3 amlcs_da.py amlcs_da_t21_4.csv   # EnKF_Sh_Binv_Stein (shrinkage inverse-cov, Stein)
python3 amlcs_da.py amlcs_da_t21_5.csv   # EnKF_Sh_Binv_DA    (shrinkage inverse-cov, DA-aware)
