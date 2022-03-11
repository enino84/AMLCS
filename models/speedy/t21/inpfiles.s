# Script to link fortran units to input files; $1 = resolution (t21, t30, ..)

SB=../data/bc/$1/clim
SC=../data/bc/$1/anom
SH=../hflux	

 cp $SB/orog_lsm_alb.${1}.grd         fort.20
 cp $SB/sst_clim.${1}.sea.grd     fort.21
 cp $SB/seaice_clim.${1}.sea.grd  fort.22
 cp $SB/surfv_st3_clim.${1}.land.grd   fort.23	
 cp $SB/sndep_clim.${1}.land.grd  fort.24
 cp $SB/veget.${1}.land.grd           fort.25
 cp $SB/soilw_clim.${1}.land.grd  fort.26

 cp    $SC/hadisst_anom_1_1.${1}.grd fort.30	
	
 cp    $SH/hflux_speedy_ver41.5_1979_2008_clim.grd  fort.31

