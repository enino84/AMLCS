      MODULE read_write_netcdf
      IMPLICIT NONE
      include 'netcdf.inc'
      CONTAINS
      
C--   Multiple steps


      SUBROUTINE write_variables(UG0,UG1,VG0,VG1,TG0,TG1,PSG0,PSG1,
     +     TRG0,TRG1,IX,IL,KX,NTR)
      INTEGER, INTENT(IN) :: IX,IL,KX,NTR
      REAL, INTENT(IN)  ::  UG0(IX,IL,KX)
      REAL, INTENT(IN)  ::  UG1(IX,IL,KX)
      REAL, INTENT(IN)  ::  VG0(IX,IL,KX)
      REAL, INTENT(IN)  ::  VG1(IX,IL,KX)
      REAL, INTENT(IN)  ::  TG0(IX,IL,KX)
      REAL, INTENT(IN)  ::  TG1(IX,IL,KX)
      REAL, INTENT(IN)  ::  PSG0(IX,IL)
      REAL, INTENT(IN)  ::  PSG1(IX,IL)
      REAL, INTENT(IN)  ::  TRG0(IX,IL,KX,NTR)
      REAL, INTENT(IN)  ::  TRG1(IX,IL,KX,NTR)

      character*(*) LEV_NAME, LAT_NAME, LON_NAME, REC_NAME
      parameter (LEV_NAME = 'level')
      parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
      parameter (REC_NAME = 'recname')

      INTEGER lon_dimid, lat_dimid, lev_dimid, rec_dimid
      INTEGER :: dim4d(4), dim3d(3), dim2d(2)
      integer nlons, nlats, nlevs

      character*(*) UG1_N, UG0_N, VG0_N, VG1_N, TG0_N, TG1_N
      character*(*) PSG0_N, PSG1_N, TRG0_N, TRG1_N
      parameter (UG0_N = 'UG0')
      parameter (UG1_N = 'UG1')
      parameter (VG0_N = 'VG0')
      parameter (VG1_N = 'VG1')
      parameter (TG0_N = 'TG0')
      parameter (TG1_N = 'TG1')
      parameter (PSG0_N = 'PSG0')
      parameter (PSG1_N = 'PSG1')
      parameter (TRG0_N = 'TRG0')
      parameter (TRG1_N = 'TRG1')

      INTEGER ncid, retval
      integer start4(4), cuent4(4), start3(3), cuent3(3)
      integer start2(2), cuent2(2)


      character*(*) FILE_NAME
      parameter (FILE_NAME='ensemble_member.nc')


      INTEGER :: lon_varid, lat_varid, lev_varid, rec_varid
      INTEGER :: ug0_id,ug1_id, vg0_id, vg1_id, tg0_id, tg1_id
      INTEGER :: psg0_id, psg1_id, trg0_id, trg1_id

      INTEGER :: dims4d, dims3d, dims2d
      parameter (dims4d=4, dims3d = 3, dims2d = 2)

      INTEGER :: dimids4d(dims4d), dimids3d(dims3d), dimids2d(dims2d)

      nlevs = KX
      nlats = IL
      nlons = IX

C--   we create the netcdf file
      retval = nf_create(FILE_NAME, NF_CLOBBER, ncid) 

C--   we define the longitude, the latitude, and the levels
      retval = nf_def_dim(ncid, lon_name, nlons, lon_dimid) 
      retval = nf_def_dim(ncid, lat_name, nlats, lat_dimid)
      retval = nf_def_dim(ncid, lev_name, nlevs, lev_dimid)
C--   undefined dimension for 4d variables
      retval = nf_def_dim(ncid, rec_name, nf_unlimited, rec_dimid)

C--   two dimensional fields
      dimids2d(1) = lon_dimid
      dimids2d(2) = lat_dimid

C--   three dimensional fields
      dimids3d(1) = lon_dimid
      dimids3d(2) = lat_dimid
      dimids3d(3) = lev_dimid

C--   four dimensional fields
      dimids4d(1) = lon_dimid
      dimids4d(2) = lat_dimid
      dimids4d(3) = lev_dimid
      dimids4d(4) = rec_dimid

      retval = nf_def_var(ncid, lon_name, nf_double, 1, lon_dimid, 
     +     lon_varid)
      retval = nf_def_var(ncid, lat_name, nf_double, 1, lat_dimid, 
     +     lat_varid)
      retval = nf_def_var(ncid, lev_name, nf_double, 1, lev_dimid, 
     +     lev_varid)

C--   defining 3d fields
      retval = nf_def_var(ncid, UG0_N, nf_double, dims3d, dimids3d, 
     +     ug0_id) 
      retval = nf_def_var(ncid, UG1_N, nf_double, dims3d, dimids3d, 
     +     ug1_id) 
      retval = nf_def_var(ncid, VG0_N, nf_double, dims3d, dimids3d, 
     +     vg0_id) 
      retval = nf_def_var(ncid, VG1_N, nf_double, dims3d, dimids3d, 
     +     vg1_id) 
      retval = nf_def_var(ncid, TG0_N, nf_double, dims3d, dimids3d, 
     +     tg0_id) 
      retval = nf_def_var(ncid, TG1_N, nf_double, dims3d, dimids3d, 
     +     tg1_id) 

C--   defining 2d fields
      retval = nf_def_var(ncid, PSG0_N, nf_double, dims2d, dimids2d, 
     +     psg0_id) 
      retval = nf_def_var(ncid, PSG1_N, nf_double, dims2d, dimids2d, 
     +     psg1_id) 

C--   defining 4d fields
      retval = nf_def_var(ncid, TRG0_N, nf_double, dims4d, dimids4d, 
     +     trg0_id) 
      retval = nf_def_var(ncid, TRG1_N, nf_double, dims4d, dimids4d, 
     +     trg1_id) 

C--   area for writting 2d fields
      start2(1) = 1
      start2(2) = 1
      cuent2(1) = nlons
      cuent2(2) = nlats

C--   cube for writting 3d fields
      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
      cuent3(1) = nlons
      cuent3(2) = nlats
      cuent3(3) = nlevs

C--   hypercube for writting 4d fields
      start4(1) = 1
      start4(2) = 1
      start4(3) = 1
      start4(4) = 1
      cuent4(1) = nlons
      cuent4(2) = nlats
      cuent4(3) = nlevs
      cuent4(4) = NTR

      print *, shape(UG1)
      print *, start3
      print *, cuent3

C--   we stop defining variables
      retval = nf_enddef(ncid)

C--   writting 3d fields
      retval = nf_put_vara_double(ncid, ug0_id, start3, cuent3, 
     +        UG0)
      retval = nf_put_vara_double(ncid, ug1_id, start3, cuent3, 
     +        UG1)

      retval = nf_put_vara_double(ncid, vg0_id, start3, cuent3, 
     +        VG0)
      retval = nf_put_vara_double(ncid, vg1_id, start3, cuent3, 
     +        VG1)

      retval = nf_put_vara_double(ncid, tg0_id, start3, cuent3, 
     +        TG0)
      retval = nf_put_vara_double(ncid, tg1_id, start3, cuent3, 
     +        TG1)

C--   writting 2d fields
      retval = nf_put_vara_double(ncid, psg0_id, start2, cuent2, 
     +        PSG0)
      retval = nf_put_vara_double(ncid, psg1_id, start2, cuent2, 
     +        PSG1)

C--   writting 4d fields
      retval = nf_put_vara_double(ncid, trg0_id, start4, cuent4, 
     +        TRG0)
      retval = nf_put_vara_double(ncid, trg1_id, start4, cuent4, 
     +        TRG1)

C--   we close de netcdf file
      retval = nf_close(ncid)

      RETURN

      END SUBROUTINE
      
      SUBROUTINE read_variables(UG0,UG1,VG0,VG1,TG0,TG1,PSG0,PSG1,
     +     TRG0,TRG1,IX,IL,KX,NTR)
      INTEGER, INTENT(IN) :: IX,IL,KX,NTR
      REAL, INTENT(INOUT)  ::  UG0(IX,IL,KX)
      REAL, INTENT(INOUT)  ::  UG1(IX,IL,KX)
      REAL, INTENT(INOUT)  ::  VG0(IX,IL,KX)
      REAL, INTENT(INOUT)  ::  VG1(IX,IL,KX)
      REAL, INTENT(INOUT)  ::  TG0(IX,IL,KX)
      REAL, INTENT(INOUT)  ::  TG1(IX,IL,KX)
      REAL, INTENT(INOUT)  ::  PSG0(IX,IL)
      REAL, INTENT(INOUT)  ::  PSG1(IX,IL)
      REAL, INTENT(INOUT)  ::  TRG0(IX,IL,KX,NTR)
      REAL, INTENT(INOUT)  ::  TRG1(IX,IL,KX,NTR)
      character*(*) FILE_NAME
      parameter (FILE_NAME='ensemble_member.nc')
      integer retval, ncid

      character*(*) UG1_N, UG0_N, VG0_N, VG1_N, TG0_N, TG1_N
      character*(*) PSG0_N, PSG1_N, TRG0_N, TRG1_N
      parameter (UG0_N = 'UG0')
      parameter (UG1_N = 'UG1')
      parameter (VG0_N = 'VG0')
      parameter (VG1_N = 'VG1')
      parameter (TG0_N = 'TG0')
      parameter (TG1_N = 'TG1')
      parameter (PSG0_N = 'PSG0')
      parameter (PSG1_N = 'PSG1')
      parameter (TRG0_N = 'TRG0')
      parameter (TRG1_N = 'TRG1')

      INTEGER :: lon_varid, lat_varid, lev_varid, rec_varid
      INTEGER :: ug0_id,ug1_id, vg0_id, vg1_id, tg0_id, tg1_id
      INTEGER :: psg0_id, psg1_id, trg0_id, trg1_id

      integer start4(4), cuent4(4), start3(3), cuent3(3)
      integer start2(2), cuent2(2)

      integer nlons, nlats, nlevs
      nlevs = KX
      nlats = IL
      nlons = IX

C--   we open thefile

      retval = nf_open(FILE_NAME, nf_nowrite, ncid)

C--   we get varids from the variables
      retval = nf_inq_varid(ncid, UG0_N, ug0_id)
      retval = nf_inq_varid(ncid, UG1_N, ug1_id)
      retval = nf_inq_varid(ncid, VG0_N, vg0_id)
      retval = nf_inq_varid(ncid, VG1_N, vg1_id)
      retval = nf_inq_varid(ncid, TG0_N, tg0_id)
      retval = nf_inq_varid(ncid, TG1_N, tg1_id)
      retval = nf_inq_varid(ncid, PSG0_N, psg0_id)
      retval = nf_inq_varid(ncid, PSG1_N, psg1_id)
      retval = nf_inq_varid(ncid, TRG0_N, trg0_id)
      retval = nf_inq_varid(ncid, TRG1_N, trg1_id)

C--   area for reading 2d fields
      start2(1) = 1
      start2(2) = 1
      cuent2(1) = nlons
      cuent2(2) = nlats

C--   cube for reading 3d fields
      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
      cuent3(1) = nlons
      cuent3(2) = nlats
      cuent3(3) = nlevs

C--   hypercube for reading 4d fields
      start4(1) = 1
      start4(2) = 1
      start4(3) = 1
      start4(4) = 1
      cuent4(1) = nlons
      cuent4(2) = nlats
      cuent4(3) = nlevs
      cuent4(4) = NTR

C--   read 2d variables
      retval = nf_get_vara_double(ncid, psg0_id, start2, cuent2,
     +        PSG0)
      retval = nf_get_vara_double(ncid, psg1_id, start2, cuent2,
     +        PSG1)

C--   read 3d variables
      retval = nf_get_vara_double(ncid, ug0_id, start3, cuent3,
     +        UG0)
      retval = nf_get_vara_double(ncid, ug1_id, start3, cuent3,
     +        UG1)
      retval = nf_get_vara_double(ncid, vg0_id, start3, cuent3,
     +        VG0)
      retval = nf_get_vara_double(ncid, vg1_id, start3, cuent3,
     +        VG1)
      retval = nf_get_vara_double(ncid, tg0_id, start3, cuent3,
     +        TG0)
      retval = nf_get_vara_double(ncid, tg1_id, start3, cuent3,
     +        TG1)

c--  read 4d variables
      retval = nf_get_vara_double(ncid, trg0_id, start4, cuent4,
     +        TRG0)
      retval = nf_get_vara_double(ncid, trg1_id, start4, cuent4,
     +        TRG1)

c--   close the netcdf file
      retval = nf_close(ncid)


      print *,'Hi from read ',retval

      RETURN

      END SUBROUTINE

	
      END MODULE




