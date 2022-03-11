
      SUBROUTINE RESTART (JDAY)
      use read_write_netcdf
C--
C--   SUBROUTINE RESTART (JDAY)
C--
C--   Purpose : read or write a restart file
C--   Input :   JDAY  = 0 : read model variables from a restart file
C--                   > 0 : write model variables  to a restart file
C--                         at selected dates and at the end of run 
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_date.h"
      include "com_tsteps.h"

      include "com_dynvar.h"

      include "com_lflags.h"

      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

      COMPLEX DUMC(MX,NX,3)

      REAL DUMR(IX,IL,3)

      REAL UG(IX,IL,KX), VG(IX,IL,KX), TG(IX,IL,KX),
     *     TRG(IX,IL,KX,NTR), PSG(IX,IL)
      REAL UG1(IX,IL,KX), VG1(IX,IL,KX), TG1(IX,IL,KX),
     *     TRG1(IX,IL,KX,NTR), PSG1(IX,IL)
      REAL*4 UGO(NGP,KX), VGO(NGP,KX), TGO(NGP,KX),
     *     TRGO(NGP,KX,NTR), PSGO(NGP)

      CHARACTER (len=3) cexp        ! experiment identifier

C      print *,'BLOCKHOURS', BLOCKHOURS

      IF (JDAY.EQ.0) THEN

         IF (IYEAR.EQ.IYEAR0.AND.IMONTH.EQ.IMONT0) THEN

           print*, '**READING THE RESTART FILE ElíasN - NO TR**'


      CALL read_variables(UG,UG1,VG,VG1,TG,TG1,PSG,PSG1,
     +     TRG,TRG1,IX,IL,KX,NTR)

C     Convert back to spectral tendencies
C
      DO K=1,KX

C--   Step(1)

        CALL VDSPEC(UG(1,1,K),VG(1,1,K),VOR(1,1,K,1),
     *              DIV(1,1,K,1),2)

		if (ix.eq.iy*4) then
		  call TRUNCT (VOR(:,:,K,1))
		  call TRUNCT (DIV(:,:,K,1))
		endif

        CALL SPEC(TG(1,1,K),T(1,1,K,1))
c
        DO  ITR=1,NTR
          CALL SPEC(TRG(1,1,K,ITR),TR(1,1,K,1,ITR))
        ENDDO

C--  Step(2)
        CALL VDSPEC(UG1(1,1,K),VG1(1,1,K),VOR(1,1,K,2),
     *              DIV(1,1,K,2),2)

		if (ix.eq.iy*4) then
		  call TRUNCT (VOR(:,:,K,2))
		  call TRUNCT (DIV(:,:,K,2))
		endif

        CALL SPEC(TG1(1,1,K),T(1,1,K,2))
c
        DO  ITR=1,NTR
          CALL SPEC(TRG1(1,1,K,ITR),TR(1,1,K,2,ITR))
        ENDDO
      ENDDO
      CALL SPEC(PSG(1,1),PS(1,1,1))
      CALL SPEC(PSG1(1,1),PS(1,1,2))

         CALL REST_LAND (0)

         CALL REST_SEA (0)


C--   Forward half step


C       STOP 8888






         ENDIF

C     Check for write-up dates

      ELSE IF (JDAY.eq.NDAYTOT) THEN
      
C      ELSE IF ((mod(IMONTH-1,NMONRS).eq.0)) THEN

c--   spectral tendencies are mapped to grid points
      DO  K=1,KX

C--   Step(1)
        CALL GRID(T(1,1,K,1),TG(1,1,K),1)
        DO  ITR=1,NTR
          CALL GRID(TR(1,1,K,1,ITR),TRG(1,1,K,ITR),1)
        ENDDO

        CALL UVSPEC(VOR(1,1,K,1),DIV(1,1,K,1),
     *              DUMC(1,1,1),DUMC(1,1,2))
        CALL GRID(DUMC(1,1,2),VG(1,1,K),2)
        CALL GRID(DUMC(1,1,1),UG(1,1,K),2)
C--   Step(2)

        CALL GRID(T(1,1,K,2),TG1(1,1,K),1)
        DO  ITR=1,NTR
          CALL GRID(TR(1,1,K,2,ITR),TRG1(1,1,K,ITR),1)
        ENDDO

        CALL UVSPEC(VOR(1,1,K,2),DIV(1,1,K,2),
     *              DUMC(1,1,1),DUMC(1,1,2))
        CALL GRID(DUMC(1,1,2),VG1(1,1,K),2)
        CALL GRID(DUMC(1,1,1),UG1(1,1,K),2)

      ENDDO

      CALL GRID(PS(1,1,1),PSG(1,1),1)
      CALL GRID(PS(1,1,2),PSG1(1,1),1)



C      print *, '+wr VOR',VOR(5,5,1,2),VOR(15,15,2,2),VOR(15,19,3,2)
C      print *, '+wr  DIV',DIV(5,5,1,2),DIV(15,15,2,2),DIV(15,19,3,2)
C      print *, '+wr T',T(5,5,1,2),T(15,15,2,2),T(15,19,3,2)
C      print *, '+wr  PS',PS(5,5,2),PS(15,15,2),PS(15,19,2)

      CALL write_variables(UG,UG1,VG,VG1,TG,TG1,PSG,PSG1,
     +     TRG,TRG1,IX,IL,KX,NTR)

C      CALL write_variables(UG,UG,VG,VG,TG,TG,PSG,PSG,
C     +     TRG,TRG,VORG,VORG,DIVG,DIVG,IX,IL,KX,NTR)

C      CALL write_variables_ss(UG,VG,TG,PSG,
C     +     TRG,VORG,DIVG,IX,IL,KX,NTR)
     
      print *,'Finishing writting the restart file'

C       print *, '+rd Computed'
C       print *, '+rd VOR',VOR(5,5,1,2),VOR(15,15,2,2),VOR(15,19,3,2)
C       print *, '+rd  DIV',DIV(5,5,1,2),DIV(15,15,2,2),DIV(15,19,3,2)
C       print *, '+rd T',T(5,5,1,2),T(15,15,2,2),T(15,19,3,2)
C       print *, '+rd  PS',PS(5,5,2),PS(15,15,2),PS(15,19,2)

C         WRITE (10) VOR
C         WRITE (10) DIV
C         WRITE (10) T
C         WRITE (10) PS
C         WRITE (10) TR

         CALL REST_LAND (1)

         CALL REST_SEA (1)

      ENDIF
C--
      RETURN

C--   4. Stop integration if restart file is not found

  200 CONTINUE

      print*, ' No restart dataset for the specified initial date'

      STOP 'invalid restart'

C--
      END

