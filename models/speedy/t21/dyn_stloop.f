      SUBROUTINE STLOOP (ISTEP)
C--
C--   SUBROUTINE STLOOP (ISTEP)
C--
C--   Purpose: Perform a series of time steps calling 
C--            post-processing/output routines at selected steps
C--   Input/output : ISTEP = time step index
C--   Updated common block : LFLAG2
C-- 

      include "com_tsteps.h"
      include "com_date.h"
      include "com_lflags.h"
 
      iitest=0

      DO J=1,NSTEPS

      rday = float(J)/float(NSTEPS)

        if (iitest.eq.1) print*, 'STLOOP: calling step ', istep

C       Set logical flags

        LRADSW = (MOD(ISTEP,NSTRAD).EQ.1)
        LRANDF = ((ISTEP.LE.NSTRDF).OR.(NSTRDF.LT.0))

C       Perform one leapfrog time step

c        print *, 'HEYYYYYY ',DELT2,ALPH,ROB,WIL

        CALL STEP (2,2,DELT2,ALPH,ROB,WIL)   

C       Do diagnostic, post-processing and I/O tasks 
 
C        CALL DIAGNS (2,ISTEP)

c        IF (MOD(ISTEP,NSTPPR).EQ.0) CALL TMINC

c        IF (NSTOUT.GT.0.AND.MOD(ISTEP,NSTOUT).EQ.0) CALL TMOUT (1)

        ISTEP=ISTEP+1
        HOUR = INT(FLOAT(ISTEP)*BLOCKHOURS)*1

C        print *,'HOUR = , HOURS = ',HOUR,HOURS

        if ((hour.gt.hours).and.(hours.gt.0)) then
        print *, 'Hours reached'
        CALL AGCM_TO_COUPLER (1)
        CALL COUPLER_TO_AGCM (1)
        call RESTART (1)
        STOP 7777
        endif

      ENDDO

      RETURN
      END
  
