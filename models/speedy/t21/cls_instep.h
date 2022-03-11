      NMONTS = 12
      NDAYSL = 0
      NSTEPS = 36
      NSTDIA = 36*5
      NSTPPR = 6
      NSTOUT = -1
      IDOUT  = 0
      NMONRS = -1
      ISEASC = 1
      IYEAR0 = 1979
      IMONT0 = 1
      NSTRAD = 3
      NSTRDF = 0
      INDRDF = 1
      ICLAND = 1
      ICSEA  = 0
      ICICE  = 1
      ISSTAN = 1
      ISSTY0 = 1870
      ISST0  = (IYEAR0-ISSTY0)*12+IMONT0
      LPPRES = .true.
      LCO2 = .false.
      
C--   EliasN
C--   Uncomment for t47 resolution 12*11
C      NMONTS = 1 
C      NDAYSL = 0
C      NSTEPS = 72

C      NSTDIA = 72*5
C      NSTPPR = 6
C      NSTOUT = -1
      
      
C--   Uncomment for t63 resolution 12*55
C      NMONTS = 1
C      NDAYSL = 0
C      NSTEPS = 96

C      NSTDIA = 96*5
C      NSTPPR = 6
      
C--   Uncomment for t106 resolution 12*55
C      NMONTS = 1
C      NDAYSL = 0
C      NSTEPS = 192

C      NSTDIA = 192*5
C      NSTPPR = 6
      
      istart = 0
      HOURS = 0
      BLOCKHOURS = 24./FLOAT(NSTEPS)
