C===================================================================
C===================================================================
C
C
C                        UGAS ROUTINES:
C                        --------------
C
C       SUBROUTINES REPRESENTING THE TRANSPORT PROPERTIES
C       -------------------------------------------------
C
C                     OF EQUILIBRIUM AIR
C                     ------------------
C
C              (by  S. Srinivasan and J.C. Tannehill, 
C               Iowa State University, December, 1987)
C
C                 (transcribed and modified by
C              C.P.T. Groth, U.T.I.A.S., July 1989)
C
C    Subroutines:
C
C       1) UGAS1  - Dynamic viscosity for equilibrium air,
C                   returns viscosity given the equilibrium
C                   temperature and density
C
C       2) UGAS2  - Prandtl number for equilibrium air,
C                   returns Prandtl number given the equilibrium
C                   temperature and density
C
C       3) UGAS4  - Coefficient of thermal conductivity for
C                   equilibrium air,
C                   returns thermal conductivity coefficient
C                   given the equilibrium density and specific
C                   internal energy
C
C
C*******************************************************************
C
C
      SUBROUTINE UGAS1(T,RHO,MU,IERROR)
C
C        Subroutine UGAS1 describes the dynamic viscosity of
C     equilibrium air as a function of the temperature and density.
C     This routine uses curve fits (Srinivasan and Tannehill, 1987)
C     that are based on the calculations of Peng and Pindroh (1962)
C     and constructed from Grabau-type transition functions and
C     bicubic polynomials to model the transport properties of air
C     in a piecewise manner.  The curve fits are valid for
C     temperatures up to 15,000 K and density ratios (RHO/RH0) from
C     1.0E-05 to 10.
C
C     Variable description:
C
C     (call statement parameter list)
C
C     T     Equilibrium temperature (K).
C
C     RHO   Equilibrium density (kg/m**3).
C
C     MU    Dynamic visosity (kg/m-s).
C
C     IERROR Error indicator:
C            If value of 0, then no error.
C            If value of -999 is returned by UGAS1 the
C            equilibrium state is outside the range of validity
C            for the curve fits.
C
C
C     Begin subroutine UGAS1.
C
      IMPLICIT REAL(A-H,M,O-Z)
      REAL MU
      DATA RHO0/1.243/
      DATA ARGMAX/85.00/,EXMAX/1.00E37/,EXMIN/1.00E-37/,
     &     POWMAX/33.00/
C
      IERROR=0
C
      X=T/1000.0E00
      Y=LOG10(RHO/RHO0)
      IF ((Y.LT.-5.0E00).OR.(Y.GT.1.0E00).OR.(X.GT.1.5E01)) THEN
         IERROR=-999
      END IF
C
      IF (T.GT.300.0E00) GO TO 5
      MU=1.462E-06*SQRT(T)/(1.0E00+112.0E00/T)
      RETURN
C
    5 IF (T.GT.600.0E00) GO TO 10
      GAS1=4.13906E-01+2.16606E00*X
      GAS2=(1.30718E-05+7.44367E-05*X)*Y
      GAS3=(-5.45043E-02-1.74550E-04*Y-1.15324E-01*X)*X*X
      GAS4=(2.43199E-05-2.14485E-05*X+3.03976E-06*Y)*Y*Y
      Z=GAS1+GAS2+GAS3+GAS4
      GO TO 100
C
   10 IF (T.GT.5500.0E00) GO TO 20
      GAS1=4.8653102E-01+2.1053953E00*X
      GAS2=(-4.4502862E-02+5.0622325E-02*X)*Y
      GAS3=(-2.3327267E-01-1.019074E-02*Y+1.9685295E-02*X)*X*X
      GAS4=(2.7564680E-04+2.7222564E-03*X+8.6649903E-04*Y)*Y*Y
      Z=GAS1+GAS2+GAS3+GAS4
      GO TO 100
C
   20 IF (T.GT.1.05E04) GO TO 40
      IF (Y.GT.-2.5E00) GO TO 30
      GAS1=5.993881E01-1.698837E01*X
      GAS2=(2.113989E01-5.130287E00*X)*Y
      GAS3=(1.742364E00+2.778705E-01*Y-5.201473E-02*X)*X*X
      GAS4=(1.637675E00-2.460623E-01*X+1.671319E-02*Y)*Y*Y
      GAS5=4.438706E02-5.640044E01*X
      GAS6=(7.553438E01-3.177507E00*X)*Y
      GAS7=(2.078559E00-7.210869E-02*Y-1.669889E-02*X)*X*X
      GAS8=(5.507033E00+1.196004E-02*X+2.199922E-01*Y)*Y*Y
      XXX=107.00E00-7.40E00*X+11.50E00*Y-0.41*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 90
C
   30 GAS1=3.53316E00+4.93425E-01*X
      GAS2=(-4.06143E-01+2.10671E-01*X)*Y
      GAS3=(7.34934E-02-2.63952E-02*Y-1.26658E-03*X)*X*X
      GAS4=(-1.30624E-01+2.74790E-02*X-2.79567E-03*Y)*Y*Y
      GAS5=2.00538E01-6.67992E00*X
      GAS6=(1.03098E01-2.48845E00*X)*Y
      GAS7=(7.57575E-01+1.52736E-01*Y-2.91622E-02*X)*X*X
      GAS8=(1.40864E00-1.80169E-01*X+4.70062E-02*Y)*Y*Y
      XXX=63.75E00-7.976E00*X+5.357E-01*Y+8.333E-01*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 90
C
   40 IF (T.GT.13.0E03) GO TO 60
      IF (Y.GT.-4.00E00) GO TO 50
      GAS1=3.24885E02-5.46359E01*X
      GAS2=(2.86103E01+1.11967E00*X)*Y
      GAS3=(4.10141E00-2.38424E-02*Y-1.08784E-01*X)*X*X
      GAS4=(4.44362E00+2.24907E-01*X+4.59308E-01*Y)*Y*Y
      GAS5=-4.50893E02+3.61004E01*X
      GAS6=(-1.46489E02+1.39297E01*X)*Y
      GAS7=(3.63567E-01-1.33686E-01*Y-2.80207E-02*X)*X*X
      GAS8=(-8.36798E00+8.40047E-01*X+1.77782E-01*Y)*Y*Y
      XXX=2.447E02-1.874E01*X+4.856E01*Y-3.723E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.00-GAS9)
      GO TO 100
C
   50 IF (Y.GT.-2.5E00) GO TO 80
      GAS1=4.74364E01-2.52946E00*X
      GAS2=(-3.40953E01+4.33761E00*X)*Y
      GAS3=(4.19920E-02-1.12842E-01*Y+9.09826E-05*X)*X*X
      GAS4=(-6.95452E00+3.95061E-01*X-3.33072E-01*Y)*Y*Y
      GAS5=-3.45758E02+6.54812E01*X
      GAS6=(-3.77086E01+4.97501E00*X)*Y
      GAS7=(-4.17677E00-1.59916E-01*Y+9.0358E-02*X)*X*X
      GAS8=(2.01908E00-8.41274E-02*X+2.61401E-01*Y)*Y*Y
      XXX=-197.0E00+14.6E00*X-41.2E00*Y+2.85E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 90
C
   60 IF (Y.GT.-4.00E00) GO TO 70
      GAS1=4.53184E02-5.27482E01*X
      GAS2=(1.29609E02-9.91921E00*X)*Y
      GAS3=(2.07504E00+1.90185E-01*Y-2.76109E-02*X)*X*X
      GAS4=(1.26755E01-4.88955E-01*X+4.09428E-01*Y)*Y*Y
      GAS5=-1.15162E02-4.02569E00*X
      GAS6=(-8.84578E01-2.17376E01*X)*Y
      GAS7=(-3.8511E00-9.30247E-01*Y-2.00626E-02*X)*X*X
      GAS8=(-4.78794E01-4.61508E00*X-7.34409E00*Y)*Y*Y
      XXX=76.82-2.29*X+15.08*Y-0.4475*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.00-GAS9)
      GO TO 100
C
   70 IF (Y.GT.-2.50E00) GO TO 80
      GAS1=4.52289E02-5.50932E01*X
      GAS2=(1.18987E02-9.05427E00*X)*Y
      GAS3=(2.34316E00+1.857E-01*Y-3.40138E-02*X)*X*X
      GAS4=(1.13829E01-3.91996E-01*X+4.075E-01*Y)*Y*Y
      GAS5=-7.93738E02+1.27668E02*X
      GAS6=(-1.76755E02+2.93067E01*X)*Y
      GAS7=(-6.18382E00-7.86123E-01*Y+9.67472E-02*X)*X*X
      GAS8=(2.67668E01+1.08936E00*X+7.20878E00*Y)*Y*Y
      XXX=-63.33E00+3.33E00*X-16.67E00*Y+0.667E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 90
C
   80 GAS1=5.05519E02-6.60139E01*X
      GAS2=(1.20621E02-9.43025E00*X)*Y
      GAS3=(3.07622E00+1.93366E-01*Y-5.10125E-02*X)*X*X
      GAS4=(1.10339E01-4.07027E-01*X+3.59506E-01*Y)*Y*Y
      GAS5=-5.14505E02+6.91016E01*X
      GAS6=(-1.15237E02+8.06602E00*X)*Y
      GAS7=(-3.12671E00-1.13473E-01*Y+4.8532E-02*X)*X*X
      GAS8=(-8.75453E00+1.69286E-01*X-2.57493E-01*Y)*Y*Y
      XXX=-156.10E00+9.58E00*X-32.3E00*Y+1.64E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
C
   90 Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.00+GAS9)
  100 MU=1.058E-06*16.5273E00*Z
C
C     End subroutine UGAS1.
C
      RETURN
      END
C
C
C*******************************************************************
C
C
      SUBROUTINE UGAS2(T,RHO,PR,IERROR)
C
C        Subroutine UGAS2 describes the Prandtl number of
C     equilibrium air as a function of the temperature and density.
C     This routine uses curve fits (Srinivasan and Tannehill, 1987)
C     that are based on the calculations of Peng and Pindroh (1962)
C     and constructed from Grabau-type transition functions and
C     bicubic polynomials to model the transport properties of air
C     in a piecewise manner.  The curve fits are valid for
C     temperatures up to 15,000 K and density ratios (RHO/RH0) from
C     1.0E-05 to 10.
C
C     Variable description:
C
C     (call statement parameter list)
C
C     T     Equilibrium temperature (K).
C
C     RHO   Equilibrium density (kg/m**3).
C
C     PR    Prandtl number.
C
C     IERROR Error indicator:
C            If value of 0, then no error.
C            If value of -999 is returned by UGAS2 the
C            equilibrium state is outside the range of validity
C            for the curve fits.
C
C
C     Begin subroutine UGAS2.
C
      IMPLICIT REAL(A-H,M,O-Z)
      DATA RHO0/1.243/
      DATA ARGMAX/85.00/,EXMAX/1.00E37/,EXMIN/1.00E-37/,
     &     POWMAX/33.00/
C
      IERROR=0
C
      X=T/1000.0E00
      Y=LOG10(RHO/RHO0)
      IF ((Y.LT.-5.0E00).OR.(Y.GT.1.0E00).OR.(X.GT.1.5E01)) THEN
         IERROR=-999
      END IF
C
      IF (T.GT.500.0E00) GO TO 10
      GAS1=7.16321E-01+1.1135E00*X
      GAS2=(5.58243E-06-7.16815E-05*X)*Y
      GAS3=(-7.72911E00+2.25827E-04*Y+1.44166E01*X)*X*X
      GAS4=(-1.47156E-07-2.28926E-07*X-2.88338E-08*Y)*Y*Y
      GAS5=-1.4099E-01-3.35055E-01*X
      GAS6=(-2.55975E-05+1.5853E-04*X)*Y
      GAS7=(6.09194E00-3.18345E-04*Y-1.32747E01*X)*X*X
      GAS8=(1.3742E-06-1.29479E-06*X+1.48302E-07*Y)*Y*Y
      XXX=8.636-3.03E01*X
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   10 IF (T.GT.2.0E03) GO TO 20
      GAS1=6.766E-01+5.33391E-02*X
      GAS2=(-2.01021E-02+4.04905E-03*X)*X*X
      Z=GAS1+GAS2
      GO TO 110
C
   20 IF (T.GT.4.0E03) GO TO 30
      GAS1=5.35204E-01+1.64262E-01*X
      GAS2=(-6.72637E-02+3.42314E-02*X)*Y
      GAS3=(-3.88497E-02-3.16248E-03*Y+3.05280E-03*X)*X*X
      GAS4=(-7.81832E-03+1.84389E-03*X-3.46855E-04*Y)*Y*Y
      Z=GAS1+GAS2+GAS3+GAS4
      GO TO 110
C
   30 IF (T.GT.6.5E03) GO TO 40
      GAS1=-2.39283E00+1.28399E00*X
      GAS2=(-7.675E-01+1.89502E-01*X)*Y
      GAS3=(-1.79581E-01-1.20286E-02*Y+8.30322E-03*X)*X*X
      GAS4=(-7.09301E-02+7.19471E-03*X-2.78371E-03*Y)*Y*Y
      GAS5=3.06018E00-1.20461E00*X
      GAS6=(6.77677E-01-1.43868E-01*X)*Y
      GAS7=(1.62407E-01+7.85746E-03*Y-7.39086E-03*X)*X*X
      GAS8=(4.14157E-02-8.3635E-04*X-3.70369E-04*Y)*Y*Y
      XXX=-26.39E00+2.969E00*X-5.042E00*Y-0.112E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   40 IF (T.GE.9.40E03) GO TO 50
      GAS1=6.13473E00-1.54169E00*X
      GAS2=(1.08128E00-2.04154E-01*X)*Y
      GAS3=(1.43737E-01+9.91640E-03*Y-4.54467E-03*X)*X*X
      GAS4=(6.19987E-02-5.05808E-03*X+1.56791E-03*Y)*Y*Y
      GAS5=-5.44445E00+1.58459E00*X
      GAS6=(-1.10792E00+2.13203E-01*X)*Y
      GAS7=(-1.51000E-01-1.00257E-02*Y+4.72964E-03*X)*X*X
      GAS8=(-7.80793E-02+7.29918E-03*X-2.29357E-03*Y)*Y*Y
      XXX=13.39E00-4.258E00*X+2.298E00*Y-1.233E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   50 IF (T.GT.11.5E03) GO TO 70
      IF (Y.GT.-2.5E00) GO TO 60
      GAS1=-3.24776E01+8.72772E00*X
      GAS2=(-2.58872E00+3.59002E-01*X)*Y
      GAS3=(-7.61542E-01+2.10923E-02*Y+2.67953E-02*X)*X*X
      GAS4=(-1.9321E-01+8.94685E-02*X+5.64303E-02*Y)*Y*Y
      GAS5=3.99935E01-9.68334E00*X
      GAS6=(6.78337E00-9.07345E-01*X)*Y
      GAS7=(7.87932E-01-3.99108E-03*Y-2.64764E-02*X)*X*X
      GAS8=(6.97742E-01-1.25709E-01*X-4.08833E-02*Y)*Y*Y
      XXX=105.8E00-11.67E00*X+31.67E00*Y-3.33E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   60 GAS1=-2.80755E01+6.80406E00*X
      GAS2=(-2.63243E00+4.0185E-01*X)*Y
      GAS3=(-5.45283E-01-1.63614E-02*Y+1.45424E-02*X)*X*X
      GAS4=(2.12026E-02-3.62386E-03*X+3.50018E-03*Y)*Y*Y
      GAS5=2.82604E01-6.62279E00*X
      GAS6=(2.06694E00-2.89135E-01*X)*Y
      GAS7=(5.26582E-01+1.13732E-02*Y-1.40944E-02*X)*X*X
      GAS8=(-1.31445E-01+1.50468E-02*X-9.13033E-03*Y)*Y*Y
      XXX=-35.41E00+2.148E00*X-1.481E00*Y-0.3704E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   70 IF (Y.GT.-2.5E00) GO TO 90
      IF (T.GT.13.5E03) GO TO 80
      GAS1=6.08811E01-9.88231E00*X
      GAS2=(9.51872E00-9.95583E-01*X)*Y
      GAS3=(5.48699E-01+2.67619E-02*Y-1.03794E-02*X)*X*X
      GAS4=(5.0199E-01-2.55257E-02*X+8.45834E-03*Y)*Y*Y
      GAS5=-7.1667E01+1.2239E01*X
      GAS6=(-1.09959E01+1.23865E00*X)*Y
      GAS7=(-7.0869E-01-3.61494E-02*Y+1.38445E-02*X)*X*X
      GAS8=(-3.64368E-01+1.89712E-02*X+4.04684E-03*Y)*Y*Y
      XXX=-71.89E00+2.248E00*X+0.4746E00*Y-0.9469E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   80 GAS1=2.99485E01-4.12112E00*X
      GAS2=(5.92879E00-5.27093E-01*X)*Y
      GAS3=(1.93397E-01+1.19371E-02*Y-3.08939E-03*X)*X*X
      GAS4=(4.09472E-01-1.78772E-02*X+9.49505E-03*Y)*Y*Y
      GAS5=-2.66557E01+3.05342E00*X
      GAS6=(-9.53775E00+8.98359E-01*X)*Y
      GAS7=(-9.53141E-02-1.98247E-02*Y+4.69853E-04*X)*X*X
      GAS8=(-7.62232E-01+4.34126E-02*X-5.10053E-03*Y)*Y*Y
      XXX=-540.2E00+34.3E00*X-146.4E00*Y+9.148E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   90 GAS1=-3.18666E00+8.08818E-01*X
      GAS2=(-4.00164E-01+3.59959E-02*X)*Y
      GAS3=(-6.06519E-02-1.04205E-03*Y+1.45243E-03*X)*X*X
      GAS4=(1.6658E-02-4.36487E-03*X-1.86593E-03*Y)*Y*Y
      GAS5=2.68501E00-4.32123E-01*X
      GAS6=(1.36103E-01+2.5886E-02*X)*Y
      GAS7=(2.32842E-02-1.82391E-03*Y-4.09433E-04*X)*X*X
      GAS8=(-5.37705E-02+1.0074E-02*X+2.05852E-03*Y)*Y*Y
      XXX=-31.16E00+1.633E00*X+2.395E00*Y-0.8707E00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
C
  100 Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.00+GAS9)
  110 PR=Z
C
C     End subroutine UGAS2.
C
      RETURN
      END
C
C
C*******************************************************************
C
C
      SUBROUTINE UGAS4(E,RHO,K,IERROR)
C
C        Subroutine UGAS4 describes the coefficient of thermal
C     conductivity of equilibrium air as a function of the internal
C     energy and density.  This routine uses curve fits (Srinivasan
C     and Tannehill, 1987) that are based on the calculations of
C     Peng and Pindroh (1962) and constructed from Grabau-type
C     transition functions and bicubic polynomials to model the
C     transport properties of air in a piecewise manner.  The curve
C     fits are valid for temperatures up to 15,000 K and density
C     ratios (RHO/RH0) from 1.0E-05 to 10.
C
C     Variable description:
C
C     (call statement parameter list)
C
C     E     Equilibrium specific internal energy (J/kg).
C
C     RHO   Equilibrium density (kg/m**3).
C
C     K     Coefficient of thermal conductivity (J/K-m-s).
C
C     IERROR Error indicator:
C            If value of 0, then no error.
C            If value of -999 is returned by UGAS4 the
C            equilibrium state is outside the range of validity
C            for the curve fits.
C
C
C     Begin subroutine UGAS4.
C
      IMPLICIT REAL(A-H,M,O-Z)
      REAL K
      DATA RHO0/1.243E00/,E0/78408.4E00/
      DATA ARGMAX/85.00/,EXMAX/1.00E37/,EXMIN/1.00E-37/,
     &     POWMAX/33.00/
C
      IERROR=0
C
      Z=LOG10(E/E0)
      Y=LOG10(RHO/RHO0)
      IF (Y.LT.-5.0E00.OR.Y.GT.1.0E00) THEN
         IERROR=-999
      END IF
C
      IF (Z.GT.0.44E00) GO TO 5
      T=0.4E00*E/287.06E00
      K=1.994E-03*SQRT(T)/(1.0E00+112.0E00/T)
      RETURN
C
    5 IF (Z.GT.0.65E00) GO TO 10
      GAS1=1.8100369E-01+4.8126802E00*Z
      GAS2=(-2.7231116E-02+1.2691337E-01*Z)*Y
      GAS3=(-8.9913034E00-1.2624085E-01*Y+8.9649105E00*Z)*Z*Z
      GAS4=(-4.7198236E-03+9.2328079E-03*Z-2.9488327E-04*Y)*Y*Y
      F=GAS1+GAS2+GAS3+GAS4
      GO TO 200
C
   10 IF (Y.GT.-1.00E00) GO TO 130
      IF (Y.GT.-3.00E00) GO TO 70
      IF (Z.GT.1.25E00) GO TO 20
      GAS1=-1.05935E04+2.31470E04*Z
      GAS2=(-7.41294E02+1.21724E03*Z)*Y
      GAS3=(-1.67601E04-4.43184E02*Y+4.06631E03*Z)*Z*Z
      GAS4=(1.35105E01+4.94914E00*Z+1.55386E00*Y)*Y*Y
      GAS5=1.06032E04-2.31560E04*Z
      GAS6=(7.46951E02-1.22465E03*Z)*Y
      GAS7=(1.67604E04+4.45919E02*Y-4.06258E03*Z)*Z*Z
      GAS8=(-1.28615E01-5.32398E00*Z-1.52956E00*Y)*Y*Y
      XXX=-4.219E01-4.687E00*Y+2.812E01*Z+3.125E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0-GAS9)
      GO TO 200
C
   20 IF (Z.GT.1.775E00) GO TO 30
      GAS1=3.79375E03-7.40351E03*Z
      GAS2=(3.29698E02-3.55916E02*Z)*Y
      GAS3=(4.77122E03+1.00241E02*Y-1.00740E03*Z)*Z*Z
      GAS4=(1.97061E01-8.42554E00*Z+4.80494E-01*Y)*Y*Y
      GAS5=-4.53603E03+9.05605E03*Z
      GAS6=(-4.95870E02+6.33563E02*Z)*Y
      GAS7=(-5.95317E03-2.05442E02*Y+1.28945E03*Z)*Z*Z
      GAS8=(-2.00087E01+1.18851E01*Z-1.71735E-01*Y)*Y*Y
      XXX=-3.318E01+3.158E-01*Y+1.863E01*Z-1.035E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   30 IF (Z.GT.1.93E00) GO TO 40
      GAS1=2.06651875E05-3.165645E05*Z
      GAS2=(-3.07322021E02+4.57036377E02*Z)*Y
      GAS3=(1.61824937E05-1.55508453E02*Y-2.7603957E04*Z)*Z*Z
      GAS4=(1.92260265E00-2.24788094E00*Z-3.06226015E-01*Y)*Y*Y
      GAS5=-2.06564312E05+3.18191312E05*Z
      GAS6=(2.17542285E03-2.46670776E03*Z)*Y
      GAS7=(-1.63597062E05+7.16753174E02*Y+2.80926367E04*Z)*Z*Z
      GAS8=(3.39526825E01-7.53846645E00*Z+1.91214371E00*Y)*Y*Y
      XXX=-3.924E02-5.206E01*Y+2.054E02*Z+2.679E01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   40 IF (Z.GT.2.60E00) GO TO 50
      GAS1=7.1572625E04-9.2471625E04*Z
      GAS2=(1.9646323E03-2.0280527E03*Z)*Y
      GAS3=(3.9446105E04+4.5673853E02*Y-5.5728672E03*Z)*Z*Z
      GAS4=(-9.2131958E01+1.2724541E01*Z-5.0568476E00*Y)*Y*Y
      GAS5=-3.2910781E04+4.2551211E04*Z
      GAS6=(1.4566331E03-2.2653745E03*Z)*Y
      GAS7=(-1.9476277E04+8.4370288E02*Y+3.2389702E03*Z)*Z*Z
      GAS8=(-1.3324594E02+1.0591533E02*Z+5.8639469E00*Y)*Y*Y
      XXX=4.917E01+2.415E01*Y-2.455E01*Z-1.181E01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   50 IF (Z.GT.2.69E00) GO TO 60
      GAS1=1.145683E06-1.237525E06*Z
      GAS2=(1.4024508E04-9.3467227E03*Z)*Y
      GAS3=(4.4593056E05+1.533074E03*Y-5.3608352E04*Z)*Z*Z
      GAS4=(2.8485107E02-1.0968916E02*Z-1.0955791E00*Y)*Y*Y
      GAS5=-1.752087E06+1.79675E06*Z
      GAS6=(-1.3278737E05+9.8215562E04*Z)*Y
      GAS7=(-6.0791744E05-1.811943E04*Y+6.7709875E04*Z)*Z*Z
      GAS8=(-1.3384084E03+5.2707324E02*Z+2.5904894E00*Y)*Y*Y
      XXX=-1.798E02+7.371E00*Y+6.731E01*Z-3.205E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   60 GAS1=-8.5499625E04+1.1739656E05*Z
      GAS2=(6.4563168E04-3.9551203E04*Z)*Y
      GAS3=(-4.8170254E04+6.0816055E03*Y+6.2052031E03*Z)*Z*Z
      GAS4=(2.3473167E-01+1.8871567E01*Z+4.0757723E00*Y)*Y*Y
      GAS5=5.8546883E04-9.4634875E04*Z
      GAS6=(-6.6513812E04+4.0899945E04*Z)*Y
      GAS7=(4.2127227E04-6.3717305E03*Y-5.7495195E03*Z)*Z*Z
      GAS8=(-1.0260344E00-5.343277E01*Z-1.1017392E01*Y)*Y*Y
      XXX=5.411E00+1.162E01*Y-1.082E00*Z-3.391E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0-GAS9)
      GO TO 200
C
   70 IF (Z.GT.1.29E00) GO TO 80
      GAS1=-1.22493E04+2.41071E04*Z
      GAS2=(-1.61829E03+2.22535E03*Z)*Y
      GAS3=(-1.59261E04-7.53213E02*Y+3.53376E03*Z)*Z*Z
      GAS4=(1.98026E00+5.18483E00*Z+1.47851E00*Y)*Y*Y
      GAS5=1.22486E04-2.41023E04*Z
      GAS6=(1.61810E03-2.22571E03*Z)*Y
      GAS7=(1.59235E04+7.53746E02*Y-3.53168E03*Z)*Z*Z
      GAS8=(-2.15482E00-5.05115E00*Z-1.48795E00*Y)*Y*Y
      XXX=-3.111E01-4.444E00*Y+1.944E01*Z+2.778E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0-GAS9)
      GO TO 200
C
   80 IF (Z.GT.1.85E00) GO TO 90
      GAS1=3.18060E03-6.69664E03*Z
      GAS2=(4.33382E01-2.14649E02*Z)*Y
      GAS3=(4.41377E03+9.41359E01*Y-9.29758E02*Z)*Z*Z
      GAS4=(-3.62190E01+1.15538E01*Z-2.14621E00*Y)*Y*Y
      GAS5=-5.98764E03+1.29243E04*Z
      GAS6=(-2.72261E02+5.42378E02*Z)*Y
      GAS7=(-9.03293E03-2.11787E02*Y+2.07831E03*Z)*Z*Z
      GAS8=(2.74179E01-5.68578E00*Z+1.91217E00*Y)*Y*Y
      XXX=-1.854E01+7.11E00*Y+1.068E01*Z-5.449E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   90 IF (Z.GT.2.0E00) GO TO 100
      GAS1=5.14024E04-7.52733E04*Z
      GAS2=(-3.30889E02+3.11550E02*Z)*Y
      GAS3=(3.66539E04-7.41227E01*Y-5.93015E03*Z)*Z*Z
      GAS4=(-4.84164E01+2.23133E01*Z-9.19118E-01*Y)*Y*Y
      GAS5=-1.80898E05+2.82532E05*Z
      GAS6=(-1.01053E03+9.75576E02*Z)*Y
      GAS7=(-1.47220E05-2.33631E02*Y+2.55940E04*Z)*Z*Z
      GAS8=(3.28681E00-1.76588E00*Z-1.54962E-01*Y)*Y*Y
      XXX=-4.104E01+6.507E01*Y+2.083E01*Z-3.472E01*Z*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  100 IF (Z.GT.2.58E00) GO TO 110
      GAS1=5.1131824E04-6.664875E04*Z
      GAS2=(2.02171E03-1.9306292E03*Z)*Y
      GAS3=(2.8762395E04+4.3353467E02*Y-4.1064609E03*Z)*Z*Z
      GAS4=(-8.4970047E01+1.7925919E01*Z-6.2576542E00*Y)*Y*Y
      GAS5=-6.2768156E04+8.6015875E04*Z
      GAS6=(-1.0002036E03+6.2537280E02*Z)*Y
      GAS7=(-3.957827E04-3.8467377E01*Y+6.12953E03*Z)*Z*Z
      GAS8=(-1.0591702E02+7.636142E01*Z+5.938859E00*Y)*Y*Y
      XXX=-3.901E00+2.418E01*Y+1.374E00*Z-1.145E01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  110 IF (Z.GT.2.73E00) GO TO 120
      GAS1=1.0088046E06-1.086321E06*Z
      GAS2=(1.3844801E04-9.7268516E03*Z)*Y
      GAS3=(3.8985325E05+1.7091665E03*Y-4.6621066E04*Z)*Z*Z
      GAS4=(1.4840726E02-5.2645004E01*Z-1.5477133E-01*Y)*Y*Y
      GAS5=-1.073351E06+1.14571E06*Z
      GAS6=(-1.9343957E04+1.3366211E04*Z)*Y
      GAS7=(-4.0670987E05-2.2955198E03*Y+4.7999871E04*Z)*Z*Z
      GAS8=(-4.1016724E02+1.4994148E02*Z-1.9779787E00*Y)*Y*Y
      XXX=-1.026E02+6.302E01*Y+3.819E01*Z-2.431E01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  120 GAS1=-9.6638500E04+1.3206488E04*Z
      GAS2=(-4.7458105E04+2.3596875E04*Z)*Y
      GAS3=(1.8602773E04-2.306802E03*Y-4.0413552E03*Z)*Z*Z
      GAS4=(-5.3564258E03+2.2433904E03*Z+2.5188145E02*Y)*Y*Y
      GAS5=1.0962581E05-2.990116E04*Z
      GAS6=(4.7883496E04-2.3785383E04*Z)*Y
      GAS7=(-1.1753969E04+2.2905522E03*Y+3.1304399E03*Z)*Z*Z
      GAS8=(5.473418E03-2.3208018E03*Z-2.6570068E02*Y)*Y*Y
      XXX=-3.107E01+1.082E01*Y+1.047E01*Z-3.047E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0-GAS9)
      GO TO 200
C
  130 IF (Z.GT.1.40E00) GO TO 140
      GAS1=-1.58386E03+3.49223E03*Z
      GAS2=(-8.39834E02+1.09565E03*Z)*Y
      GAS3=(-2.56175E03-3.56197E02*Y+6.25145E02*Z)*Z*Z
      GAS4=(-1.22407E01+7.65634E00*Z+2.58235E-01*Y)*Y*Y
      GAS5=1.58025E03-3.47664E03*Z
      GAS6=(8.39588E02-1.09490E03*Z)*Y
      GAS7=(2.54682E03+3.55674E02*Y-6.18504E02*Z)*Z*Z
      GAS8=(1.20843E01-7.44857E00*Z-2.91202E-01*Y)*Y*Y
      XXX=-2.171E01-4.342E00*Y+1.316E01*Z+2.632E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0-GAS9)
      GO TO 200
C
  140 IF (Z.GT.1.91E00) GO TO 150
      GAS1=7.89255E02-1.91743E03*Z
      GAS2=(3.59227E02-4.44070E02*Z)*Y
      GAS3=(1.39463E03+1.34083E02*Y-3.13446E02*Z)*Z*Z
      GAS4=(1.90681E01-1.09285E01*Z+4.24933E-02*Y)*Y*Y
      GAS5=-1.31401E03+3.13134E03*Z
      GAS6=(-5.18755E02+6.80268E02*Z)*Y
      GAS7=(-2.32493E03-2.21393E02*Y+5.52563E02*Z)*Z*Z
      GAS8=(-3.32001E01+2.11819E01*Z-4.75163E-01*Y)*Y*Y
      XXX=-5.025E01-8.412E00*Y+2.982E01*Z+3.509E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  150 IF (Z.GT.2.05E00) GO TO 160
      GAS1=3.58691E04-5.16852E04*Z
      GAS2=(-6.30189E02+6.63314E02*Z)*Y
      GAS3=(2.47471E04-1.73538E02*Y-3.93167E03*Z)*Z*Z
      GAS4=(-4.23871E01+2.08048E01*Z-1.05512E00*Y)*Y*Y
      GAS5=-1.10522E05+1.67591E05*Z
      GAS6=(4.61877E03-4.94930E03*Z)*Y
      GAS7=(-8.46558E04+1.32441E03*Y+1.42438E04*Z)*Z*Z
      GAS8=(2.25065E01-1.10316E01*Z+9.62887E-01*Y)*Y*Y
      XXX=-1.681E02+7.063E01*Y+8.75E01*Z-3.75E01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  160 IF (Z.GT.2.57E00) GO TO 170
      GAS1=3.1899562E04-4.2186664E04*Z
      GAS2=(2.3055603E03-1.9897017E03*Z)*Y
      GAS3=(1.849998E04+4.2561816E02*Y-2.6808696E03*Z)*Z*Z
      GAS4=(-1.6195114E01+5.8640623E00*Z-3.6172504E00*Y)*Y*Y
      GAS5=-5.7594039E04+7.9328437E04*Z
      GAS6=(-1.9275989E03+1.6730544E03*Z)*Y
      GAS7=(-3.6473008E04-3.6100732E02*Y+5.597543E03*Z)*Z*Z
      GAS8=(-7.920808E01+4.0542084E01*Z+2.1495867E00*Y)*Y*Y
      XXX=-5.733E01+2.088E01*Y+2.592E01*Z-9.793E00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  170 IF (Z.GT.2.75E00) GO TO 180
      GAS1=7.0838087E05-7.5619919E05*Z
      GAS2=(3.9503091E03-2.7381802E03*Z)*Y
      GAS3=(2.6888181E05+4.7728687E02*Y-3.183816E04*Z)*Z*Z
      GAS4=(-1.2532251E02+4.7734787E01*Z-4.0148029E00*Y)*Y*Y
      GAS5=-2.5216325E05+2.1727769E05*Z
      GAS6=(9.2882383E03-7.780918E03*Z)*Y
      GAS7=(-5.6539297E04+1.6120212E03*Y+3.9419248E03*Z)*Z*Z
      GAS8=(1.8537296E02-7.1010757E01*Z+1.1307096E00*Y)*Y*Y
      XXX=-1.786E02+2.18E-01*Y+6.714E01*Z-4.739E-01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  180 GAS1=3.1855037E05-3.3041156E05*Z
      GAS2=(2.2983352E04-1.6623461E04*Z)*Y
      GAS3=(1.13848E05+3.0098223E03*Y-1.3020133E04*Z)*Z*Z
      GAS4=(-1.8599039E02+6.9840683E01*Z-7.7371645E00*Y)*Y*Y
      F=GAS1+GAS2+GAS3+GAS4
      GO TO 200
C
  190 F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0+GAS9)
  200 K=1.87915E-02*F
C
C     End subroutine UGAS4.
C
      RETURN
      END
