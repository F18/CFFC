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
C     1.0D-05 to 10.
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
      IMPLICIT REAL*8(A-H,M,O-Z)
      REAL*8 MU
      DATA RHO0/1.243D00/
      DATA ARGMAX/690.00D00/,EXMAX/1.00D300/,EXMIN/1.00D-300/,
     &     POWMAX/299.00D00/
C
      IERROR=0
C
      X=T/1000.00D00
      Y=LOG10(RHO/RHO0)
      IF ((Y.LT.-5.0D00).OR.(Y.GT.1.0D00).OR.(X.GT.1.5D01)) THEN
         IERROR=-999
      END IF
C
      IF (T.GT.300.0D00) GO TO 5
      MU=1.462D-06*SQRT(T)/(1.0D00+112.0D00/T)
      RETURN
C
    5 IF (T.GT.600.0D00) GO TO 10
      GAS1=4.13906D-01+2.16606D00*X
      GAS2=(1.30718D-05+7.44367D-05*X)*Y
      GAS3=(-5.45043D-02-1.74550D-04*Y-1.15324D-01*X)*X*X
      GAS4=(2.43199D-05-2.14485D-05*X+3.03976D-06*Y)*Y*Y
      Z=GAS1+GAS2+GAS3+GAS4
      GO TO 100
C
   10 IF (T.GT.5500.0D00) GO TO 20
      GAS1=4.8653102D-01+2.1053953D00*X
      GAS2=(-4.4502862D-02+5.0622325D-02*X)*Y
      GAS3=(-2.3327267D-01-1.019074D-02*Y+1.9685295D-02*X)*X*X
      GAS4=(2.7564680D-04+2.7222564D-03*X+8.6649903D-04*Y)*Y*Y
      Z=GAS1+GAS2+GAS3+GAS4
      GO TO 100
C
   20 IF (T.GT.1.05D04) GO TO 40
      IF (Y.GT.-2.5D00) GO TO 30
      GAS1=5.993881D01-1.698837D01*X
      GAS2=(2.113989D01-5.130287D00*X)*Y
      GAS3=(1.742364D00+2.778705D-01*Y-5.201473D-02*X)*X*X
      GAS4=(1.637675D00-2.460623D-01*X+1.671319D-02*Y)*Y*Y
      GAS5=4.438706D02-5.640044D01*X
      GAS6=(7.553438D01-3.177507D00*X)*Y
      GAS7=(2.078559D00-7.210869D-02*Y-1.669889D-02*X)*X*X
      GAS8=(5.507033D00+1.196004D-02*X+2.199922D-01*Y)*Y*Y
      XXX=107.00D00-7.40D00*X+11.50D00*Y-0.41*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 90
C
   30 GAS1=3.53316D00+4.93425D-01*X
      GAS2=(-4.06143D-01+2.10671D-01*X)*Y
      GAS3=(7.34934D-02-2.63952D-02*Y-1.26658D-03*X)*X*X
      GAS4=(-1.30624D-01+2.74790D-02*X-2.79567D-03*Y)*Y*Y
      GAS5=2.00538D01-6.67992D00*X
      GAS6=(1.03098D01-2.48845D00*X)*Y
      GAS7=(7.57575D-01+1.52736D-01*Y-2.91622D-02*X)*X*X
      GAS8=(1.40864D00-1.80169D-01*X+4.70062D-02*Y)*Y*Y
      XXX=63.75D00-7.976D00*X+5.357D-01*Y+8.333D-01*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 90
C
   40 IF (T.GT.13.0D03) GO TO 60
      IF (Y.GT.-4.00D00) GO TO 50
      GAS1=3.24885D02-5.46359D01*X
      GAS2=(2.86103D01+1.11967D00*X)*Y
      GAS3=(4.10141D00-2.38424D-02*Y-1.08784D-01*X)*X*X
      GAS4=(4.44362D00+2.24907D-01*X+4.59308D-01*Y)*Y*Y
      GAS5=-4.50893D02+3.61004D01*X
      GAS6=(-1.46489D02+1.39297D01*X)*Y
      GAS7=(3.63567D-01-1.33686D-01*Y-2.80207D-02*X)*X*X
      GAS8=(-8.36798D00+8.40047D-01*X+1.77782D-01*Y)*Y*Y
      XXX=2.447D02-1.874D01*X+4.856D01*Y-3.723D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/
     &  (1.00D00-GAS9)
      GO TO 100
C
   50 IF (Y.GT.-2.5D00) GO TO 80
      GAS1=4.74364D01-2.52946D00*X
      GAS2=(-3.40953D01+4.33761D00*X)*Y
      GAS3=(4.19920D-02-1.12842D-01*Y+9.09826D-05*X)*X*X
      GAS4=(-6.95452D00+3.95061D-01*X-3.33072D-01*Y)*Y*Y
      GAS5=-3.45758D02+6.54812D01*X
      GAS6=(-3.77086D01+4.97501D00*X)*Y
      GAS7=(-4.17677D00-1.59916D-01*Y+9.0358D-02*X)*X*X
      GAS8=(2.01908D00-8.41274D-02*X+2.61401D-01*Y)*Y*Y
      XXX=-197.0D00+14.6D00*X-41.2D00*Y+2.85D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 90
C
   60 IF (Y.GT.-4.00D00) GO TO 70
      GAS1=4.53184D02-5.27482D01*X
      GAS2=(1.29609D02-9.91921D00*X)*Y
      GAS3=(2.07504D00+1.90185D-01*Y-2.76109D-02*X)*X*X
      GAS4=(1.26755D01-4.88955D-01*X+4.09428D-01*Y)*Y*Y
      GAS5=-1.15162D02-4.02569D00*X
      GAS6=(-8.84578D01-2.17376D01*X)*Y
      GAS7=(-3.8511D00-9.30247D-01*Y-2.00626D-02*X)*X*X
      GAS8=(-4.78794D01-4.61508D00*X-7.34409D00*Y)*Y*Y
      XXX=76.82-2.29*X+15.08*Y-0.4475*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/
     &  (1.00D00-GAS9)
      GO TO 100
C
   70 IF (Y.GT.-2.50D00) GO TO 80
      GAS1=4.52289D02-5.50932D01*X
      GAS2=(1.18987D02-9.05427D00*X)*Y
      GAS3=(2.34316D00+1.857D-01*Y-3.40138D-02*X)*X*X
      GAS4=(1.13829D01-3.91996D-01*X+4.075D-01*Y)*Y*Y
      GAS5=-7.93738D02+1.27668D02*X
      GAS6=(-1.76755D02+2.93067D01*X)*Y
      GAS7=(-6.18382D00-7.86123D-01*Y+9.67472D-02*X)*X*X
      GAS8=(2.67668D01+1.08936D00*X+7.20878D00*Y)*Y*Y
      XXX=-63.33D00+3.33D00*X-16.67D00*Y+0.667D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 90
C
   80 GAS1=5.05519D02-6.60139D01*X
      GAS2=(1.20621D02-9.43025D00*X)*Y
      GAS3=(3.07622D00+1.93366D-01*Y-5.10125D-02*X)*X*X
      GAS4=(1.10339D01-4.07027D-01*X+3.59506D-01*Y)*Y*Y
      GAS5=-5.14505D02+6.91016D01*X
      GAS6=(-1.15237D02+8.06602D00*X)*Y
      GAS7=(-3.12671D00-1.13473D-01*Y+4.8532D-02*X)*X*X
      GAS8=(-8.75453D00+1.69286D-01*X-2.57493D-01*Y)*Y*Y
      XXX=-156.10D00+9.58D00*X-32.3D00*Y+1.64D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
C
   90 Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/
     &  (1.00D00+GAS9)
  100 MU=1.058D-06*16.5273D00*Z
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
C     1.0D-05 to 10.
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
      IMPLICIT REAL*8(A-H,M,O-Z)
      DATA RHO0/1.243D00/
      DATA ARGMAX/690.00D00/,EXMAX/1.00D300/,EXMIN/1.00D-300/,
     &     POWMAX/299.00D00/
C
      IERROR=0
C
      X=T/1000.00D00
      Y=LOG10(RHO/RHO0)
      IF ((Y.LT.-5.0D00).OR.(Y.GT.1.0D00).OR.(X.GT.1.5D01)) THEN
         IERROR=-999
      END IF
C
      IF (T.GT.500.0D00) GO TO 10
      GAS1=7.16321D-01+1.1135D00*X
      GAS2=(5.58243D-06-7.16815D-05*X)*Y
      GAS3=(-7.72911D00+2.25827D-04*Y+1.44166D01*X)*X*X
      GAS4=(-1.47156D-07-2.28926D-07*X-2.88338D-08*Y)*Y*Y
      GAS5=-1.4099D-01-3.35055D-01*X
      GAS6=(-2.55975D-05+1.5853D-04*X)*Y
      GAS7=(6.09194D00-3.18345D-04*Y-1.32747D01*X)*X*X
      GAS8=(1.3742D-06-1.29479D-06*X+1.48302D-07*Y)*Y*Y
      XXX=8.636-3.03D01*X
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   10 IF (T.GT.2.0D03) GO TO 20
      GAS1=6.766D-01+5.33391D-02*X
      GAS2=(-2.01021D-02+4.04905D-03*X)*X*X
      Z=GAS1+GAS2
      GO TO 110
C
   20 IF (T.GT.4.0D03) GO TO 30
      GAS1=5.35204D-01+1.64262D-01*X
      GAS2=(-6.72637D-02+3.42314D-02*X)*Y
      GAS3=(-3.88497D-02-3.16248D-03*Y+3.05280D-03*X)*X*X
      GAS4=(-7.81832D-03+1.84389D-03*X-3.46855D-04*Y)*Y*Y
      Z=GAS1+GAS2+GAS3+GAS4
      GO TO 110
C
   30 IF (T.GT.6.5D03) GO TO 40
      GAS1=-2.39283D00+1.28399D00*X
      GAS2=(-7.675D-01+1.89502D-01*X)*Y
      GAS3=(-1.79581D-01-1.20286D-02*Y+8.30322D-03*X)*X*X
      GAS4=(-7.09301D-02+7.19471D-03*X-2.78371D-03*Y)*Y*Y
      GAS5=3.06018D00-1.20461D00*X
      GAS6=(6.77677D-01-1.43868D-01*X)*Y
      GAS7=(1.62407D-01+7.85746D-03*Y-7.39086D-03*X)*X*X
      GAS8=(4.14157D-02-8.3635D-04*X-3.70369D-04*Y)*Y*Y
      XXX=-26.39D00+2.969D00*X-5.042D00*Y-0.112D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   40 IF (T.GE.9.40D03) GO TO 50
      GAS1=6.13473D00-1.54169D00*X
      GAS2=(1.08128D00-2.04154D-01*X)*Y
      GAS3=(1.43737D-01+9.91640D-03*Y-4.54467D-03*X)*X*X
      GAS4=(6.19987D-02-5.05808D-03*X+1.56791D-03*Y)*Y*Y
      GAS5=-5.44445D00+1.58459D00*X
      GAS6=(-1.10792D00+2.13203D-01*X)*Y
      GAS7=(-1.51000D-01-1.00257D-02*Y+4.72964D-03*X)*X*X
      GAS8=(-7.80793D-02+7.29918D-03*X-2.29357D-03*Y)*Y*Y
      XXX=13.39D00-4.258D00*X+2.298D00*Y-1.233D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   50 IF (T.GT.11.5D03) GO TO 70
      IF (Y.GT.-2.5D00) GO TO 60
      GAS1=-3.24776D01+8.72772D00*X
      GAS2=(-2.58872D00+3.59002D-01*X)*Y
      GAS3=(-7.61542D-01+2.10923D-02*Y+2.67953D-02*X)*X*X
      GAS4=(-1.9321D-01+8.94685D-02*X+5.64303D-02*Y)*Y*Y
      GAS5=3.99935D01-9.68334D00*X
      GAS6=(6.78337D00-9.07345D-01*X)*Y
      GAS7=(7.87932D-01-3.99108D-03*Y-2.64764D-02*X)*X*X
      GAS8=(6.97742D-01-1.25709D-01*X-4.08833D-02*Y)*Y*Y
      XXX=105.8D00-11.67D00*X+31.67D00*Y-3.33D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   60 GAS1=-2.80755D01+6.80406D00*X
      GAS2=(-2.63243D00+4.0185D-01*X)*Y
      GAS3=(-5.45283D-01-1.63614D-02*Y+1.45424D-02*X)*X*X
      GAS4=(2.12026D-02-3.62386D-03*X+3.50018D-03*Y)*Y*Y
      GAS5=2.82604D01-6.62279D00*X
      GAS6=(2.06694D00-2.89135D-01*X)*Y
      GAS7=(5.26582D-01+1.13732D-02*Y-1.40944D-02*X)*X*X
      GAS8=(-1.31445D-01+1.50468D-02*X-9.13033D-03*Y)*Y*Y
      XXX=-35.41D00+2.148D00*X-1.481D00*Y-0.3704D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   70 IF (Y.GT.-2.5D00) GO TO 90
      IF (T.GT.13.5D03) GO TO 80
      GAS1=6.08811D01-9.88231D00*X
      GAS2=(9.51872D00-9.95583D-01*X)*Y
      GAS3=(5.48699D-01+2.67619D-02*Y-1.03794D-02*X)*X*X
      GAS4=(5.0199D-01-2.55257D-02*X+8.45834D-03*Y)*Y*Y
      GAS5=-7.1667D01+1.2239D01*X
      GAS6=(-1.09959D01+1.23865D00*X)*Y
      GAS7=(-7.0869D-01-3.61494D-02*Y+1.38445D-02*X)*X*X
      GAS8=(-3.64368D-01+1.89712D-02*X+4.04684D-03*Y)*Y*Y
      XXX=-71.89D00+2.248D00*X+0.4746D00*Y-0.9469D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   80 GAS1=2.99485D01-4.12112D00*X
      GAS2=(5.92879D00-5.27093D-01*X)*Y
      GAS3=(1.93397D-01+1.19371D-02*Y-3.08939D-03*X)*X*X
      GAS4=(4.09472D-01-1.78772D-02*X+9.49505D-03*Y)*Y*Y
      GAS5=-2.66557D01+3.05342D00*X
      GAS6=(-9.53775D00+8.98359D-01*X)*Y
      GAS7=(-9.53141D-02-1.98247D-02*Y+4.69853D-04*X)*X*X
      GAS8=(-7.62232D-01+4.34126D-02*X-5.10053D-03*Y)*Y*Y
      XXX=-540.2D00+34.3D00*X-146.4D00*Y+9.148D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 100
C
   90 GAS1=-3.18666D00+8.08818D-01*X
      GAS2=(-4.00164D-01+3.59959D-02*X)*Y
      GAS3=(-6.06519D-02-1.04205D-03*Y+1.45243D-03*X)*X*X
      GAS4=(1.6658D-02-4.36487D-03*X-1.86593D-03*Y)*Y*Y
      GAS5=2.68501D00-4.32123D-01*X
      GAS6=(1.36103D-01+2.5886D-02*X)*Y
      GAS7=(2.32842D-02-1.82391D-03*Y-4.09433D-04*X)*X*X
      GAS8=(-5.37705D-02+1.0074D-02*X+2.05852D-03*Y)*Y*Y
      XXX=-31.16D00+1.633D00*X+2.395D00*Y-0.8707D00*X*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
C
  100 Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/
     &  (1.00D00+GAS9)
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
C     ratios (RHO/RH0) from 1.0D-05 to 10.
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
      IMPLICIT REAL*8(A-H,M,O-Z)
      REAL*8 K
      DATA RHO0/1.243D00/,D0/78408.4D00/
      DATA ARGMAX/690.00D00/,EXMAX/1.00D300/,EXMIN/1.00D-300/,
     &     POWMAX/299.00D00/
C
      IERROR=0
C
      Z=LOG10(E/E0)
      Y=LOG10(RHO/RHO0)
      IF (Y.LT.-5.0D00.OR.Y.GT.1.0D00) THEN
         IERROR=-999
      END IF
C
      IF (Z.GT.0.44D00) GO TO 5
      T=0.4D00*E/287.06D00
      K=1.994D-03*SQRT(T)/(1.0D00+112.0D00/T)
      RETURN
C
    5 IF (Z.GT.0.65D00) GO TO 10
      GAS1=1.8100369D-01+4.8126802D00*Z
      GAS2=(-2.7231116D-02+1.2691337D-01*Z)*Y
      GAS3=(-8.9913034D00-1.2624085D-01*Y+8.9649105D00*Z)*Z*Z
      GAS4=(-4.7198236D-03+9.2328079D-03*Z-2.9488327D-04*Y)*Y*Y
      F=GAS1+GAS2+GAS3+GAS4
      GO TO 200
C
   10 IF (Y.GT.-1.00D00) GO TO 130
      IF (Y.GT.-3.00D00) GO TO 70
      IF (Z.GT.1.25D00) GO TO 20
      GAS1=-1.05935D04+2.31470D04*Z
      GAS2=(-7.41294D02+1.21724D03*Z)*Y
      GAS3=(-1.67601D04-4.43184D02*Y+4.06631D03*Z)*Z*Z
      GAS4=(1.35105D01+4.94914D00*Z+1.55386D00*Y)*Y*Y
      GAS5=1.06032D04-2.31560D04*Z
      GAS6=(7.46951D02-1.22465D03*Z)*Y
      GAS7=(1.67604D04+4.45919D02*Y-4.06258D03*Z)*Z*Z
      GAS8=(-1.28615D01-5.32398D00*Z-1.52956D00*Y)*Y*Y
      XXX=-4.219D01-4.687D00*Y+2.812D01*Z+3.125D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/
     &  (1.00D00-GAS9)
      GO TO 200
C
   20 IF (Z.GT.1.775D00) GO TO 30
      GAS1=3.79375D03-7.40351D03*Z
      GAS2=(3.29698D02-3.55916D02*Z)*Y
      GAS3=(4.77122D03+1.00241D02*Y-1.00740D03*Z)*Z*Z
      GAS4=(1.97061D01-8.42554D00*Z+4.80494D-01*Y)*Y*Y
      GAS5=-4.53603D03+9.05605D03*Z
      GAS6=(-4.95870D02+6.33563D02*Z)*Y
      GAS7=(-5.95317D03-2.05442D02*Y+1.28945D03*Z)*Z*Z
      GAS8=(-2.00087D01+1.18851D01*Z-1.71735D-01*Y)*Y*Y
      XXX=-3.318D01+3.158D-01*Y+1.863D01*Z-1.035D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   30 IF (Z.GT.1.93D00) GO TO 40
      GAS1=2.06651875D05-3.165645D05*Z
      GAS2=(-3.07322021D02+4.57036377D02*Z)*Y
      GAS3=(1.61824937D05-1.55508453D02*Y-2.7603957D04*Z)*Z*Z
      GAS4=(1.92260265D00-2.24788094D00*Z-3.06226015D-01*Y)*Y*Y
      GAS5=-2.06564312D05+3.18191312D05*Z
      GAS6=(2.17542285D03-2.46670776D03*Z)*Y
      GAS7=(-1.63597062D05+7.16753174D02*Y+2.80926367D04*Z)*Z*Z
      GAS8=(3.39526825D01-7.53846645D00*Z+1.91214371D00*Y)*Y*Y
      XXX=-3.924D02-5.206D01*Y+2.054D02*Z+2.679D01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   40 IF (Z.GT.2.60D00) GO TO 50
      GAS1=7.1572625D04-9.2471625D04*Z
      GAS2=(1.9646323D03-2.0280527D03*Z)*Y
      GAS3=(3.9446105D04+4.5673853D02*Y-5.5728672D03*Z)*Z*Z
      GAS4=(-9.2131958D01+1.2724541D01*Z-5.0568476D00*Y)*Y*Y
      GAS5=-3.2910781D04+4.2551211D04*Z
      GAS6=(1.4566331D03-2.2653745D03*Z)*Y
      GAS7=(-1.9476277D04+8.4370288D02*Y+3.2389702D03*Z)*Z*Z
      GAS8=(-1.3324594D02+1.0591533D02*Z+5.8639469D00*Y)*Y*Y
      XXX=4.917D01+2.415D01*Y-2.455D01*Z-1.181D01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   50 IF (Z.GT.2.69D00) GO TO 60
      GAS1=1.145683D06-1.237525D06*Z
      GAS2=(1.4024508D04-9.3467227D03*Z)*Y
      GAS3=(4.4593056D05+1.533074D03*Y-5.3608352D04*Z)*Z*Z
      GAS4=(2.8485107D02-1.0968916D02*Z-1.0955791D00*Y)*Y*Y
      GAS5=-1.752087D06+1.79675D06*Z
      GAS6=(-1.3278737D05+9.8215562D04*Z)*Y
      GAS7=(-6.0791744D05-1.811943D04*Y+6.7709875D04*Z)*Z*Z
      GAS8=(-1.3384084D03+5.2707324D02*Z+2.5904894D00*Y)*Y*Y
      XXX=-1.798D02+7.371D00*Y+6.731D01*Z-3.205D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   60 GAS1=-8.5499625D04+1.1739656D05*Z
      GAS2=(6.4563168D04-3.9551203D04*Z)*Y
      GAS3=(-4.8170254D04+6.0816055D03*Y+6.2052031D03*Z)*Z*Z
      GAS4=(2.3473167D-01+1.8871567D01*Z+4.0757723D00*Y)*Y*Y
      GAS5=5.8546883D04-9.4634875D04*Z
      GAS6=(-6.6513812D04+4.0899945D04*Z)*Y
      GAS7=(4.2127227D04-6.3717305D03*Y-5.7495195D03*Z)*Z*Z
      GAS8=(-1.0260344D00-5.343277D01*Z-1.1017392D01*Y)*Y*Y
      XXX=5.411D00+1.162D01*Y-1.082D00*Z-3.391D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0-GAS9)
      GO TO 200
C
   70 IF (Z.GT.1.29D00) GO TO 80
      GAS1=-1.22493D04+2.41071D04*Z
      GAS2=(-1.61829D03+2.22535D03*Z)*Y
      GAS3=(-1.59261D04-7.53213D02*Y+3.53376D03*Z)*Z*Z
      GAS4=(1.98026D00+5.18483D00*Z+1.47851D00*Y)*Y*Y
      GAS5=1.22486D04-2.41023D04*Z
      GAS6=(1.61810D03-2.22571D03*Z)*Y
      GAS7=(1.59235D04+7.53746D02*Y-3.53168D03*Z)*Z*Z
      GAS8=(-2.15482D00-5.05115D00*Z-1.48795D00*Y)*Y*Y
      XXX=-3.111D01-4.444D00*Y+1.944D01*Z+2.778D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/
     &  (1.00D00-GAS9)
      GO TO 200
C
   80 IF (Z.GT.1.85D00) GO TO 90
      GAS1=3.18060D03-6.69664D03*Z
      GAS2=(4.33382D01-2.14649D02*Z)*Y
      GAS3=(4.41377D03+9.41359D01*Y-9.29758D02*Z)*Z*Z
      GAS4=(-3.62190D01+1.15538D01*Z-2.14621D00*Y)*Y*Y
      GAS5=-5.98764D03+1.29243D04*Z
      GAS6=(-2.72261D02+5.42378D02*Z)*Y
      GAS7=(-9.03293D03-2.11787D02*Y+2.07831D03*Z)*Z*Z
      GAS8=(2.74179D01-5.68578D00*Z+1.91217D00*Y)*Y*Y
      XXX=-1.854D01+7.11D00*Y+1.068D01*Z-5.449D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
   90 IF (Z.GT.2.0D00) GO TO 100
      GAS1=5.14024D04-7.52733D04*Z
      GAS2=(-3.30889D02+3.11550D02*Z)*Y
      GAS3=(3.66539D04-7.41227D01*Y-5.93015D03*Z)*Z*Z
      GAS4=(-4.84164D01+2.23133D01*Z-9.19118D-01*Y)*Y*Y
      GAS5=-1.80898D05+2.82532D05*Z
      GAS6=(-1.01053D03+9.75576D02*Z)*Y
      GAS7=(-1.47220D05-2.33631D02*Y+2.55940D04*Z)*Z*Z
      GAS8=(3.28681D00-1.76588D00*Z-1.54962D-01*Y)*Y*Y
      XXX=-4.104D01+6.507D01*Y+2.083D01*Z-3.472D01*Z*Y
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  100 IF (Z.GT.2.58D00) GO TO 110
      GAS1=5.1131824D04-6.664875D04*Z
      GAS2=(2.02171D03-1.9306292D03*Z)*Y
      GAS3=(2.8762395D04+4.3353467D02*Y-4.1064609D03*Z)*Z*Z
      GAS4=(-8.4970047D01+1.7925919D01*Z-6.2576542D00*Y)*Y*Y
      GAS5=-6.2768156D04+8.6015875D04*Z
      GAS6=(-1.0002036D03+6.2537280D02*Z)*Y
      GAS7=(-3.957827D04-3.8467377D01*Y+6.12953D03*Z)*Z*Z
      GAS8=(-1.0591702D02+7.636142D01*Z+5.938859D00*Y)*Y*Y
      XXX=-3.901D00+2.418D01*Y+1.374D00*Z-1.145D01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  110 IF (Z.GT.2.73D00) GO TO 120
      GAS1=1.0088046D06-1.086321D06*Z
      GAS2=(1.3844801D04-9.7268516D03*Z)*Y
      GAS3=(3.8985325D05+1.7091665D03*Y-4.6621066D04*Z)*Z*Z
      GAS4=(1.4840726D02-5.2645004D01*Z-1.5477133D-01*Y)*Y*Y
      GAS5=-1.073351D06+1.14571D06*Z
      GAS6=(-1.9343957D04+1.3366211D04*Z)*Y
      GAS7=(-4.0670987D05-2.2955198D03*Y+4.7999871D04*Z)*Z*Z
      GAS8=(-4.1016724D02+1.4994148D02*Z-1.9779787D00*Y)*Y*Y
      XXX=-1.026D02+6.302D01*Y+3.819D01*Z-2.431D01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  120 GAS1=-9.6638500D04+1.3206488D04*Z
      GAS2=(-4.7458105D04+2.3596875D04*Z)*Y
      GAS3=(1.8602773D04-2.306802D03*Y-4.0413552D03*Z)*Z*Z
      GAS4=(-5.3564258D03+2.2433904D03*Z+2.5188145D02*Y)*Y*Y
      GAS5=1.0962581D05-2.990116D04*Z
      GAS6=(4.7883496D04-2.3785383D04*Z)*Y
      GAS7=(-1.1753969D04+2.2905522D03*Y+3.1304399D03*Z)*Z*Z
      GAS8=(5.473418D03-2.3208018D03*Z-2.6570068D02*Y)*Y*Y
      XXX=-3.107D01+1.082D01*Y+1.047D01*Z-3.047D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0-GAS9)
      GO TO 200
C
  130 IF (Z.GT.1.40D00) GO TO 140
      GAS1=-1.58386D03+3.49223D03*Z
      GAS2=(-8.39834D02+1.09565D03*Z)*Y
      GAS3=(-2.56175D03-3.56197D02*Y+6.25145D02*Z)*Z*Z
      GAS4=(-1.22407D01+7.65634D00*Z+2.58235D-01*Y)*Y*Y
      GAS5=1.58025D03-3.47664D03*Z
      GAS6=(8.39588D02-1.09490D03*Z)*Y
      GAS7=(2.54682D03+3.55674D02*Y-6.18504D02*Z)*Z*Z
      GAS8=(1.20843D01-7.44857D00*Z-2.91202D-01*Y)*Y*Y
      XXX=-2.171D01-4.342D00*Y+1.316D01*Z+2.632D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/
     &  (1.00D00-GAS9)
      GO TO 200
C
  140 IF (Z.GT.1.91D00) GO TO 150
      GAS1=7.89255D02-1.91743D03*Z
      GAS2=(3.59227D02-4.44070D02*Z)*Y
      GAS3=(1.39463D03+1.34083D02*Y-3.13446D02*Z)*Z*Z
      GAS4=(1.90681D01-1.09285D01*Z+4.24933D-02*Y)*Y*Y
      GAS5=-1.31401D03+3.13134D03*Z
      GAS6=(-5.18755D02+6.80268D02*Z)*Y
      GAS7=(-2.32493D03-2.21393D02*Y+5.52563D02*Z)*Z*Z
      GAS8=(-3.32001D01+2.11819D01*Z-4.75163D-01*Y)*Y*Y
      XXX=-5.025D01-8.412D00*Y+2.982D01*Z+3.509D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  150 IF (Z.GT.2.05D00) GO TO 160
      GAS1=3.58691D04-5.16852D04*Z
      GAS2=(-6.30189D02+6.63314D02*Z)*Y
      GAS3=(2.47471D04-1.73538D02*Y-3.93167D03*Z)*Z*Z
      GAS4=(-4.23871D01+2.08048D01*Z-1.05512D00*Y)*Y*Y
      GAS5=-1.10522D05+1.67591D05*Z
      GAS6=(4.61877D03-4.94930D03*Z)*Y
      GAS7=(-8.46558D04+1.32441D03*Y+1.42438D04*Z)*Z*Z
      GAS8=(2.25065D01-1.10316D01*Z+9.62887D-01*Y)*Y*Y
      XXX=-1.681D02+7.063D01*Y+8.75D01*Z-3.75D01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  160 IF (Z.GT.2.57D00) GO TO 170
      GAS1=3.1899562D04-4.2186664D04*Z
      GAS2=(2.3055603D03-1.9897017D03*Z)*Y
      GAS3=(1.849998D04+4.2561816D02*Y-2.6808696D03*Z)*Z*Z
      GAS4=(-1.6195114D01+5.8640623D00*Z-3.6172504D00*Y)*Y*Y
      GAS5=-5.7594039D04+7.9328437D04*Z
      GAS6=(-1.9275989D03+1.6730544D03*Z)*Y
      GAS7=(-3.6473008D04-3.6100732D02*Y+5.597543D03*Z)*Z*Z
      GAS8=(-7.920808D01+4.0542084D01*Z+2.1495867D00*Y)*Y*Y
      XXX=-5.733D01+2.088D01*Y+2.592D01*Z-9.793D00*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  170 IF (Z.GT.2.75D00) GO TO 180
      GAS1=7.0838087D05-7.5619919D05*Z
      GAS2=(3.9503091D03-2.7381802D03*Z)*Y
      GAS3=(2.6888181D05+4.7728687D02*Y-3.183816D04*Z)*Z*Z
      GAS4=(-1.2532251D02+4.7734787D01*Z-4.0148029D00*Y)*Y*Y
      GAS5=-2.5216325D05+2.1727769D05*Z
      GAS6=(9.2882383D03-7.780918D03*Z)*Y
      GAS7=(-5.6539297D04+1.6120212D03*Y+3.9419248D03*Z)*Z*Z
      GAS8=(1.8537296D02-7.1010757D01*Z+1.1307096D00*Y)*Y*Y
      XXX=-1.786D02+2.18D-01*Y+6.714D01*Z-4.739D-01*Y*Z
      IF (ABS(XXX).LT.ARGMAX) THEN
         GAS9=EXP(XXX)
      ELSE IF (XXX.GT.0.00D00) THEN
         GAS9=EXMAX
      ELSE
         GAS9=EXMIN
      END IF
      GO TO 190
C
  180 GAS1=3.1855037D05-3.3041156D05*Z
      GAS2=(2.2983352D04-1.6623461D04*Z)*Y
      GAS3=(1.13848D05+3.0098223D03*Y-1.3020133D04*Z)*Z*Z
      GAS4=(-1.8599039D02+6.9840683D01*Z-7.7371645D00*Y)*Y*Y
      F=GAS1+GAS2+GAS3+GAS4
      GO TO 200
C
  190 F=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/
     &  (1.00D00+GAS9)
  200 K=1.87915D-02*F
C
C     End subroutine UGAS4.
C
      RETURN
      END
