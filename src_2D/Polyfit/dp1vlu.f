*DECK DP1VLU
      SUBROUTINE DP1VLU (L, NDER, X, YFIT, YP, A)
C***BEGIN PROLOGUE  DP1VLU
C***PURPOSE  Use the coefficients generated by DPOLFT to evaluate the
C            polynomial fit of degree L, along with the first NDER of
C            its derivatives, at a specified point.
C***LIBRARY   SLATEC
C***CATEGORY  K6
C***TYPE      DOUBLE PRECISION (PVALUE-S, DP1VLU-D)
C***KEYWORDS  CURVE FITTING, LEAST SQUARES, POLYNOMIAL APPROXIMATION
C***AUTHOR  Shampine, L. F., (SNLA)
C           Davenport, S. M., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     The subroutine  DP1VLU  uses the coefficients generated by  DPOLFT
C     to evaluate the polynomial fit of degree  L , along with the first
C     NDER  of its derivatives, at a specified point.  Computationally
C     stable recurrence relations are used to perform this task.
C
C     The parameters for  DP1VLU  are
C
C     Input -- ALL TYPE REAL variables are DOUBLE PRECISION
C         L -      the degree of polynomial to be evaluated.  L  may be
C                  any non-negative integer which is less than or equal
C                  to  NDEG , the highest degree polynomial provided
C                  by  DPOLFT .
C         NDER -   the number of derivatives to be evaluated.  NDER
C                  may be 0 or any positive value.  If NDER is less
C                  than 0, it will be treated as 0.
C         X -      the argument at which the polynomial and its
C                  derivatives are to be evaluated.
C         A -      work and output array containing values from last
C                  call to  DPOLFT .
C
C     Output -- ALL TYPE REAL variables are DOUBLE PRECISION
C         YFIT -   value of the fitting polynomial of degree  L  at  X
C         YP -     array containing the first through  NDER  derivatives
C                  of the polynomial of degree  L .  YP  must be
C                  dimensioned at least  NDER  in the calling program.
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DP1VLU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I,IC,ILO,IN,INP1,IUP,K1,K1I,K2,K3,K3P1,K3PN,K4,K4P1,K4PN,
     * KC,L,LM1,LP1,MAXORD,N,NDER,NDO,NDP1,NORD
      DOUBLE PRECISION A(*),CC,DIF,VAL,X,YFIT,YP(*)
      CHARACTER*8 XERN1, XERN2
C***FIRST EXECUTABLE STATEMENT  DP1VLU
c      IF (L .LT. 0) GO TO 12
      NDO = MAX(NDER,0)
      NDO = MIN(NDO,L)
      MAXORD = A(1) + 0.5D0
      K1 = MAXORD + 1
      K2 = K1 + MAXORD
      K3 = K2 + MAXORD + 2
      NORD = A(K3) + 0.5D0
      IF (L .GT. NORD) GO TO 11
      K4 = K3 + L + 1
      IF (NDER .LT. 1) GO TO 2
      DO 1 I = 1,NDER
 1      YP(I) = 0.0D0
 2    IF (L .GE. 2) GO TO 4
      IF (L .EQ. 1) GO TO 3
C
C L IS 0
C
      VAL = A(K2+1)
      GO TO 10
C
C L IS 1
C
 3    CC = A(K2+2)
      VAL = A(K2+1) + (X-A(2))*CC
      IF (NDER .GE. 1) YP(1) = CC
      GO TO 10
C
C L IS GREATER THAN 1
C
 4    NDP1 = NDO + 1
      K3P1 = K3 + 1
      K4P1 = K4 + 1
      LP1 = L + 1
      LM1 = L - 1
      ILO = K3 + 3
      IUP = K4 + NDP1
      DO 5 I = ILO,IUP
 5      A(I) = 0.0D0
      DIF = X - A(LP1)
      KC = K2 + LP1
      A(K4P1) = A(KC)
      A(K3P1) = A(KC-1) + DIF*A(K4P1)
      A(K3+2) = A(K4P1)
C
C EVALUATE RECURRENCE RELATIONS FOR FUNCTION VALUE AND DERIVATIVES
C
      DO 9 I = 1,LM1
        IN = L - I
        INP1 = IN + 1
        K1I = K1 + INP1
        IC = K2 + IN
        DIF = X - A(INP1)
        VAL = A(IC) + DIF*A(K3P1) - A(K1I)*A(K4P1)
        IF (NDO .LE. 0) GO TO 8
        DO 6 N = 1,NDO
          K3PN = K3P1 + N
          K4PN = K4P1 + N
 6        YP(N) = DIF*A(K3PN) + N*A(K3PN-1) - A(K1I)*A(K4PN)
C
C SAVE VALUES NEEDED FOR NEXT EVALUATION OF RECURRENCE RELATIONS
C
        DO 7 N = 1,NDO
          K3PN = K3P1 + N
          K4PN = K4P1 + N
          A(K4PN) = A(K3PN)
 7        A(K3PN) = YP(N)
 8      A(K4P1) = A(K3P1)
 9      A(K3P1) = VAL
C
C NORMAL RETURN OR ABORT DUE TO ERROR
C
 10   YFIT = VAL
      RETURN
C
   11 WRITE (XERN1, '(I8)') L
      WRITE (XERN2, '(I8)') NORD
      WRITE(*,*) 'THE ORDER OF POLYNOMIAL EVALUATION, L = ' // XERN1 //
     * ' REQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD = ' // XERN2 //
     * ', COMPUTED BY DPOLFT -- EXECUTION TERMINATED.'
c      CALL XERMSG ('SLATEC', 'DP1VLU',
c     *   'THE ORDER OF POLYNOMIAL EVALUATION, L = ' // XERN1 //
c     *   ' REQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD = ' // XERN2 //
c     *   ', COMPUTED BY DPOLFT -- EXECUTION TERMINATED.', 8, 2)
      RETURN
C
c   12 CALL XERMSG ('SLATEC', 'DP1VLU',
c     +   'INVALID INPUT PARAMETER.  ORDER OF POLYNOMIAL EVALUATION ' //
c     +   'REQUESTED IS NEGATIVE.', 2, 2)
c      RETURN
      END
