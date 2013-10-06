	SUBROUTINE EM(N, X1, X2, P, XMEAN, XSIGMA, E1, E2, MAXITS,
     *  COV, NOBS, K, IFAULT)
C
C	  ALGORITHM AS 138 APPL. STATIST. (1979) VOL.28, NO.2
C
C	  COMPUTE THE MAXIMUM LIKELIHOOD ESTIMATES OF THE MEAN
C	  AND STANDARD DEVIATION OF A SINGLE NORMAL POPULATION,
C	  THE DATA MAY CONTAIN CENSORED OR CONFINED OBSERVATIONS.
C
	implicit double precision (a-h,o-z)
	DIMENSION X1(N), X2(N), COV(2, 2)
	INTEGER P(N), NOBS(4)
C
	DATA C /0.39894228D0/, ONEPLS /1.0001D0/, TOL /0.00001D0/,
     *  TOLL /0.00001D0/
C
	IFAULT = -2
	IF (N .LT. 2) RETURN
C
C	  INITIALIZE COUNTERS
C
	K = 0
	SUM = 0.0D0
	SUM2 = 0.0D0
	SUMG = 0.0D0
	SUMG2 = 0.0D0
	IP = 0
	IQ = 0
	IR = 0
	IS = 0
C
C	  EXACTLY SPECIFIED OBSERVATIONS ARE REMOVED,
C	  THE REMAINING DATA PACKED INTO FIRST PART OF ARRAY X
C
	IFAULT = -4
	DO 200 I = 1, N
	IPT = P(I)
	IF (IPT .EQ. 0) GOTO 100
	IF (IPT .EQ. 2 .AND. DABS(X1(I) - X2(I)) .LE. DABS(X1(I) * TOLL))
     *  GOTO 100
C
C	  OBSERVATION NOT EXACTLY SPECIFIED
C
	IS = IS + 1
	P(IS) = IPT
	X1(IS) = X1(I)
C
C	  HANDLE GROUPED DATA
C
	IF (IPT .NE. 2) GOTO 50
	IQ = IQ + 1
	IF (X1(I) .GT. X2(I)) RETURN
	X2(IS) = X2(I)
	XTEMP = 0.5 * (X1(I) + X2(I))
	SUMG = SUMG + XTEMP
	SUMG2 = SUMG2 + XTEMP ** 2
	GOTO 200
C	  ACCUMULATE NUMBER OF OBSERVATIONS CENSORED ON THE RIGHT
C
50	IF (IPT .EQ. 1) IR = IR + 1
	GOTO 200
C
C	  HANDLE EXACTLY-SPECIFIED OBSERVATIONS
C
100	IP = IP + 1
	XTEMP = X1(I)
	SUM = SUM + XTEMP
	SUM2 = SUM2 + XTEMP ** 2
200	CONTINUE
C
C	  INITIAL PASS THROUGH DATA COMPLETED
C
	NOBS(1) = IP
	NOBS(2) = IR
	NOBS(3) = N - IP - IR - IQ
	NOBS(4) = IQ
	RIM = IP + IQ
	IF (IP .EQ. N) GOTO 230
	IF (XSIGMA .GT. 0.0) GOTO 350
	IF (RIM .GT. ONEPLS) GOTO 250
C
C	  AT MOST ONE OBSERVATION HAS BEEN EXACTLY
C	  SPECIFIED OR CONFINED
C
	XMEAN = 1.0D0
	XSIGMA = 1.0D0
	GOTO 350
C
C	  ALL OBSERVATIONS EXACTLY SPECIFIED
C
230	XMEAN = SUM / RIM
	XSIGMA = DSQRT((SUM2 - RIM * XMEAN ** 2) / RIM)
	COV(1, 1) = XSIGMA ** 2 / RIM
	COV(2, 2) = COV(1, 1) * 0.5D0
	COV(1, 2) = 0.0D0
	COV(2, 1) = 0.0D0
C
C	  NORMAL RETURN
C
240	IFAULT = 0
	RETURN
C
C	  OBTAIN INITIAL ESTIMATES
C
250	XMEAN = (SUM + SUMG) / RIM
	XSIGMA = DSQRT((SUM2 + SUMG2 - RIM * XMEAN ** 2) / RIM)
C
C	  INITIALIZE BEFORE STARTING FIRST ITERATION
C
350	RP = IP
	RN = N
C
C	  START OF ITERATION CYCLE,
C	  ESTIMATE CONDITIONAL EXPECTATION OF CONFINED AND
C	  CENSORED OBSERVATIONS
C
	IFAULT = -3
400	TS = SUM
	SUMG2 = SUM2
	TD = RP
	DO 610 I = 1, IS
	YS = (X1(I) - XMEAN) / XSIGMA
	IF (P(I) - 1) 500, 450, 550
C
C	  OBSERVATION CENSORED ON THE RIGHT
C
450	CALL RMILLS(YS, F, TOL)
	W = XMEAN + XSIGMA * F
	TD = TD + F * (F - YS)
	GOTO 600
C
C	  OBSERVATION CENSORED ON THE LEFT
C
500	CALL RMILLS(-YS, F, TOL)
	W = XMEAN - XSIGMA * F
	TD = TD + F * (F + YS)
	GOTO 600
C
C	  CONFINED OBSERVATION.
C	  USE MILLS RATIO RECIPROCAL TO COMPUTE PROBABILITY
C	  INTEGRALS THAT ARE REQUIRED,
C	  AS IN ORIGINAL ALGORITHM ASSUMING X1(I) IS
C	  NEVER GREATER THAN X2(I) FOR CONFINED OBSERVATIONS
C
550	YN = DEXP(-0.5 * YS ** 2) * C
	CALL RMILLS(YS, F, TOL)
	YQ = YN / F
	YSU = (X2(I) - XMEAN) / XSIGMA
	YNU = DEXP(-0.5 * YSU ** 2) * C
	CALL RMILLS(YSU, FU, TOL)
	YQU = YNU / FU
	YD = YQ - YQU
C
C	  IF INTEGRAL NOT EQUAL TO ZERO, CARRY ON
C
	IF (YD .LT. TOLL) RETURN
	A = (YN - YNU) / YD
	W = XMEAN + XSIGMA * A
	TD = TD + (A ** 2 + (YSU * YNU - YS * YN) / YD)
C
600	TS = TS + W
	SUMG2 = SUMG2 + W ** 2
610	CONTINUE
C
C	  CALCULATE NEW ESTIMATES
C
	XNEW = TS / RN
	YNEW = DSQRT((SUMG2 + RN * XMEAN ** 2 - 2.0 * TS * XMEAN) / TD)
	K = K + 1
	IF (DABS(XNEW - XMEAN) .LT. E1 .AND. DABS(YNEW - XSIGMA) .LT. E2)
     *  GOTO 700
	IF (K .GE. MAXITS) GOTO 650
	XMEAN = XNEW
	XSIGMA = YNEW
	GOTO 400
C
C	  MAXIMUM NUMBER OF ITERATIONS EXCEEDED
C
650	IFAULT = -1
	COV(1, 1) = 0.0D0
	COV(2, 2) = 0.0D0
	COV(1, 2) = XNEW - XMEAN
	COV(2, 1) = YNEW - XSIGMA
	RETURN
C
C	  CONVERGENCE OBTAINED
C
700	XMEAN = XNEW
	XSIGMA = YNEW
	XSIG2 = XSIGMA ** 2
C
C	  CALCULATE VARIANCE-COVARIANCE MATRIX
C
	X11 = RP
	X12 = (SUM - RP * XMEAN) / XSIGMA
	X22 = RP + (SUM2 + RP * XMEAN ** 2 - 2.0 * SUM * XMEAN) / XSIG2
	DO 800 I = 1, IS
	YS = (X1(I) - XMEAN) / XSIGMA
	IF (P(I) - 1) 740, 710, 770
710	CALL RMILLS(YS, F, TOL)
C
C	  OBSERVATION CENSORED ON THE RIGHT
C
	FL = F * (F - YS)
730	X11 = X11 + FL
	FL = FL * YS
	X12 = X12 + FL
	FL = FL * YS
	X22 = X22 + FL
	GOTO 800
C
740	CALL RMILLS(-YS, F, TOL)
C
C	  OBSERVATION CENSORED ON THE LEFT
C
	FL = F * (F + YS)
	GOTO 730
C
770	CALL RMILLS(YS, F, TOL)
C
C	  OBSERVATION CONFINED BETWEEN 2 FINITE LIMITS
C
	YN = DEXP(-0.5D0 * YS ** 2) * C
	YQ = YN / F
	YSU = (X2(I) - XMEAN) / XSIGMA
	CALL RMILLS(YSU, FU, TOL)
	YNU = DEXP(-0.5D0 * YSU ** 2) * C
	YQU = YNU / FU
	YD = YQ - YQU
	A = (YN - YNU) / YD
	B = (YNU * YSU - YN * YS) / YD
	X11 = X11 + A ** 2 + B
	B1 = (YS ** 2 * YN - YSU ** 2 * YNU) / YD
	X12 = X12 - A * B - B1
	B1 = (YS ** 3 * YN - YSU ** 3 * YNU) / YD
	X22 = X22 - B1 + B ** 2
800	CONTINUE
	CONST = XSIG2 / (X11 * X22 - X12 * X12)
	COV(1, 1) = CONST * X22
	COV(2, 2) = CONST * X11
	COV(1, 2) = -CONST * X12
	COV(2, 1) = COV(1, 2)
	GOTO 240
	END
C
	SUBROUTINE RMILLS(X, FUNC, TOL)
C
C	  ALGORITHM AS 138.1 APPL. STATIST. (1979) VOL.28 NO.2
C
C	  COMPUTE THE RECIPROCAL OF MILLS RATIO
C
	implicit double precision (a-h,o-z)
	DATA FPI /1.2533141D0/, FPII /0.7978846D0/
C
	FUNC = 0.0D0
	IF (X .LT. -10.0D0) RETURN
	FUNC = FPII
	Y = DABS(X)
	IF (Y .LT. 0.000001D0) RETURN
	SGN = 1.0
	IF (X .LT. 0.0D0) SGN = -1.0D0
	IF (Y .GT. 2.0D0) GOTO 100
	S = 0.0D0
	A = 1.0D0
	T = Y
	R = Y
	B = Y ** 2
40	A = A + 2.0
	S = T
	R = R * B / A
	T = T + R
	IF (R .GT. TOL) GOTO 40
	FUNC = 1.0D0 / (FPI * DEXP(0.5D0 * B) - SGN * T)
	RETURN
100	A = 2.0
	B1 = Y
	S = Y
	A1 = Y ** 2 + 1.0
	A2 = Y * (A1 + 2.0)
	B2 = A1 + 1.0
	T = A2 / B2
140	A = A + 1.0
	A0 = A1
	A1 = A2
	A2 = Y * A1 + A * A0
	B0 = B1
	B1 = B2
	B2 = Y * B1 + A * B0
	R = S
	S = T
	T = A2 / B2
	IF (T - R .GT. TOL .OR. T - S .GT. TOL) GOTO 140
	FUNC = T
	IF (SGN .LT. 0.0D0) FUNC =
     *  T / (2.0D0 * FPI * DEXP(0.5D0 * Y ** 2) * T - 1.0D0)
	RETURN
	END

      SUBROUTINE REGRES(N, Y1, Y2, P, MPLONE, X, ROWX, COLX, W, LENW,
     *    VCOV, WORK, LENWRK, ALPHA, TOL, MAXITS, IFAULT)
C***** NOTE: this routine uses the auxiliary routine
C***** AS7  with the modified argument lists suggested by Freeman 
C***** Vol 31 No 3 p336--339.
C
C  as r31 Vol 29 No 2 1980 p228 -- incorporated
C  as r32 Vol 30 No 1 1981 p105 -- incorporated
C  as r91 Vol 42 No 3 1993 p583 -- incorporated
C
C        ALGORITHM AS139 APPL. STATIST. (1979) VOL.28, NO.2
C
C        COMPUTE MAXIMUM LIKELIHOOD ESTIMATES
C        FROM A LINEAR MODEL WITH NORMAL HETEROGENOUS VARIANCE.
C        THE DESIGN MATRIX MUST BE NON-SINGULAR.  THE DEPENDENT
C        VARIABLE MAY INCLUDE OBSERVATIONS CENSORED IN EITHER TAIL
C        AND/OR OBSERVATIONS CONFINED BETWEEN FINITE LIMITS.
C
C     Auxiliary routine required: RMILLS which is the last routine in
C     AS138.
C
      IMPLICIT NONE

      INTEGER N, MPLONE, ROWX, COLX, P(N), LENW, LENWRK, MAXITS, IFAULT
      DOUBLE PRECISION X(ROWX,COLX), TOL(MPLONE), Y1(N), Y2(N),
     *    ALPHA(MPLONE)
      DOUBLE PRECISION VCOV(LENWRK), WORK(LENWRK), W(LENW)
C
C     Local variables
C
      INTEGER M, I, J, K, II, NDIMC, NUL, JJ, IIK, IJ, IPT, NJ, IIJ,
     *      JJI, KK
      DOUBLE PRECISION QLIMIT, RLIMIT, C, ZERO, HALF, ONE
      PARAMETER (QLIMIT = 0.00001D0, RLIMIT=0.00001D0, C=0.39894228D0)
      PARAMETER (zero=0.d0, half=0.5D0, one=1.D0)
      DOUBLE PRECISION TEMP, SUM2, DEMP, R, R2, XSIG, TD, YMEAN, F, A,
     *      FU, YN, YQ, TMPU, YNU, YQU, TINT, TEMP2, YS, YSU, B
C
C        INITIALIZATION
C
      M = MPLONE-1
C
C        CHECK ARRAY SIZES, ETC
C
      IFAULT = -7
      IF(ROWX.LT.N) RETURN
      IFAULT = -8
      IF(COLX.LT.M) RETURN
      IFAULT = -9
      IF(LENW.LT.(MPLONE+N)) RETURN
      IFAULT = -10
      IF(LENWRK.LT.(M*N)) RETURN
C
C        COMPUTE X'X IN LOWER TRIANGULAR FORM
C
      II = 0
      DO 53 I = 1, M
        DO 50 J = 1, I
          TEMP = ZERO
          DO 40 K = 1, N
   40     TEMP = TEMP + X(K,I)*X(K,J)
          II = II + 1
          VCOV(II) = TEMP
   50   CONTINUE
   53 CONTINUE
      NDIMC = M*(M+1)/2
      CALL SYMINV(VCOV, M, NDIMC, WORK, W, NUL, IFAULT)
      IF(IFAULT.NE.0) GOTO 60
      IF(NUL.EQ.0) GOTO 70
   60 VCOV(2) = IFAULT
      VCOV(1) = NUL
      IFAULT = -5
      RETURN
C
C        MATRIX NON-SINGULAR AND INVERSE OBTAINED
C        COMPUTE (X'X)INVERSE*X
C        FOLLOWING SCHEME USED TO REDUCE NUMBER OF STORAGE ARRAYS
C        NEEDED.   EXPAND FROM TRIANGULAR TO SQUARE MATRIX
C
   70 CALL UNPACK(WORK, M, LENWRK)
C
C        DO MULTIPLICATION - ONE ROW AT A TIME - STARTING
C        WITH THE LAST ONE
C
      JJ = N*M
      II = M*M
      DO 220 I = 1, M
        II = II - M
        DO 200 J = 1, N
          TEMP = ZERO
          DO 170 K = 1, M
            IIK = II + K
            TEMP = TEMP + WORK(IIK)*X(J,K)
  170     CONTINUE
          W(J) = TEMP
  200   CONTINUE
        DO 210 J = 1, N
          IJ = N+1-J
          WORK(JJ) = W(IJ)
          JJ = JJ-1
  210   CONTINUE
  220 CONTINUE
C
      XSIG = ALPHA(MPLONE)
      IF(XSIG.GT.ZERO) GOTO 500
C
C        NO ACCEPTABLE INITIAL VALUE FOR SIGMA HAS BEEN INPUT,
C        OBTAIN INITIAL ESTIMATES FROM EXACTLY SPECIFIED
C        OBSERVATIONS ONLY (ALTHOUGH MATRIX BASED ON ALL
C        OBSERVATIONS) AND CONFINED OBSERVATIONS
C
      II = -N
      DO 300 I = 1, M
        II = II + N
        TEMP = ZERO
        DO 280 J = 1, N
          IIJ = II + J
          IPT = P(J)
          IF(IPT.EQ.0) GOTO 270
          IF(IPT.EQ.2) TEMP = TEMP + WORK(IIJ)*(Y1(J) + Y2(J))*HALF
          GOTO 280
  270     TEMP = TEMP + WORK(IIJ)*Y1(J)
  280   CONTINUE
        ALPHA(I) = TEMP
  300 CONTINUE
C
C        CALCULATE INITIAL ESTIMATE OF SIGMA
C
      SUM2 = ZERO
      TEMP = ZERO
      DO 350 I = 1, N
        IPT = P(I)
        IF(IABS(IPT).EQ.1) GOTO 350
        DEMP = Y1(I)
        IF(IPT.EQ.2) DEMP = (DEMP + Y2(I))*HALF
        DO 320 J = 1, M
  320   DEMP = DEMP - ALPHA(J)*X(I, J)
        SUM2 = SUM2 + DEMP**2
        TEMP = TEMP + ONE
  350 CONTINUE
      XSIG = SQRT(SUM2/TEMP)
C
C        COMPUTE SOME CONSTANTS NEEDED THROUGHOUT
C
  500 R = ZERO
      R2 = ZERO
      IFAULT = -2
      DO 600 I = 1, N
        IPT = P(I)
        IF(IPT.EQ.0) GOTO 550
        IF(IPT.EQ.2 .AND. ABS(Y1(I)-Y2(I)).LE. QLIMIT*ABS(Y1(I)))
     *                     GOTO 540
        IF(IPT.NE.2) GOTO 600
        R2 = R2 + ONE
        IF(Y1(I).LT.Y2(I)) GOTO 600
        RETURN
  540   P(I) = 0
  550   R = R + ONE
        W(I) = Y1(I)
  600 CONTINUE
      I = R + R2 + 0.01
      IFAULT = -4
      IF(I.LT.MPLONE) RETURN
      IFAULT = 0
C
C        START OF ITERATION PROCEDURE
C
  620 TD = R
      SUM2 = ZERO
C
C        COMPLETE W-VECTOR
C
      DO 1000 I = 1, N
        IPT = P(I)
        YMEAN = ZERO
        DO 650 J = 1, M
  650   YMEAN = YMEAN + ALPHA(J)*X(I,J)
        IF(IPT.EQ.0) GOTO 990
C
C        OBSERVATION NOT EXACTLY SPECIFIED
C
        TEMP = (Y1(I) - YMEAN)/XSIG
        IF(IPT-1)750, 700, 800
C
C        OBSERVATION CENSORED FROM ABOVE - LOWER BOUND KNOWN
C
  700   CALL RMILLS(TEMP, F, RLIMIT)
        W(I) = YMEAN + XSIG*F
        TD = TD + F*(F-TEMP)
        GOTO 990
C
C        OBSERVATION CENSORED FROM BELOW - UPPER BOUND KNOWN
C
  750   CALL RMILLS(-TEMP, F, RLIMIT)
        W(I) = YMEAN - XSIG*F
        TD = TD + F*(F+TEMP)
        GOTO 990
C
C        OBSERVATION CONFINED TO LIE BETWEEN TWO FINITE LIMITS
C
  800   YN = EXP(-HALF*TEMP**2)*C
        CALL RMILLS(TEMP, F, RLIMIT)
        YQ = YN/F
        TMPU = (Y2(I) - YMEAN)/XSIG
        YNU = EXP(-HALF*TMPU**2)*C
        CALL RMILLS(TMPU, FU, RLIMIT)
        YQU = YNU/FU
        TINT = YQ - YQU
        IF(TINT.GE.QLIMIT) GOTO 820
C
C        AFTER STANDARDIZING, UPPER AND LOWER LIMITS RESULT IN
C        SAME PROBABILITY INTEGRAL
C
        IFAULT = -3
        RETURN
  820   A = (YN - YNU)/TINT
        W(I) = YMEAN + XSIG*A
        TD = TD + (A**2 + (TMPU*YNU - TEMP*YN)/TINT)
C
C        CALCULATE RESIDUAL SUM OF SQUARES
C
  990   SUM2 = SUM2 + (W(I) - YMEAN)**2
 1000 CONTINUE
C
C        UPDATE PARAMETER ESTIMATES - STORE IN END OF W-VECTOR
C
      JJ = -N
      DO 1200 J = 1, M
        JJ = JJ + N
        TEMP = ZERO
        DO 1100 I = 1, N
          JJI = JJ + I
          TEMP = TEMP + WORK(JJI)*W(I)
 1100   CONTINUE
        NJ = N + J
        W(NJ) = TEMP
 1200 CONTINUE
      NJ = N + MPLONE
      W(NJ) = SQRT(SUM2/TD)
C
C        TEST FOR CONVERGENCE
C
      DO 1300 J = 1, MPLONE
        NJ = N + J
        IF(ABS(ALPHA(J)-W(NJ)).GE.TOL(J)) GOTO 1400
 1300 CONTINUE
C
C        IF WE REACH HERE, CONVERGENCE OBTAINED
C
      IJ = IFAULT
      IFAULT = -1
C
C        UPDATE VALUES
C
 1400 DO 1450 J = 1, MPLONE
        NJ = N + J
        ALPHA(J) = W(NJ)
 1450 CONTINUE
      XSIG = ALPHA(MPLONE)
      IFAULT = IFAULT + 1
      IF(IFAULT.EQ.0) GOTO 1600
      IF(IFAULT.LE.MAXITS) GOTO 620
      IFAULT = -1
      RETURN
C
C        CONVERGENCE OBTAINED - COMPUTE VARIANCE-COVARIANCE
C        MATRIX. INITIALIZE WORK ARRAY
C
 1600 II = MPLONE*(MPLONE + 1)/2
      DO 1650 I = 1, II
 1650 WORK(I) = ZERO
      DO 2500 I = 1, N
        IPT = P(I)
        YS = Y1(I)
        DO 1680 J = 1, M
 1680   YS = YS - ALPHA(J)*X(I,J)
        YS = YS/XSIG
        JJ = 0
        IF(IPT.NE.0) GOTO 1900
C
C        EXACTLY SPECIFIED OBSERVATION
C
        DO 1750 K = 1, M
          DO 1720 J = 1, K
            JJ = JJ + 1
            WORK(JJ) = WORK(JJ) + X(I,K)*X(I,J)
 1720     CONTINUE
          KK = II - MPLONE + K
          WORK(KK) = WORK(KK) + YS*X(I,K)
 1750   CONTINUE
        WORK(II) = WORK(II) + ONE + YS**2
        GOTO 2500
 1900   IF(IPT-1) 2100, 2000, 2300
C
C        OBSERVATION CONSORED FROM ABOVE - LOWER BOUND KNOWN
C
 2000   CALL RMILLS(YS, F, RLIMIT)
        TEMP = F*(F - YS)
        GOTO 2150
C
C        OBSERVATION CONSORED FROM BELOW - UPPER BOUND KNOWN
C
 2100   CALL RMILLS(-YS, F, RLIMIT)
        TEMP = F*(F + YS)
C
C        ROUTINE FOR CENSORED OBSERVATIONS
C
 2150   DO 2190 K = 1, M
          DO 2170 J = 1, K
            JJ = JJ + 1
            WORK(JJ) = WORK(JJ) + X(I,J)*X(I,K)*TEMP
 2170     CONTINUE
          KK = II - MPLONE + K
          WORK(KK) = WORK(KK) + YS*X(I,K)*TEMP
 2190   CONTINUE
        WORK(II) = WORK(II) + YS**2*TEMP
        GOTO 2500
C
C        OBSERVATION CONFINED BETWEEN TWO FINITE LIMITS
C
 2300   YN = EXP(-HALF*YS**2)*C
        CALL RMILLS(YS, F, RLIMIT)
        YQ = YN/F
        YSU = YS + (Y2(I) - Y1(I))/XSIG
        CALL RMILLS(YSU, FU, RLIMIT)
        YNU = EXP(-HALF*YSU**2)*C
        YQU = YNU/FU
        TINT = YQ - YQU
        A = (YN - YNU)/TINT
        B = (YNU*YSU - YN*YS)/TINT
        TEMP = A**2 + B
        TEMP2 = A*B + (YS**2 * YN - YSU**2 * YNU) / TINT
        DO 2350 K = 1, M
          DO 2330 J = 1, K
            JJ = JJ + 1
            WORK(JJ) = WORK(JJ) + X(I,J)*X(I,K)*TEMP
 2330     CONTINUE
          KK = II - MPLONE + K
          WORK(KK) = WORK(KK) - X(I,K) * TEMP2
 2350   CONTINUE
        TEMP = (YS**3*YN - YSU**3*YNU)/TINT
        WORK(II) = WORK(II) - TEMP + B**2
 2500 CONTINUE
C
C        INVERT THE MATRIX
C
      CALL SYMINV(WORK, MPLONE, II, VCOV, W, NUL, IFAULT)
      IF(IFAULT.EQ.0 .AND. NUL.EQ.0) GOTO 2550
      VCOV(2) = IFAULT
      VCOV(1) = NUL
      IFAULT = -6
      RETURN
C
C        RESTORE ITERATION COUNTER
C
 2550 IFAULT = IJ
C
C        MULTIPLY BY SIGMA-SQUARED
C
      TEMP = XSIG**2
      DO 2580 I = 1, II
 2580 VCOV(I) = VCOV(I)*TEMP
C
C        UNPACK THE MATRIX
C
      CALL UNPACK(VCOV, MPLONE, LENWRK)
      RETURN
      END
*****************************************************************
      SUBROUTINE UNPACK(X, N, LENX)
C
C        ALGORITHM AS139.1 APPL. STATIST. (1979) VOL.28 NO.2
C
C        THIS SUBROUTINE EXPANDS A SYMMETRIC MATRIX STORED IN LOWER
C        TRIANGULAR FORM IN THE FIRST N*(N+1)/2 POSITIONS OF X
C        INTO A MATRIX USING THE FIRST N*N POSITIONS
C
C        LENX = LENGTH OF VECTOR X - MUST BE LESS THAN N*N
C
      INTEGER N, LENX
      DOUBLE PRECISION X(LENX)
C
C     Local variables
C
      INTEGER NSQ, II, JJ, I, IJ, KK, J

      NSQ = N*N
      II = NSQ
      JJ = N*(N+1)/2
C
C        STORE LAST ROW
C
      DO 10 I = 1, N
        X(II) = X(JJ)
        II = II-1
        JJ = JJ-1
   10 CONTINUE
      DO 80 I = 2, N
C
C        OBTAIN UPPER PART OF MATRIX FROM PART ALREADY SHIFTED
C
        IJ = I - 1
        KK = NSQ+1-I
        DO 50 J = 1, IJ
          X(II) = X(KK)
          II = II - 1
          KK = KK - N
   50   CONTINUE
C
C        OBTAIN LOWER PART OF MATRIX FROM
C        ORIGINAL TRIANGULAR STORAGE
C
        IJ = N - IJ
        DO 70 J = 1, IJ
          X(II) = X(JJ)
          II = II - 1
          JJ = JJ - 1
   70   CONTINUE
   80 CONTINUE
      RETURN
      END

        subroutine syminv(a, n, nn, c, w, nullty, ifault)
c
c       Algorithm AS7, Applied Statistics, vol.17, 1968, p.198.
c
c       Forms in c( ) as lower triangle, a generalised inverse
c       of the positive semi-definite symmetric matrix a( )
c       order n, stored as lower triangle.
c
c       arguments:-
c       a()     = input, the symmetric matrix to be inverted, stored in
c                 lower triangular form
c       n       = input, order of the matrix
c       nn      = input, the size of the a and c arrays     n*(n+1)/2
c       c()     = output, the inverse of a (a generalized inverse if c is
c                 singular), also stored in lower triangular.
c                 c and a may occupy the same locations.
c       w()     = workspace, dimension at least n.
c       nullty  = output, the rank deficiency of a.
c       ifault  = output, error indicator
c                       = 1 if n < 1
c                       = 2 if a is not +ve semi-definite
c                       = 3 if nn < n*(n+1)/2
c                       = 0 otherwise
c
c***************************************************************************
c
        double precision a(nn), c(nn), w(n), x, zero, one
c
        data zero, one /0.0d0, 1.0d0/
c
c       cholesky factorization of a, result in c
c
        call chol(a, n, nn, c, nullty, ifault)
        if(ifault.ne.0) return
c
c       invert c & form the product (cinv)'*cinv, where cinv is the inverse
c       of c, row by row starting with the last row.
c       irow = the row number, ndiag = location of last element in the row.
c
        irow=n
        ndiag=nn
   10   l=ndiag
        if (c(ndiag) .eq. zero) goto 60
        do 20 i=irow,n
          w(i)=c(l)
          l=l+i
   20   continue
        icol=n
        jcol=nn
        mdiag=nn
   30   l=jcol
        x=zero
        if(icol.eq.irow) x=one/w(irow)
        k=n
   40   if(k.eq.irow) go to 50
        x=x-w(k)*c(l)
        k=k-1
        l=l-1
        if(l.gt.mdiag) l=l-k+1
        go to 40
   50   c(l)=x/w(irow)
        if(icol.eq.irow) go to 80
        mdiag=mdiag-icol
        icol=icol-1
        jcol=jcol-1
        go to 30
   60   do 70 j=irow,n
          c(l)=zero
          l=l+j
   70   continue
   80   ndiag=ndiag-irow
        irow=irow-1
        if(irow.ne.0) go to 10
        return
        end



      SUBROUTINE CHOL (A, N, NN, U, NULLTY, IFAULT)
C 
C       Algorithm AS6, Applied Statistics, vol.17, (1968)
C 
C       Given a symmetric matrix order n as lower triangle in a( )
C       calculates an upper triangle, u( ), such that uprime * u = a.
C       a must be positive semi-definite.  eta is set to multiplying
C       factor determining effective zero for pivot.
C 
C       arguments:-
C       a()     = input, a +ve definite matrix stored in lower-triangula
C                 form.
C       n       = input, the order of a
C       nn      = input, the size of the a and u arrays      n*(n+1)/2
C       u()     = output, a lower triangular matrix such that u*u' = a.
C                 a & u may occupy the same locations.
C       nullty  = output, the rank deficiency of a.
C       ifault  = output, error indicator
C                       = 1 if n < 1
C                       = 2 if a is not +ve semi-definite
C                       = 3 if nn < n*(n+1)/2
C                       = 0 otherwise
C 
C***********************************************************************
C 
      DOUBLE PRECISION A(NN), U(NN), ETA, ETA2, X, W, ZERO
C 
C       The value of eta will depend on the word-length of the
C       computer being used.  See introductory text.
C 
      DATA ETA, ZERO/1.D-9, 0.0D0/
C 
      IFAULT = 1
      IF (N.LE.0) RETURN
      IFAULT = 3
      IF (NN.LT.N*(N+1)/2) RETURN
      IFAULT = 2
      NULLTY = 0
      J = 1
      K = 0
      ETA2 = ETA*ETA
      II = 0
C 
C       Factorize column by column, icol = column no.
C 
      DO 80 ICOL = 1,N
        II = II+ICOL
        X = ETA2*A(II)
        L = 0
        KK = 0
C 
C       IROW = row number within column ICOL
C 
        DO 40 IROW = 1,ICOL
          KK = KK+IROW
          K = K+1
          W = A(K)
          M = J
          DO 10 I = 1,IROW
            L = L+1
            IF (I.EQ.IROW) GO TO 20
            W = W-U(L)*U(M)
            M = M+1
 10       CONTINUE
 20       IF (IROW.EQ.ICOL) GO TO 50
          IF (U(L).EQ.ZERO) GO TO 30
          U(K) = W/U(L)
          GO TO 40
 30       IF (W*W.GT.ABS(X*A(KK))) RETURN
          U(K) = ZERO
 40     CONTINUE
 50     IF (ABS(W).LE.ABS(ETA*A(K))) GO TO 60
        IF (W.LT.ZERO) RETURN
        U(K) = SQRT(W)
        GO TO 70
 60     U(K) = ZERO
        NULLTY = NULLTY+1
 70     J = J+ICOL
 80   CONTINUE
      IFAULT = 0
      END
