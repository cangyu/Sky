C  This file is part of Construct2D.

C  Construct2D is free software: you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation, either version 3 of the License, or
C  (at your option) any later version.

C  Construct2D is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with Construct2D.  If not, see <http://www.gnu.org/licenses/>.

C Contains xfoil paneling subroutine (pangen) and dependencies
C Only minor changes made to xfoil source code

C  Copyright (C) 2013 -- 2016 Daniel Prosser (this modified version of 
C  XFoil code)
C  Original copyright (C) 2000 Mark Drela (original XFoil code)

      SUBROUTINE SPLIND(X,XS,S,N,XS1,XS2)
      DIMENSION X(N),XS(N),S(N)
      PARAMETER (NMAX=1000)
      DIMENSION  A(NMAX),B(NMAX),C(NMAX)
C-------------------------------------------------------
C     Calculates spline coefficients for X(S).          |
C     Specified 1st derivative and/or usual zero 2nd    |
C     derivative end conditions are used.               |
C     To evaluate the spline at some value of S,        |
C     use SEVAL and/or DEVAL.                           |
C                                                       |
C     S        independent variable array (input)       |
C     X        dependent variable array   (input)       |
C     XS       dX/dS array                (calculated)  |
C     N        number of points           (input)       |
C     XS1,XS2  endpoint derivatives       (input)       |
C              If = 999.0, then usual zero second       |
C              derivative end condition(s) are used     |
C              If = -999.0, then zero third             |
C              derivative end condition(s) are used     |
C                                                       |
C-------------------------------------------------------
      IF(N.GT.NMAX) STOP 'SPLIND: array overflow, increase NMAX'
C     
      DO 1 I=2, N-1
        DSM = S(I) - S(I-1)
        DSP = S(I+1) - S(I)
        B(I) = DSP
        A(I) = 2.0*(DSM+DSP)
        C(I) = DSM
        XS(I) = 3.0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
    1 CONTINUE
C
      IF(XS1.EQ.999.0) THEN
C----- set zero second derivative end condition
       A(1) = 2.0
       C(1) = 1.0
       XS(1) = 3.0*(X(2)-X(1)) / (S(2)-S(1))
      ELSE IF(XS1.EQ.-999.0) THEN
C----- set zero third derivative end condition
       A(1) = 1.0
       C(1) = 1.0
       XS(1) = 2.0*(X(2)-X(1)) / (S(2)-S(1))
      ELSE
C----- set specified first derivative end condition
       A(1) = 1.0
       C(1) = 0.
       XS(1) = XS1
      ENDIF
C
      IF(XS2.EQ.999.0) THEN
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ELSE IF(XS2.EQ.-999.0) THEN
       B(N) = 1.0
       A(N) = 1.0
       XS(N) = 2.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ELSE
       A(N) = 1.0
       B(N) = 0.
       XS(N) = XS2
      ENDIF
C
      IF(N.EQ.2 .AND. XS1.EQ.-999.0 .AND. XS2.EQ.-999.0) THEN
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ENDIF
C
C---- solve for derivative array XS
      CALL TRISOL(A,B,C,XS,N)
C
      RETURN
      END ! SPLIND

      FUNCTION D2VAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates d2X/dS2(SS)                       |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      D2VAL = (6.*T-4.)*CX1 + (6.*T-2.0)*CX2
      D2VAL = D2VAL/DS**2
      RETURN
      END ! D2VAL

      SUBROUTINE SCALC(X,Y,S,N)
      DIMENSION X(N), Y(N), S(N)
C----------------------------------------
C     Calculates the arc length array S  |
C     for a 2-D array of points (X,Y).   |
C----------------------------------------
C
      S(1) = 0.
      DO 10 I=2, N
        S(I) = S(I-1) + SQRT((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
   10 CONTINUE
C
      RETURN
      END ! SCALC

      SUBROUTINE SEGSPL(X,XS,S,N)
C-----------------------------------------------
C     Splines X(S) array just like SPLINE,      |
C     but allows derivative discontinuities     |
C     at segment joints.  Segment joints are    |
C     defined by identical successive S values. |
C-----------------------------------------------
      DIMENSION X(N), XS(N), S(N)
C
      IF(S(1).EQ.S(2)  ) STOP 'SEGSPL:  First input point duplicated'
      IF(S(N).EQ.S(N-1)) STOP 'SEGSPL:  Last  input point duplicated'
C
      ISEG0 = 1
      DO 10 ISEG=2, N-2
        IF(S(ISEG).EQ.S(ISEG+1)) THEN
         NSEG = ISEG - ISEG0 + 1
         CALL SPLIND(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG,-999.0,-999.0)
         ISEG0 = ISEG+1
        ENDIF
   10 CONTINUE
C
      NSEG = N - ISEG0 + 1
      CALL SPLIND(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG,-999.0,-999.0)
C
      RETURN
      END ! SEGSPL

      FUNCTION CURV(SS,X,XS,Y,YS,S,N)
      DIMENSION X(N), XS(N), Y(N), YS(N), S(N)
C-----------------------------------------------
C     Calculates curvature of splined 2-D curve |
C     at S = SS                                 |
C                                               |
C     S        arc length array of curve        |
C     X, Y     coordinate arrays of curve       |
C     XS,YS    derivative arrays                |
C              (calculated earlier by SPLINE)   |
C-----------------------------------------------
C     
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
C
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      XD = X(I) - X(I-1) + (1.0-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.0)*CX2
      XDD = (6.0*T-4.0)*CX1 + (6.0*T-2.0)*CX2
C
      CY1 = DS*YS(I-1) - Y(I) + Y(I-1)
      CY2 = DS*YS(I)   - Y(I) + Y(I-1)
      YD = Y(I) - Y(I-1) + (1.0-4.0*T+3.0*T*T)*CY1 + T*(3.0*T-2.0)*CY2
      YDD = (6.0*T-4.0)*CY1 + (6.0*T-2.0)*CY2
C 
      SD = SQRT(XD*XD + YD*YD)
      SD = MAX(SD,0.001*DS)
C
      CURV = (XD*YDD - YD*XDD) / SD**3
C
      RETURN
      END ! CURV

      SUBROUTINE LEFIND(SLE,X,XP,Y,YP,S,N)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C------------------------------------------------------
C     Locates leading edge spline-parameter value SLE
C
C     The defining condition is
C         
C      (X-XTE,Y-YTE) . (X',Y') = 0     at  S = SLE
C
C     i.e. the surface tangent is normal to the chord
C     line connecting X(SLE),Y(SLE) and the TE point.
C------------------------------------------------------
C
C---- convergence tolerance
      DSEPS = (S(N)-S(1)) * 1.0E-5
C
C---- set trailing edge point coordinates
      XTE = 0.5*(X(1) + X(N))
      YTE = 0.5*(Y(1) + Y(N))
C
C---- get first guess for SLE
      DO 10 I=3, N-2
        DXTE = X(I) - XTE
        DYTE = Y(I) - YTE
        DX = X(I+1) - X(I)
        DY = Y(I+1) - Y(I)
        DOTP = DXTE*DX + DYTE*DY
        IF(DOTP .LT. 0.0) GO TO 11
   10 CONTINUE
C
   11 SLE = S(I)
C
C---- check for sharp LE case
      IF(S(I) .EQ. S(I-1)) THEN
ccc        WRITE(*,*) 'Sharp LE found at ',I,SLE
        RETURN
      ENDIF
C
C---- Newton iteration to get exact SLE value
      DO 20 ITER=1, 50
        XLE  = SEVAL(SLE,X,XP,S,N)
        YLE  = SEVAL(SLE,Y,YP,S,N)
        DXDS = DEVAL(SLE,X,XP,S,N)
        DYDS = DEVAL(SLE,Y,YP,S,N)
        DXDD = D2VAL(SLE,X,XP,S,N)
        DYDD = D2VAL(SLE,Y,YP,S,N)
C
        XCHORD = XLE - XTE
        YCHORD = YLE - YTE
C
C------ drive dot product between chord line and LE tangent to zero
        RES  = XCHORD*DXDS + YCHORD*DYDS
        RESS = DXDS  *DXDS + DYDS  *DYDS
     &       + XCHORD*DXDD + YCHORD*DYDD
C
C------ Newton delta for SLE 
        DSLE = -RES/RESS
C
        DSLE = MAX( DSLE , -0.02*ABS(XCHORD+YCHORD) )
        DSLE = MIN( DSLE ,  0.02*ABS(XCHORD+YCHORD) )
        SLE = SLE + DSLE
        IF(ABS(DSLE) .LT. DSEPS) RETURN
   20 CONTINUE
      WRITE(*,*) 'LEFIND:  LE point not found.  Continuing...'
      SLE = S(I)
      RETURN
      END

      FUNCTION SEVAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates X(SS)                             |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      SEVAL = T*X(I) + (1.0-T)*X(I-1) + (T-T*T)*((1.0-T)*CX1 - T*CX2)
      RETURN
      END ! SEVAL

      FUNCTION DEVAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates dX/dS(SS)                         |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      DEVAL = X(I) - X(I-1) + (1.-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.)*CX2
      DEVAL = DEVAL/DS
      RETURN
      END ! DEVAL

      SUBROUTINE TRISOL(A,B,C,D,KK)
      DIMENSION A(KK),B(KK),C(KK),D(KK)
C-----------------------------------------
C     Solves KK long, tri-diagonal system |
C                                         |
C             A C          D              |
C             B A C        D              |
C               B A .      .              |
C                 . . C    .              |
C                   B A    D              |
C                                         |
C     The righthand side D is replaced by |
C     the solution.  A, C are destroyed.  |
C-----------------------------------------
C
      DO 1 K=2, KK
        KM = K-1
        C(KM) = C(KM) / A(KM)
        D(KM) = D(KM) / A(KM)
        A(K) = A(K) - B(K)*C(KM)
        D(K) = D(K) - B(K)*D(KM)
    1 CONTINUE
C
      D(KK) = D(KK)/A(KK)
C
      DO 2 K=KK-1, 1, -1
        D(K) = D(K) - C(K)*D(K+1)
    2 CONTINUE
C
      RETURN
      END ! TRISOL

      SUBROUTINE PANGEN(XNEW,YNEW,NPAN,XBIN,YBIN,NB,CVPAR,CTERAT,
     &                  CTRRAT,XSREF1,XSREF2,XPREF1,XPREF2)
C---------------------------------------------------
C     Set paneling distribution from buffer airfoil
C     geometry, thus creating current airfoil.
C 
C     If REFINE=True, bunch points at x=XSREF on
C     top side and at x=XPREF on bottom side
C     by setting a fictitious local curvature of
C     CTRRAT*(LE curvature) there.
C---------------------------------------------------
      PARAMETER (IQX=1000)
      PARAMETER (IBX=4*IQX)
      PARAMETER (IWX=IQX/8+2)
      PARAMETER (IZX=IQX+IWX)

      INTEGER, INTENT(IN) :: NPAN, NB
      REAL*8, DIMENSION(NB), INTENT(IN) :: XBIN, YBIN
      REAL*8, DIMENSION(NPAN), INTENT(OUT) :: XNEW, YNEW
      REAL*8, INTENT(IN) :: CVPAR, CTERAT, CTRRAT, XSREF1, XSREF2,
     &                      XPREF1, XPREF2
      REAL*8 W1(6*IQX),W2(6*IQX),W3(6*IQX),W4(6*IQX),
     &     W5(6*IQX),W6(6*IQX)
      REAL*8 SB(IBX), SNEW(5*IBX), S(IZX)
      REAL*8 XB(IBX), YB(IBX), XBP(IBX), YBP(IBX), X(IZX), Y(IZX)

C     Populate XB and YB arrays
      XB(1:NB) = XBIN
      YB(1:NB) = YBIN
C
      IF(NB.LT.2) THEN
       WRITE(*,*) 'PANGEN: Buffer airfoil not available.'
       N = 0
       RETURN
      ENDIF
C
C---- Number of temporary nodes for panel distribution calculation
C       exceeds the specified panel number by factor of IPFAC.
      IPFAC = 3
      IPFAC = 5
C
C---- number of airfoil panel points
      N = NPAN

cC---- number of wake points
c      NW = NPAN/8 + 2
c      IF(NW.GT.IWX) THEN
c       WRITE(*,*)
c     &  'Array size (IWX) too small.  Last wake point index reduced.'
c       NW = IWX
c      ENDIF
C
C---- set arc length spline parameter
      CALL SCALC(XB,YB,SB,NB)
C
C---- spline raw airfoil coordinates
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
C---- normalizing length (~ chord)
      SBREF = 0.5*(SB(NB)-SB(1))
C
C---- set up curvature array
      DO I = 1, NB
        W5(I) = ABS( CURV(SB(I),XB,XBP,YB,YBP,SB,NB) ) * SBREF
      ENDDO
C
C---- locate LE point arc length value and the normalized curvature there
      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
      CVLE = ABS( CURV(SBLE,XB,XBP,YB,YBP,SB,NB) ) * SBREF
C
C---- check for doubled point (sharp corner) at LE
      IBLE = 0
      DO I = 1, NB-1
        IF(SBLE.EQ.SB(I) .AND. SBLE.EQ.SB(I+1)) THEN
         IBLE = I
         WRITE(*,*)
         WRITE(*,*) 'Sharp leading edge'
         GO TO 21
        ENDIF
      ENDDO
 21   CONTINUE
C
C---- set LE, TE points
      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XBTE = 0.5*(XB(1)+XB(NB))
      YBTE = 0.5*(YB(1)+YB(NB))
      CHBSQ = (XBTE-XBLE)**2 + (YBTE-YBLE)**2
C
C---- set average curvature over 2*NK+1 points within Rcurv of LE point
      NK = 3
      CVSUM = 0.
      DO K = -NK, NK
        FRAC = FLOAT(K)/FLOAT(NK)
        SBK = SBLE + FRAC*SBREF/MAX(CVLE,20.0)
        CVK = ABS( CURV(SBK,XB,XBP,YB,YBP,SB,NB) ) * SBREF
        CVSUM = CVSUM + CVK
      ENDDO
      CVAVG = CVSUM/FLOAT(2*NK+1)
C
C---- dummy curvature for sharp LE
      IF(IBLE.NE.0) CVAVG = 10.0
C
C---- set curvature attraction coefficient actually used
      CC = 6.0 * CVPAR
C
C---- set artificial curvature at TE to bunch panels there
      CVTE = CVAVG * CTERAT
      W5(1)  = CVTE
      W5(NB) = CVTE
C
C
C**** smooth curvature array for smoother panel size distribution  ****
C
CCC      CALL ASKR('Enter curvature smoothing length/c^',SMOOL)
CCC      SMOOL = 0.010
C
C---- set smoothing length = 1 / averaged LE curvature, but 
C-    no more than 5% of chord and no less than 1/4 average panel spacing
      SMOOL = MAX( 1.0/MAX(CVAVG,20.0) , 0.25 /FLOAT(NPAN/2) )
C
      SMOOSQ = (SMOOL*SBREF) ** 2
C
C---- set up tri-diagonal system for smoothed curvatures
      W2(1) = 1.0
      W3(1) = 0.0
      DO I=2, NB-1
        DSM = SB(I) - SB(I-1)
        DSP = SB(I+1) - SB(I)
        DSO = 0.5*(SB(I+1) - SB(I-1))
C
        IF(DSM.EQ.0.0 .OR. DSP.EQ.0.0) THEN
C------- leave curvature at corner point unchanged
         W1(I) = 0.0
         W2(I) = 1.0
         W3(I) = 0.0
        ELSE
         W1(I) =  SMOOSQ * (         - 1.0/DSM) / DSO
         W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
        ENDIF
      ENDDO
C
      W1(NB) = 0.0
      W2(NB) = 1.0
C
C---- fix curvature at LE point by modifying equations adjacent to LE
      DO I=2, NB-1
        IF(SB(I).EQ.SBLE .OR. I.EQ.IBLE .OR. I.EQ.IBLE+1) THEN
C------- if node falls right on LE point, fix curvature there
         W1(I) = 0.
         W2(I) = 1.0
         W3(I) = 0.
         W5(I) = CVLE
        ELSE IF(SB(I-1).LT.SBLE .AND. SB(I).GT.SBLE) THEN
C------- modify equation at node just before LE point
         DSM = SB(I-1) - SB(I-2)
         DSP = SBLE    - SB(I-1)
         DSO = 0.5*(SBLE - SB(I-2))
C
         W1(I-1) =  SMOOSQ * (         - 1.0/DSM) / DSO
         W2(I-1) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I-1) =  0.
         W5(I-1) = W5(I-1) + SMOOSQ*CVLE/(DSP*DSO)
C
C------- modify equation at node just after LE point
         DSM = SB(I) - SBLE
         DSP = SB(I+1) - SB(I)
         DSO = 0.5*(SB(I+1) - SBLE)
         W1(I) =  0.
         W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
         W5(I) = W5(I) + SMOOSQ*CVLE/(DSM*DSO)
C
         GO TO 51
        ENDIF
      ENDDO
   51 CONTINUE
C
C---- set artificial curvature at bunching points and fix it there
      DO I=2, NB-1
C------ chord-based x/c coordinate
        XOC = (  (XB(I)-XBLE)*(XBTE-XBLE)
     &         + (YB(I)-YBLE)*(YBTE-YBLE) ) / CHBSQ
C
        IF(SB(I).LT.SBLE) THEN
C------- check if top side point is in refinement area
         IF(XOC.GT.XSREF1 .AND. XOC.LT.XSREF2) THEN
          W1(I) = 0.
          W2(I) = 1.0
          W3(I) = 0.
          W5(I) = CVLE*CTRRAT
         ENDIF
        ELSE
C------- check if bottom side point is in refinement area
         IF(XOC.GT.XPREF1 .AND. XOC.LT.XPREF2) THEN
          W1(I) = 0.
          W2(I) = 1.0
          W3(I) = 0.
          W5(I) = CVLE*CTRRAT
         ENDIF
        ENDIF
      ENDDO
C
C---- solve for smoothed curvature array W5
      IF(IBLE.EQ.0) THEN
       CALL TRISOL(W2,W1,W3,W5,NB)
      ELSE
       I = 1
       CALL TRISOL(W2(I),W1(I),W3(I),W5(I),IBLE)
       I = IBLE+1
       CALL TRISOL(W2(I),W1(I),W3(I),W5(I),NB-IBLE)
      ENDIF
C
C---- find max curvature
      CVMAX = 0.
      DO I=1, NB
        CVMAX = MAX( CVMAX , ABS(W5(I)) )
      ENDDO
C
C---- normalize curvature array
      DO I=1, NB
        W5(I) = W5(I) / CVMAX
      ENDDO
C
C---- spline curvature array
      CALL SEGSPL(W5,W6,SB,NB)
C
C---- Set initial guess for node positions uniform in s.
C     More nodes than specified (by factor of IPFAC) are 
C     temporarily used  for more reliable convergence.
      NN = IPFAC*(N-1)+1
C
C---- ratio of lengths of panel at TE to one away from the TE
      RDSTE = 0.667
      RTF = (RDSTE-1.0)*2.0 + 1.0
C
      IF(IBLE.EQ.0) THEN
C
       DSAVG = (SB(NB)-SB(1))/(FLOAT(NN-3) + 2.0*RTF)
       SNEW(1) = SB(1)
       DO I=2, NN-1
         SNEW(I) = SB(1) + DSAVG * (FLOAT(I-2) + RTF)
       ENDDO
       SNEW(NN) = SB(NB)
C      DP mod: to suppress compiler warning
       NN1 = 0     
C
      ELSE
C
       NFRAC1 = (N * IBLE) / NB
C
       NN1 = IPFAC*(NFRAC1-1)+1
       DSAVG1 = (SBLE-SB(1))/(FLOAT(NN1-2) + RTF)
       SNEW(1) = SB(1)
       DO I=2, NN1
         SNEW(I) = SB(1) + DSAVG1 * (FLOAT(I-2) + RTF)
       ENDDO
C
       NN2 = NN - NN1 + 1
       DSAVG2 = (SB(NB)-SBLE)/(FLOAT(NN2-2) + RTF)
       DO I=2, NN2-1
         SNEW(I-1+NN1) = SBLE + DSAVG2 * (FLOAT(I-2) + RTF)
       ENDDO
       SNEW(NN) = SB(NB)
C
      ENDIF
C
C---- Newton iteration loop for new node positions
      DO 10 ITER=1, 20
C
C------ set up tri-diagonal system for node position deltas
        CV1  = SEVAL(SNEW(1),W5,W6,SB,NB)
        CV2  = SEVAL(SNEW(2),W5,W6,SB,NB)
        CVS1 = DEVAL(SNEW(1),W5,W6,SB,NB)
        CVS2 = DEVAL(SNEW(2),W5,W6,SB,NB)
C
        CAVM = SQRT(CV1**2 + CV2**2)
        IF(CAVM .EQ. 0.0) THEN
          CAVM_S1 = 0.
          CAVM_S2 = 0.
        ELSE
          CAVM_S1 = CVS1 * CV1/CAVM
          CAVM_S2 = CVS2 * CV2/CAVM
        ENDIF
C
        DO 110 I=2, NN-1
          DSM = SNEW(I) - SNEW(I-1)
          DSP = SNEW(I) - SNEW(I+1)
          CV3  = SEVAL(SNEW(I+1),W5,W6,SB,NB)
          CVS3 = DEVAL(SNEW(I+1),W5,W6,SB,NB)
C
          CAVP = SQRT(CV3**2 + CV2**2)
          IF(CAVP .EQ. 0.0) THEN
            CAVP_S2 = 0.
            CAVP_S3 = 0.
          ELSE
            CAVP_S2 = CVS2 * CV2/CAVP
            CAVP_S3 = CVS3 * CV3/CAVP
          ENDIF
C
          FM = CC*CAVM + 1.0
          FP = CC*CAVP + 1.0
C
          REZ = DSP*FP + DSM*FM
C
C-------- lower, main, and upper diagonals
          W1(I) =      -FM  +  CC*               DSM*CAVM_S1
          W2(I) =  FP + FM  +  CC*(DSP*CAVP_S2 + DSM*CAVM_S2)
          W3(I) = -FP       +  CC* DSP*CAVP_S3
C
C-------- residual, requiring that
C         (1 + C*curv)*deltaS is equal on both sides of node i
          W4(I) = -REZ
C
          CV1 = CV2
          CV2 = CV3
          CVS1 = CVS2
          CVS2 = CVS3
          CAVM    = CAVP
          CAVM_S1 = CAVP_S2
          CAVM_S2 = CAVP_S3
  110   CONTINUE
C
C------ fix endpoints (at TE)
        W2(1) = 1.0
        W3(1) = 0.0
        W4(1) = 0.0
        W1(NN) = 0.0
        W2(NN) = 1.0
        W4(NN) = 0.0
C
        IF(RTF .NE. 1.0) THEN
C------- fudge equations adjacent to TE to get TE panel length ratio RTF
C
         I = 2
         W4(I) = -((SNEW(I) - SNEW(I-1)) + RTF*(SNEW(I) - SNEW(I+1)))
         W1(I) = -1.0
         W2(I) =  1.0 + RTF
         W3(I) =      - RTF
C
         I = NN-1
         W4(I) = -((SNEW(I) - SNEW(I+1)) + RTF*(SNEW(I) - SNEW(I-1)))
         W3(I) = -1.0
         W2(I) =  1.0 + RTF
         W1(I) =      - RTF
        ENDIF
C
C
C------ fix sharp LE point
        IF(IBLE.NE.0) THEN
         I = NN1
         W1(I) = 0.0
         W2(I) = 1.0
         W3(I) = 0.0
         W4(I) = SBLE - SNEW(I)
        ENDIF
C
C------ solve for changes W4 in node position arc length values
        CALL TRISOL(W2,W1,W3,W4,NN)
C
C------ find under-relaxation factor to keep nodes from changing order
        RLX = 1.0
        DMAX = 0.0
        DO I=1, NN-1
          DS  = SNEW(I+1) - SNEW(I)
          DDS = W4(I+1) - W4(I)
          DSRAT = 1.0 + RLX*DDS/DS
          IF(DSRAT.GT.4.0) RLX = (4.0-1.0)*DS/DDS
          IF(DSRAT.LT.0.2) RLX = (0.2-1.0)*DS/DDS
          DMAX = MAX(ABS(W4(I)),DMAX)
        ENDDO
C
C------ update node position
        DO I=2, NN-1
          SNEW(I) = SNEW(I) + RLX*W4(I)
        ENDDO
C
CCC        IF(RLX.EQ.1.0) WRITE(*,*) DMAX
CCC        IF(RLX.NE.1.0) WRITE(*,*) DMAX,'    RLX =',RLX
        IF(ABS(DMAX).LT.1.E-3) GO TO 11
   10 CONTINUE
      WRITE(*,*) 'Paneling convergence failed.  Continuing anyway...'
C
   11 CONTINUE
C
C---- set new panel node coordinates
      DO I=1, N
        IND = IPFAC*(I-1) + 1
        S(I) = SNEW(IND)
        X(I) = SEVAL(SNEW(IND),XB,XBP,SB,NB)
        Y(I) = SEVAL(SNEW(IND),YB,YBP,SB,NB)
      ENDDO
C
C
C---- go over buffer airfoil again, checking for corners (double points)
      NCORN = 0
      DO 25 IB=1, NB-1
        IF(SB(IB) .EQ. SB(IB+1)) THEN
C------- found one !
C
         NCORN = NCORN+1
         XBCORN = XB(IB)
         YBCORN = YB(IB)
         SBCORN = SB(IB)
C
C------- find current-airfoil panel which contains corner
         DO 252 I=1, N
C
C--------- keep stepping until first node past corner
           IF(S(I) .LE. SBCORN) GO TO 252
C
C---------- move remainder of panel nodes to make room for additional node
            DO 2522 J=N, I, -1
              X(J+1) = X(J)
              Y(J+1) = Y(J)
              S(J+1) = S(J)
 2522       CONTINUE
            N = N+1
C
            IF(N .GT. IQX-1)
     &       STOP 'PANEL: Too many panels. Increase IQX in XFOIL.INC'
C
            X(I) = XBCORN
            Y(I) = YBCORN
            S(I) = SBCORN
C
C---------- shift nodes adjacent to corner to keep panel sizes comparable
            IF(I-2 .GE. 1) THEN
             S(I-1) = 0.5*(S(I) + S(I-2))
             X(I-1) = SEVAL(S(I-1),XB,XBP,SB,NB)
             Y(I-1) = SEVAL(S(I-1),YB,YBP,SB,NB)
            ENDIF
C
            IF(I+2 .LE. N) THEN
             S(I+1) = 0.5*(S(I) + S(I+2))
             X(I+1) = SEVAL(S(I+1),XB,XBP,SB,NB)
             Y(I+1) = SEVAL(S(I+1),YB,YBP,SB,NB)
            ENDIF
C
C---------- go on to next input geometry point to check for corner
            GO TO 25
C
  252    CONTINUE
        ENDIF
   25 CONTINUE

C     Output new geometry
      XNEW(1:NPAN) = X(1:NPAN)
      YNEW(1:NPAN) = Y(1:NPAN)
C
      RETURN
      END ! PANGEN
