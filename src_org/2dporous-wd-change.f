      PARAMETER(MX=1501,MY=1001,MZ=3,NNN=50000,MSTEP=320000)
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      DIMENSION UUIN(MSTEP),VVIN(MSTEP)
      DIMENSION CD(NNN),CL(NNN),TIM(NNN)
      DIMENSION PP(0:MX+2,0:MY+2),RHS(MX,MY)
      DIMENSION UU(0:MX+2,0:MY+2),VV(0:MX+2,0:MY+2)
      DIMENSION UO(0:MX+2,0:MY+2),VO(0:MX+2,0:MY+2)
      DIMENSION UAVE(MX,MY),VAVE(MX,MY),PAVE(MX,MY)
      DIMENSION U(MX,MY),V(MX,MY),P(MX,MY),ZETA(MX,MY),PSI(MX,MY)
C     DIMENSION HU(MX,MY),HV(MX,MY),HP(MX,MY),HS(MX,MY),HZ(MX,MY)
      DIMENSION HUAVE(MX,MY,MZ),HVAVE(MX,MY,MZ),HWAVE(MX,MY,MZ)
      DIMENSION HS(MX,MY,MZ),HZ(MX,MY,MZ),HP(MX,MY,MZ),A(2)
      DIMENSION HU(MX,MY,MZ),HV(MX,MY,MZ),HW(MX,MY,MZ)
      DIMENSION X(MX,MY),Y(MX,MY),IFLAG(MX,MY)
      DIMENSION DFUNC(MX,MY),PFUNC(MY)
C
      OPEN(10,FILE='3d-inst-vis.dat',FORM='UNFORMATTED')
C
C....<<SOR METHOD PARAMETERS>>
      OMEGA=.9D0
      EPS=1.E-3
      ITER=100
C
C....<<TIME STEP ETC>>
      ICON=0
      IPROG=800
      NSTEP=320000
      NSOR=100
      NQ=50
      TIME=0.D0
      DT=1.E-3
C
C....<<REYNOLDS NUMBER>>
CC    RE=1.E+2
CC    RE=1.E+3
      RE=1.E+4
      REI=1.D0/RE
C
C....<<DIVIDED LENGTH>>
      NX=1501
      NY=1001
      NZ=3
      NOXY=(NX-1)*(NY-1)
      XLEN=15.D0
      YLEN=10.D0
      DX=XLEN/DFLOAT(NX-1)
      DY=YLEN/DFLOAT(NY-1)
      DXI=1.D0/DX
      DYI=1.D0/DY
      DDXI=DXI/DX
      DDYI=DYI/DY
C
C....<<MAKING 2D-GRIDS>>
      DO J=1,NY
      DO I=1,NX
      X(I,J)=XLEN/DFLOAT(NX-1)*DFLOAT(I-1)-5.D0
      Y(I,J)=YLEN/DFLOAT(NY-1)*DFLOAT(J-1)-5.D0
      END DO
      END DO
C
C....<<POROUS DISK MODEL>>
      IA=500
      IB=502
CCC   JA=451
CCC   JB=551
      JA=426
      JB=576
C
      CDD=13.D0
C
      OPEN(1,FILE='pfunc.dat')
      READ(1,'(F15.7)')(PFUNC(J),J=JA,JB)
      CLOSE(1)
      DO J=1,NY
      DO I=1,NX
      DFUNC(I,J)=0.D0
      END DO
      END DO
      DO J=JA,JB
      DO I=IA,IB
      DFUNC(I,J)=PFUNC(J)*CDD
      END DO
      END DO
C
C....<<INFLOW DATA>>
      OPEN(1,FILE='time-hist-inflow_5.0deg.dat',FORM='FORMATTED')
      READ(1,'(2F20.15)')(UUIN(ISTEP),VVIN(ISTEP),ISTEP=1,NSTEP)
      CLOSE(1)
c     DO ISTEP=1,NSTEP
c     UUIN(ISTEP)=1.D0
c     VVIN(ISTEP)=0.D0
c     END DO
C
C....<<ZERO CLEAR>>
      DO J=1,NY
      DO I=1,NX
      UU(I,J)=UUIN(1)
      VV(I,J)=VVIN(1)
      PP(I,J)=0.D0
      UAVE(I,J)=0.D0
      VAVE(I,J)=0.D0
      PAVE(I,J)=0.D0
      END DO
      END DO
C
C....<<TO BE CONTINUED>>
      IF(ICON.EQ.1)THEN
      OPEN(1,FILE='uvps.dat')
      READ(1,*)TIME,RE,NX,NY
      READ(1,*)((UU(I,J),I=0,NX+2),J=0,NY+2)
      READ(1,*)((VV(I,J),I=0,NX+2),J=0,NY+2)
      READ(1,*)((PP(I,J),I=0,NX+2),J=0,NY+2)
      CLOSE(1)
      END IF
C
C....<<TIME EVOLUTION LOOP>>
      ICOUNT=0
      DO 1000 ISTEP=1,NSTEP
      TIME=TIME+DT
C    ===========================
C    ==Solve Navier-Stokes eq.==
C    ===========================
C....<<STORE PREVIOUS VALUE>>
      DO J=0,NY+2
      DO I=0,NX+2
      UO(I,J)=UU(I,J)
      VO(I,J)=VV(I,J)
      END DO
      END DO
C
C....<<FOR U>>
      DO J=2,NY
      DO I=2,NX-1
C....<<CONVECTIVE TERM(4TH-ORDER CENTRAL SCHEME)>>
      CU1=-UO(I+1,J)+9.D0*(UO(I,J)+UO(I-1,J))-UO(I-2,J)
      CU2=-UO(I+2,J)+9.D0*(UO(I+1,J)+UO(I,J))-UO(I-1,J)
      CV1=-VO(I+2,J-1)+9.D0*(VO(I+1,J-1)+VO(I,J-1))-VO(I-1,J-1)
      CV2=-VO(I+2,J)+9.D0*(VO(I+1,J)+VO(I,J))-VO(I-1,J)
      UU1=-UO(I+1,J)+27.D0*(UO(I,J)-UO(I-1,J))+UO(I-2,J)
      UU2=-UO(I+2,J)+27.D0*(UO(I+1,J)-UO(I,J))+UO(I-1,J)
      UU3=-UO(I,J+1)+27.D0*(UO(I,J)-UO(I,J-1))+UO(I,J-2)
      UU4=-UO(I,J+2)+27.D0*(UO(I,J+1)-UO(I,J))+UO(I,J-1)
      AX1=(CU1/16.D0)*(UU1/24.D0*DXI)
      AX2=(CU2/16.D0)*(UU2/24.D0*DXI)
      AY1=(CV1/16.D0)*(UU3/24.D0*DYI)
      AY2=(CV2/16.D0)*(UU4/24.D0*DYI)
C....<<4TH-ORDER NUMERICAL VISCOSITY(K-K SCHEME TYPE)>>
      UIJK=UO(I,J)*DXI
      VIJK=.25D0*(VO(I+1,J-1)+VO(I,J-1)+VO(I+1,J)+VO(I,J))*DYI
      U2=UO(I+2,J)-4.D0*(UO(I+1,J)+UO(I-1,J))+UO(I-2,J)+6.D0*UO(I,J)
      U4=UO(I,J+2)-4.D0*(UO(I,J+1)+UO(I,J-1))+UO(I,J-2)+6.D0*UO(I,J)
      ADVX=.5D0*(AX1+AX2+AY1+AY2)+(DABS(UIJK)*U2+DABS(VIJK)*U4)/4.D0
C....<<VISCOUS TERM(2ND-ORDER CENTRAL SCHEME)>>
      VISX=((UO(I+1,J)-2.D0*UO(I,J)+UO(I-1,J))*DDXI
     &     +(UO(I,J+1)-2.D0*UO(I,J)+UO(I,J-1))*DDYI)*REI
C....<<POROUS DISK MODEL>>
      DVEL=DSQRT(UO(I,J)**2+VO(I,J)**2)
      SDRAGX=-UO(I,J)*DVEL*DFUNC(I,J)
C....<<INTERMEDIATE VELOCITY>>
      UU(I,J)=UO(I,J)+DT*(-ADVX+VISX+SDRAGX)
      END DO
      END DO
C
C....<<FOR V>>
      DO J=2,NY-1
      DO I=2,NX
C....<<CONVECTIVE TERM(4TH-ORDER CENTRAL SCHEME)>>
      CU1=-UO(I-1,J+2)+9.D0*(UO(I-1,J+1)+UO(I-1,J))-UO(I-1,J-1)
      CU2=-UO(I,J+2)+9.D0*(UO(I,J+1)+UO(I,J))-UO(I,J-1)
      CV1=-VO(I,J+1)+9.D0*(VO(I,J)+VO(I,J-1))-VO(I,J-2)
      CV2=-VO(I,J+2)+9.D0*(VO(I,J+1)+VO(I,J))-VO(I,J-1)
      VV1=-VO(I+1,J)+27.D0*(VO(I,J)-VO(I-1,J))+VO(I-2,J)
      VV2=-VO(I+2,J)+27.D0*(VO(I+1,J)-VO(I,J))+VO(I-1,J)
      VV3=-VO(I,J+1)+27.D0*(VO(I,J)-VO(I,J-1))+VO(I,J-2)
      VV4=-VO(I,J+2)+27.D0*(VO(I,J+1)-VO(I,J))+VO(I,J-1)
      AX1=(CU1/16.D0)*(VV1/24.D0*DXI)
      AX2=(CU2/16.D0)*(VV2/24.D0*DXI)
      AY1=(CV1/16.D0)*(VV3/24.D0*DYI)
      AY2=(CV2/16.D0)*(VV4/24.D0*DYI)
C....<<4TH-ORDER NUMERICAL VISCOSITY(K-K SCHEME TYPE)>>
      UIJK=.25D0*(UO(I-1,J+1)+UO(I-1,J)+UO(I,J+1)+UO(I,J))*DXI
      VIJK=VO(I,J)*DYI
      V2=VO(I+2,J)-4.D0*(VO(I+1,J)+VO(I-1,J))+VO(I-2,J)+6.D0*VO(I,J)
      V4=VO(I,J+2)-4.D0*(VO(I,J+1)+VO(I,J-1))+VO(I,J-2)+6.D0*VO(I,J)
      ADVY=.5D0*(AX1+AX2+AY1+AY2)+(DABS(UIJK)*V2+DABS(VIJK)*V4)/4.D0
C....<<VISCOUS TERM(2ND-ORDER CENTRAL SCHEME)>>
      VISY=((VO(I+1,J)-2.D0*VO(I,J)+VO(I-1,J))*DDXI
     &     +(VO(I,J+1)-2.D0*VO(I,J)+VO(I,J-1))*DDYI)*REI
C....<<POROUS DISK MODEL>>
      DVEL=DSQRT(UO(I,J)**2+VO(I,J)**2)
      SDRAGY=-VO(I,J)*DVEL*DFUNC(I,J)
C....<<INTERMEDIATE VELOCITY>>
      VV(I,J)=VO(I,J)+DT*(-ADVY+VISY+SDRAGY)
      END DO
      END DO
C
C    ==================================================
C    ==Solve Poisson eq.for Pressure and New Velocity==
C    ==================================================
C....<<RIGHT HAND SIDE>>
      DO J=2,NY
      DO I=2,NX
      RHS(I,J)=-((UU(I,J)-UU(I-1,J))*DXI
     &          +(VV(I,J)-VV(I,J-1))*DYI)/DT
      END DO
      END DO
C
C....<<SOR METHOD LOOP>>
      DO LOOP=1,ITER
      ERR=0.D0
C....<<INTERNAL POINTS>>
      DO J=2,NY
      DO I=2,NX
C....<<POISSON EQ.>>
      PAGE=(PP(I+1,J)+PP(I-1,J))*DDXI
     &    +(PP(I,J+1)+PP(I,J-1))*DDYI+RHS(I,J)
C....<<NEW PRESSURE>>
      DP=.5D0*PAGE/(DDXI+DDYI)-PP(I,J)
      ERR=ERR+DP*DP
      PP(I,J)=PP(I,J)+OMEGA*DP
      END DO
      END DO
C....<<CONVERGE OR NOT ?>>
      ERR=DSQRT(ERR/DFLOAT(NOXY))
      IF(ERR.LT.EPS)GOTO 99
      END DO
   99 CONTINUE
C
C....<<AVERAGING & SUBTRACTION>>
      PPP=0.D0
      DO J=1,NY+1
      DO I=1,NX+1
      PPP=PPP+PP(I,J)
      END DO
      END DO
      DO J=1,NY+1
      DO I=1,NX+1
      PP(I,J)=PP(I,J)-PPP/DFLOAT((NX+1)*(NY+1))
      END DO
      END DO
C
C....<<MODIFY PRESSURE B.C.>>
      DO I=1,NX+1
      PP(I,1)=PP(I,2)
      PP(I,NY+1)=PP(I,NY)
      END DO
      DO J=1,NY+1
      PP(1,J)=PP(2,J)
      PP(NX+1,J)=PP(NX,J)
      END DO
C
C....<<NEW VELOCITY>>
      DO J=2,NY
      DO I=2,NX-1
      UU(I,J)=UU(I,J)-DT*(PP(I+1,J)-PP(I,J))*DXI
      END DO
      END DO
      DO J=2,NY-1
      DO I=2,NX
      VV(I,J)=VV(I,J)-DT*(PP(I,J+1)-PP(I,J))*DYI
      END DO
      END DO
C
C....<<MODIFY VELOCITY B.C.>>
      DO J=2,NY
      UU(1,J)=UUIN(ISTEP)
      UU(0,J)=UUIN(ISTEP)
      VV(1,J)=VVIN(ISTEP)
      VV(0,J)=VVIN(ISTEP)
      UU(NX,J)=UO(NX,J)-DT*(UO(NX,J)-UO(NX-1,J))*DXI
      VV(NX+1,J)=VV(NX,J)
      VV(NX+2,J)=2.D0*VV(NX+1,J)-VV(NX,J)
      UU(NX+1,J)=2.D0*UU(NX,J)-UU(NX-1,J)
      END DO
      DO I=0,NX+2
      UU(I,1)=UUIN(ISTEP)
      UU(I,0)=UUIN(ISTEP)
      VV(I,1)=VVIN(ISTEP)
      VV(I,0)=VVIN(ISTEP)
      UU(I,NY+1)=UUIN(ISTEP)
      UU(I,NY+2)=UUIN(ISTEP)
      VV(I,NY)=VVIN(ISTEP)
      VV(I,NY+1)=VVIN(ISTEP)
      END DO
C
C....<<REGULAR ARRANGEMENT>>
      DO J=1,NY
      DO I=1,NX
      U(I,J)=.5D0*(UU(I,J)+UU(I,J+1))
      V(I,J)=.5D0*(VV(I,J)+VV(I+1,J))
CCC   P(I,J)=.25D0*(PP(I,J)+PP(I+1,J)+PP(I+1,J+1)+PP(I,J+1))
      END DO
      END DO
C
C....<<VORTICITY>>
C     DO J=2,NY-1
C     DO I=2,NX-1
C     ZETA(I,J)=.5D0*(V(I+1,J)-V(I-1,J))*DXI
C    &         -.5D0*(U(I,J+1)-U(I,J-1))*DYI
C     END DO
C     END DO
C     DO I=1,NX
C     ZETA(I,1)=2.D0*ZETA(I,2)-ZETA(I,3)
C     ZETA(I,NY)=2.D0*ZETA(I,NY-1)-ZETA(I,NY-2)
C     END DO
C     DO J=1,NY
C     ZETA(1,J)=2.D0*ZETA(2,J)-ZETA(3,J)
C     ZETA(NX,J)=2.D0*ZETA(NX-1,J)-ZETA(NX-2,J)
C     END DO
C
C....<<STREAM FUCTION>>
C     DDY=DY*(NY-1)
C     DO J=1,NY
C     PSI(1,J)=DDY/(NY-1)*(J-1)-DDY*.5D0
C     END DO
C     DO I=2,NX
C     DO J=1,NY
C     PSI(I,J)=PSI(I-1,J)-DX*V(I,J)
C     END DO
C     END DO
C
C....<<SUMMATION>>
      DO J=1,NY
      DO I=1,NX
      UAVE(I,J)=UAVE(I,J)+U(I,J)
      VAVE(I,J)=VAVE(I,J)+V(I,J)
CCC   PAVE(I,J)=PAVE(I,J)+P(I,J)
      END DO
      END DO
C
C....<<DISPLAY>>
      IF(MOD(ISTEP,NSOR).EQ.0)THEN
      WRITE(*,*)'================================================='
      WRITE(*,700)CDD,TIME,ISTEP
  700 FORMAT(' Cd=',F5.1,' Time=',F10.5,'  Istep=',I7)
      WRITE(*,710)LOOP
  710 FORMAT(' Poisson equation iteration=',I7)
      WRITE(*,*)'Root-Mean-Square error=',ERR
      END IF
C
C....<<FLOW VISUALIZATION>>
      IF(MOD(ISTEP,IPROG).EQ.0)THEN
      ICOUNT=ICOUNT+1
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
      HU(I,J,K)=U(I,J)
      HV(I,J,K)=V(I,J)
      HW(I,J,K)=0.
CCC   HP(I,J,K)=P(I,J)
CCC   HS(I,J,K)=PSI(I,J)
CCC   HZ(I,J,K)=ZETA(I,J)
      END DO
      END DO
      END DO
      WRITE(10)ISTEP,SNGL(TIME)
      WRITE(10)(SNGL(A(II)),II=1,2)
      WRITE(10)(((HU(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      WRITE(10)(((HV(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      WRITE(10)(((HW(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
CCC   WRITE(10)(((HP(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
CCC   WRITE(10)(((HS(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
CCC   WRITE(10)(((HZ(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      WRITE(*,*)'RC-Scope date=',ICOUNT,'/',NSTEP/IPROG,'(Total)'
      END IF
 1000 CONTINUE
      CLOSE(10)
C
C....<<TIME AVERAGING>>
      DO J=1,NY
      DO I=1,NX
      UAVE(I,J)=UAVE(I,J)/DFLOAT(NSTEP)
      VAVE(I,J)=VAVE(I,J)/DFLOAT(NSTEP)
CCC   PAVE(I,J)=PAVE(I,J)/DFLOAT(NSTEP)
      END DO
      END DO
      OPEN(1,FILE='5D-10D-U-profile.dat',FORM='FORMATTED')
CC    WRITE(1,'(2F15.7)')(UAVE(1001,J),UAVE(1501,J),J=301,701)
      WRITE(1,'(2F15.7)')(UAVE(1001,J),UAVE(1501,J),J=1,NY)
      CLOSE(1)
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
      HUAVE(I,J,K)=UAVE(I,J)
      HVAVE(I,J,K)=VAVE(I,J)
      HWAVE(I,J,K)=0.
      END DO
      END DO
      END DO
      open(unit=11,file='3d-tave-vis.dat',form='unformatted')
      write(11)ISTEP,SNGL(TIME)
      write(11)(SNGL(A(II)),II=1,2)
      write(11)(((HUAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      write(11)(((HVAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      write(11)(((HWAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      close(11)
C
C....<<OUTPUT>>
C     OPEN(1,FILE='uvp-s.dat')
C     WRITE(1,*)TIME,RE,NX,NY
C     WRITE(1,720)((UU(I,J),I=0,NX+2),J=0,NY+2)
C     WRITE(1,720)((VV(I,J),I=0,NX+2),J=0,NY+2)
C     WRITE(1,720)((PP(I,J),I=0,NX+2),J=0,NY+2)
C 720 FORMAT(3F15.7)
C     CLOSE(1)
C     OPEN(2,FILE='uvp-r.dat')
C     WRITE(2,*)TIME,RE,NX,NY
C     WRITE(2,730)((U(I,J),I=1,NX),J=1,NY)
C     WRITE(2,730)((V(I,J),I=1,NX),J=1,NY)
C     WRITE(2,730)((P(I,J),I=1,NX),J=1,NY)
C 730 FORMAT(3F15.7)
C     CLOSE(2)
C     OPEN(3,FILE='2dave.dat')
C     WRITE(3,*)NX,NY
C     WRITE(3,740)((UAVE(I,J),I=1,NX),J=1,NY)
C     WRITE(3,740)((VAVE(I,J),I=1,NX),J=1,NY)
C     WRITE(3,740)((PAVE(I,J),I=1,NX),J=1,NY)
C740 FORMAT(3F15.7)
C     CLOSE(3)
C     OPEN(4,FILE='cd.dat')
C     WRITE(4,*)NCP
C     WRITE(4,750)(TIM(I),CD(I),CL(I),I=1,NCP)
C 750 FORMAT(3F15.7)
C     CLOSE(4)
C
      STOP
      END
