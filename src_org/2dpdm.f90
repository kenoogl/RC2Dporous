program main
!$ use omp_lib

implicit none
integer, PARAMETER :: NX=1501, NY=1001, NZ=3, NSTEP=320000, nWatch=23

double precision, dimension(NX,NY)         :: RHS, UAVE, VAVE, U, V, P, X, Y, DFUNC
double precision, dimension(0:NX+2,0:NY+2) :: PP, UU, VV, UO, VO
real,             dimension(NX,NY,NZ)      :: HUAVE, HVAVE, HWAVE, HU, HV, HW, HP
real,             dimension(3,NX,NY,NZ)    :: HUVW
double precision, dimension(NSTEP)         :: UUIN, VVIN
double precision, dimension(NY)            :: PFUNC
double precision, dimension(2)             :: A

double precision, dimension(4,nWatch)      :: stop_watch
character(20)                              :: tm_label(nWatch)

double precision :: omega, eps, time, dt, RE, REI, ERR, CDD
double precision :: DX, DY, XLEN, YLEN, flop
integer          :: ITER, ICON, IPROG, NSOR, ICOUNT
integer          :: NOXY, ISTEP, LOOP
integer          :: IA, IB, JA, JB

!========================================================


call init_watch(nx, ny, nWatch, stop_watch, tm_label)


! Parameters

! SOR
OMEGA=.9D0
EPS  =1.E-3
ITER =100

! TIME STEP ETC
ICON  = 0      ! initial start=0, restart=1
IPROG = 800    ! Interval for vis. data output
NSOR  = 100    ! Interval for Display
TIME  = 0.D0
DT    = 1.E-3

! REYNOLDS NUMBER
RE=1.E+4
REI=1.D0/RE

! DIVIDED LENGTH
NOXY=(NX-1)*(NY-1)
XLEN= 15.D0
YLEN= 10.D0
DX  = XLEN/DFLOAT(NX-1)
DY  = YLEN/DFLOAT(NY-1)

! POROUS DISK MODEL

IA=500
IB=502
JA=426
JB=576
CDD=13.D0


call begin_watch(1, nWatch, stop_watch)
call gen_grid(nx, ny, xlen, ylen, X, Y, flop)
call end_watch(1, nWatch, stop_watch)


call begin_watch(2, nWatch, stop_watch)
call PDM_profile(nx, ny, IA, IB, JA, JB, CDD, PFUNC, DFUNC, flop)
call end_watch(2, nWatch, stop_watch)


call begin_watch(3, nWatch, stop_watch)
call Inflow_profile(NSTEP, UUIN, VVIN, flop)
call end_watch(3, nWatch, stop_watch)


call begin_watch(4, nWatch, stop_watch)
call init_Arrays(nx, ny, nz, nstep, UU, VV, PP, UAVE, VAVE, UUIN, VVIN, HW, HWAVE, flop)
call end_watch(4, nWatch, stop_watch)


! Restart
IF (ICON==1)THEN
  call begin_watch(5, nWatch, stop_watch)
  call read_rst_Data(nx, ny, time, UU, VV, PP, flop)
  call end_watch(5, nWatch, stop_watch)
END IF


!========================================================


OPEN(10,FILE='3d-inst-vis.dat',FORM='UNFORMATTED')
ICOUNT = 0

DO ISTEP=1, NSTEP
  TIME=TIME+DT

  call begin_watch(6, nWatch, stop_watch)
  !call store_vel(nx, ny, UO, VO, UU, VV, flop)
  !call store_vel2(nx, ny, UO, VO, UU, VV, flop)
  call store_vel3(nx, ny, UO, VO, UU, VV, flop)
  call end_watch(6, nWatch, stop_watch)


  call begin_watch(7, nWatch, stop_watch)
  call uflux(nx, ny, dt, rei, dx, dy, UO, VO, UU, DFUNC, flop)
  call end_watch(7, nWatch, stop_watch)


  call begin_watch(8, nWatch, stop_watch)
  call vflux(nx, ny, dt, rei, dx, dy, UO, VO, VV, DFUNC, flop)
  call end_watch(8, nWatch, stop_watch)


  call begin_watch(9, nWatch, stop_watch)
  call Poisson_RHS(nx, ny, dt, dx, dy, RHS, UU, VV, flop)
  call end_watch(9, nWatch, stop_watch)

!--------------------------------------------------------


  DO LOOP=1,ITER
    ERR=0.D0

    call begin_watch(10, nWatch, stop_watch)
    !call Poisson_AXB(nx, ny, dx, dy, omega, err, RHS, PP, flop)
    call Poisson_AXB2C(nx, ny, dx, dy, omega, err, RHS, PP, flop)
    call end_watch(10, nWatch, stop_watch)

    ERR=DSQRT(ERR/DFLOAT(NOXY))
    IF (ERR < EPS) exit
  END DO
  

  call begin_watch(11, nWatch, stop_watch)
  call Shift_Prs(nx, ny, PP, flop)
  call end_watch(11, nWatch, stop_watch)


  call begin_watch(12, nWatch, stop_watch)
  call BC_Prs(nx, ny, PP, flop)
  call end_watch(12, nWatch, stop_watch)

  call begin_watch(13, nWatch, stop_watch)
  call Prj_Vel(nx, ny, dt, dx, dy, UU, VV, PP, flop)
  call end_watch(13, nWatch, stop_watch)


  call begin_watch(14, nWatch, stop_watch)
  call BC_Vel(nx, ny, dt, dx, dy, nstep, istep, UU, VV, UO, UUIN, VVIN, flop)
  call end_watch(14, nWatch, stop_watch)
  

!--------------------------------------------------------


  call begin_watch(15, nWatch, stop_watch)
  call stg2CC(nx, ny, U, V, UU, VV, flop)
  call end_watch(15, nWatch, stop_watch)


  call begin_watch(16, nWatch, stop_watch)
  call accumVars(nx, ny, U, V, UAVE, VAVE, flop)
  call end_watch(16, nWatch, stop_watch)


! DISPLAY
  IF (MOD(ISTEP,NSOR)==0) THEN
    call begin_watch(17, nWatch, stop_watch)
    print "('=================================================')"
    print "(' Cd=',F5.1,' Time=',F10.5,'  Istep=',I7)", CDD, TIME, ISTEP
    print "(' Poisson equation iteration=',I7)", LOOP
    print "(' Root-Mean-Square error=',E14.7)", ERR
    call end_watch(17, nWatch, stop_watch)
  END IF
  

! FLOW VISUALIZATION
  IF (MOD(ISTEP,IPROG)==0) THEN

    !call begin_watch(18, nWatch, stop_watch)
    !call wrt_RCSIns(nx, ny, nz, istep, time, A, U, V, HU, HV, HW, flop)
    !call end_watch(18, nWatch, stop_watch)

    call begin_watch(19, nWatch, stop_watch)
    call wrt_SPHIns(nx, ny, nz, dx, dy, istep, time, U, V, PP, HUVW, HP, flop)
    call end_watch(19, nWatch, stop_watch)

    ICOUNT=ICOUNT+1
    print "('RC-Scope date=',I10, '/', I10,'(Total)')", ICOUNT, NSTEP/IPROG
  END IF

END DO ! ISTEP

CLOSE(10)


!========================================================
      

call begin_watch(20, nWatch, stop_watch)
call avr_Vars(nx, ny, nstep, UAVE, VAVE, flop)
call end_watch(20, nWatch, stop_watch)


call begin_watch(21, nWatch, stop_watch)
call wrt_Profile(nx, ny, UAVE, flop)
call end_watch(21, nWatch, stop_watch)


!call begin_watch(22, nWatch, stop_watch)
!call wrt_RCSAvr(nx, ny, nz, time, nstep, UAVE, VAVE, A, HUAVE, HVAVE, HWAVE, flop)
!call end_watch(22, nWatch, stop_watch)


call begin_watch(23, nWatch, stop_watch)
call wrt_SPHAvr(nx, ny, nz, dx, dy, nstep, time, UAVE, VAVE, HUVW, flop)
call end_watch(23, nWatch, stop_watch)


!========================================================


call print_watch(nWatch, stop_watch, tm_label)

END program main


! ***************************************************
  subroutine gen_grid(nx, ny, xlen, ylen, X, Y, flop)

  implicit none
  double precision, dimension(NX,NY)         :: X, Y
  integer          :: i, j, nx, ny
  double precision :: XLEN, YLEN
  double precision :: flop

  DO J=1,NY
  DO I=1,NX
    X(I,J)=XLEN/DFLOAT(NX-1)*DFLOAT(I-1)-5.D0
    Y(I,J)=YLEN/DFLOAT(NY-1)*DFLOAT(J-1)-5.D0
  END DO
  END DO

  flop = NX*NY*6.0

  return
  end subroutine gen_grid


! ***********************************************************************
  subroutine PDM_profile(nx, ny, IA, IB, JA, JB, CDD, PFUNC, DFUNC, flop)

  implicit none
  double precision, dimension(NX,NY)         :: DFUNC
  double precision, dimension(NY)            :: PFUNC
  integer          :: i, j, nx, ny, IA, IB, JA, JB
  double precision :: flop, cdd

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

  flop = 0.0

  return
  end subroutine PDM_profile


! **************************************************
  subroutine Inflow_profile(nstep, UUIN, VVIN, flop)

  implicit none
  double precision, dimension(NSTEP)         :: UUIN, VVIN
  integer          :: i, nstep
  double precision :: flop

  OPEN(1,FILE='time-hist-inflow_5.0deg.dat',FORM='FORMATTED')
  READ(1,'(2F20.15)') (UUIN(i),VVIN(i),i=1,NSTEP)
  CLOSE(1)


!   DO I=1,NSTEP
!     UUIN(I)=1.D0
!     VVIN(I)=0.D0
!   END DO

  flop = 0.0

  return
  end subroutine Inflow_profile


! **************************************************************************************************
  subroutine init_Arrays(nx, ny, nz, nstep, UU, VV, PP, UAVE, VAVE, UUIN, VVIN, HW, HWAVE, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, PP
  double precision, dimension(NX,NY)         :: UAVE, VAVE
  double precision, dimension(NSTEP)         :: UUIN, VVIN
  real,             dimension(NX,NY,NZ)      :: HW, HWAVE
  integer          :: i, j, k, nx, ny, nz, nstep
  double precision :: flop


  DO J=1,NY
  DO I=1,NX
    UU(I,J)=UUIN(1)
    VV(I,J)=VVIN(1)
    PP(I,J)=0.D0
    UAVE(I,J)=0.D0
    VAVE(I,J)=0.D0
!   PAVE(I,J)=0.D0
  END DO
  END DO

  DO K=1,NZ
  DO J=1,NY
  DO I=1,NX
    HW(I,J,K)=0.0
    HWAVE(I,J,K)=0.0
  END DO
  END DO
  END DO

  flop = 0.0

  return
  end subroutine init_Arrays


! ********************************************************
  subroutine read_rst_Data(nx, ny, time, UU, VV, PP, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, PP
  integer          :: i, j, lx, ly, nx, ny
  double precision :: dre, time
  double precision :: flop

  OPEN(1,FILE='uvps.dat')
  READ(1,*) TIME, dre, lx, ly
  READ(1,*) ((UU(I,J),I=0,NX+2),J=0,NY+2)
  READ(1,*) ((VV(I,J),I=0,NX+2),J=0,NY+2)
  READ(1,*) ((PP(I,J),I=0,NX+2),J=0,NY+2)
  CLOSE(1)

  flop = 0.0

  return
  end subroutine read_rst_Data


! **************************************************
  subroutine store_vel(nx, ny, UO, VO, UU, VV, flop)

  implicit none
  double precision, intent(in),  dimension(0:NX+2,0:NY+2) :: UU, VV
  double precision, intent(out), dimension(0:NX+2,0:NY+2) :: UO, VO
  integer          :: i, j, nx, ny
  double precision :: flop

  DO J=0,NY+2
  DO I=0,NX+2
    UO(I,J)=UU(I,J)
    VO(I,J)=VV(I,J)
  END DO
  END DO

  flop = 0.0

  return
  end subroutine store_vel


! **************************************************
  subroutine store_vel2(nx, ny, UO, VO, UU, VV, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, UO, VO
  integer          :: i, j, nx, ny
  double precision :: flop

  DO J=0,NY+2
  DO I=0,NX+2
    UO(I,J)=UU(I,J)
  END DO
  END DO

  DO J=0,NY+2
  DO I=0,NX+2
    VO(I,J)=VV(I,J)
  END DO
  END DO

  flop = 0.0

  return
  end subroutine store_vel2


! **************************************************
  subroutine store_vel3(nx, ny, UO, VO, UU, VV, flop)

  implicit none
  double precision, dimension((NX+3)*(NY+3)) :: UU, VV, UO, VO
  integer          :: i, nx, ny
  double precision :: flop

  DO I=1,(NX+3)*(NY+3)
    UO(I)=UU(I)
    VO(I)=VV(I)
  END DO

  flop = 0.0

  return
  end subroutine store_vel3


! *********************************************************************
  subroutine uflux(nx, ny, dt, rei, dx, dy, UO, VO, UU, DFUNC, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: UU, UO, VO
  double precision, dimension(NX,NY)         :: DFUNC
  integer          :: i, j, nx, ny
  double precision :: CU1, CU2, CV1, CV2
  double precision :: AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL
  double precision :: UU1, UU2, UU3, UU4, U2, U4, ADVX, VISX, SDRAGX
  double precision :: dt, rei, DX, DY, DXI, DYI, DDXI, DDYI
  double precision :: flop

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  DDXI= DXI/DX
  DDYI= DYI/DY

  DO J=2,NY
  DO I=2,NX-1

!   CONVECTIVE TERM(4TH-ORDER CENTRAL SCHEME)
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

!   4TH-ORDER NUMERICAL VISCOSITY(K-K SCHEME TYPE)
    UIJK=UO(I,J)*DXI
    VIJK=.25D0*(VO(I+1,J-1)+VO(I,J-1)+VO(I+1,J)+VO(I,J))*DYI
    U2=UO(I+2,J)-4.D0*(UO(I+1,J)+UO(I-1,J))+UO(I-2,J)+6.D0*UO(I,J)
    U4=UO(I,J+2)-4.D0*(UO(I,J+1)+UO(I,J-1))+UO(I,J-2)+6.D0*UO(I,J)
    ADVX=.5D0*(AX1+AX2+AY1+AY2)+(DABS(UIJK)*U2+DABS(VIJK)*U4)/4.D0

!   VISCOUS TERM(2ND-ORDER CENTRAL SCHEME)
    VISX=((UO(I+1,J)-2.D0*UO(I,J)+UO(I-1,J))*DDXI &
         +(UO(I,J+1)-2.D0*UO(I,J)+UO(I,J-1))*DDYI)*REI

!   POROUS DISK MODEL
    DVEL=DSQRT(UO(I,J)**2+VO(I,J)**2)
    SDRAGX=-UO(I,J)*DVEL*DFUNC(I,J)

!   INTERMEDIATE VELOCITY
    UU(I,J)=UO(I,J)+DT*(-ADVX+VISX+SDRAGX)
  END DO
  END DO

  flop = (NY-1)*(NX-2)*107.0

  return
  end subroutine uflux


! *********************************************************************
  subroutine vflux(nx, ny, dt, rei, dx, dy, UO, VO, VV, DFUNC, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: VV, UO, VO
  double precision, dimension(NX,NY)         :: DFUNC
  integer          :: i, j, nx, ny
  double precision :: CU1, CU2, CV1, CV2
  double precision :: AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL
  double precision :: VV1, VV2, VV3, VV4, V2, V4, ADVY, VISY, SDRAGY
  double precision :: dt, rei, DX, DY, DXI, DYI, DDXI, DDYI
  double precision :: flop

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  DDXI= DXI/DX
  DDYI= DYI/DY

  DO J=2,NY-1
  DO I=2,NX

!   CONVECTIVE TERM(4TH-ORDER CENTRAL SCHEME)
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

!   4TH-ORDER NUMERICAL VISCOSITY(K-K SCHEME TYPE)
    UIJK=.25D0*(UO(I-1,J+1)+UO(I-1,J)+UO(I,J+1)+UO(I,J))*DXI
    VIJK=VO(I,J)*DYI
    V2=VO(I+2,J)-4.D0*(VO(I+1,J)+VO(I-1,J))+VO(I-2,J)+6.D0*VO(I,J)
    V4=VO(I,J+2)-4.D0*(VO(I,J+1)+VO(I,J-1))+VO(I,J-2)+6.D0*VO(I,J)
    ADVY=.5D0*(AX1+AX2+AY1+AY2)+(DABS(UIJK)*V2+DABS(VIJK)*V4)/4.D0

!   VISCOUS TERM(2ND-ORDER CENTRAL SCHEME)
    VISY=((VO(I+1,J)-2.D0*VO(I,J)+VO(I-1,J))*DDXI &
         +(VO(I,J+1)-2.D0*VO(I,J)+VO(I,J-1))*DDYI)*REI

!   POROUS DISK MODEL
    DVEL=DSQRT(UO(I,J)**2+VO(I,J)**2)
    SDRAGY=-VO(I,J)*DVEL*DFUNC(I,J)

!   INTERMEDIATE VELOCITY
    VV(I,J)=VO(I,J)+DT*(-ADVY+VISY+SDRAGY)
  END DO
  END DO

  flop = (NY-2)*(NX-1)*107.0

  return
  end subroutine vflux


! *************************************************************
  subroutine Poisson_RHS(nx, ny, dt, dx, dy, RHS, UU, VV, flop)

  implicit none
  double precision, dimension(NX,NY)         :: RHS
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV
  integer          :: i, j, nx, ny
  double precision :: dt, DX, DY, DXI, DYI, ddt
  double precision :: flop

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  ddt = 1.D0/dt

  DO J=2,NY
  DO I=2,NX
    RHS(I,J)=-((UU(I,J)-UU(I-1,J))*DXI &
              +(VV(I,J)-VV(I,J-1))*DYI)*ddt
  END DO
  END DO

  flop = (NY-1)*(NX-1)*7.0

  return
  end subroutine Poisson_RHS


! *********************************************************************
  subroutine Poisson_AXB(nx, ny, dx, dy, omega, err, RHS, PP, flop)

  implicit none
  double precision, dimension(NX,NY)         :: RHS
  double precision, dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j, nx, ny
  double precision :: DX, DY, DXI, DYI, DDXI, DDYI
  double precision :: PAGE, dp, err, omega
  double precision :: flop

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  DDXI= DXI/DX
  DDYI= DYI/DY

  DO J=2,NY
  DO I=2,NX
    PAGE=(PP(I+1,J)+PP(I-1,J))*DDXI &
        +(PP(I,J+1)+PP(I,J-1))*DDYI+RHS(I,J)
    DP=.5D0*PAGE/(DDXI+DDYI)-PP(I,J)
    ERR=ERR+DP*DP
    PP(I,J)=PP(I,J)+OMEGA*DP
  END DO
  END DO

  flop = (NY-1)*(NX-1)*14.0

  return
  end subroutine Poisson_AXB


! *********************************************************************
  subroutine Poisson_AXB2C(nx, ny, dx, dy, omega, err, RHS, PP, flop)

  implicit none
  double precision, dimension(NX,NY)         :: RHS
  double precision, dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j, nx, ny, color
  double precision :: DX, DY, DXI, DYI, DDXI, DDYI
  double precision :: PAGE, dp, err, omega
  double precision :: flop

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  DDXI= DXI/DX
  DDYI= DYI/DY

  do color=0,1

    DO J=2,NY
    do i=2+mod(j+color,2), nx, 2
      PAGE=(PP(I+1,J)+PP(I-1,J))*DDXI &
          +(PP(I,J+1)+PP(I,J-1))*DDYI+RHS(I,J)
      DP=.5D0*PAGE/(DDXI+DDYI)-PP(I,J)
      ERR=ERR+DP*DP
      PP(I,J)=PP(I,J)+OMEGA*DP
    END DO
    END DO

  end do

  flop = (NY-1)*(NX-1)*14.0

  return
  end subroutine Poisson_AXB2C


! **************************************
  subroutine Shift_Prs(nx, ny, PP, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j, nx, ny
  double precision :: PPP
  double precision :: flop

  PPP=0.D0

  DO J=1,NY+1
  DO I=1,NX+1
    PPP=PPP+PP(I,J)
  END DO
  END DO

  PPP = PPP/DFLOAT((NX+1)*(NY+1))
  
  DO J=1,NY+1
  DO I=1,NX+1
    PP(I,J)=PP(I,J)-PPP
  END DO
  END DO

  flop = (NY+1)*(NX+1)*4.0

  return
  end subroutine Shift_Prs


! ***********************************
  subroutine BC_Prs(nx, ny, PP, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j, nx, ny
  double precision :: flop

  DO I=1,NX+1
    PP(I,1)   =PP(I,2)
    PP(I,NY+1)=PP(I,NY)
  END DO

  DO J=1,NY+1
    PP(1,J)   =PP(2,J)
    PP(NX+1,J)=PP(NX,J)
  END DO

  flop = 0.0

  return
  end subroutine BC_Prs


! ********************************************************
  subroutine Prj_Vel(nx, ny, dt, dx, dy, UU, VV, PP, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, PP
  integer          :: i, j, nx, ny
  double precision :: DX, DY, DXI, DYI, dt
  double precision :: flop

  DXI = 1.D0/DX
  DYI = 1.D0/DY

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

  flop = (NY-1)*(NX-2)*4.0+(NY-2)*(NX-1)*4.0

  return
  end subroutine Prj_Vel


! *********************************************************************************
  subroutine BC_Vel(nx, ny, dt, dx, dy, nstep, istep, UU, VV, UO, UUIN, VVIN, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, UO
  double precision, dimension(NSTEP)         :: UUIN, VVIN
  integer          :: i, j, istep, nx, ny, nstep
  double precision :: DX, DY, DXI, DYI, dt
  double precision :: flop

  DXI = 1.D0/DX
  DYI = 1.D0/DY

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

  flop = (NY-1)*8.0

  return
  end subroutine BC_Vel


! *********************************************
  subroutine stg2CC(nx, ny, U, V, UU, VV, flop)

  implicit none
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV
  double precision, dimension(NX,NY)         :: U, V
  integer          :: i, j, nx, ny
  double precision :: flop

  DO J=1,NY
  DO I=1,NX
    U(I,J)=.5D0*(UU(I,J)+UU(I,J+1))
    V(I,J)=.5D0*(VV(I,J)+VV(I+1,J))
  END DO
  END DO

  flop = NX*NY*4.0

  return
  end subroutine stg2CC


! ****************************************************
  subroutine accumVars(nx, ny, U, V, UAVE, VAVE, flop)

  implicit none
  double precision, dimension(NX,NY)         :: UAVE, VAVE
  double precision, dimension(NX,NY)         :: U, V
  integer          :: i, j, nx, ny
  double precision :: flop

  DO J=1,NY
  DO I=1,NX
    UAVE(I,J)=UAVE(I,J)+U(I,J)
    VAVE(I,J)=VAVE(I,J)+V(I,J)
  END DO
  END DO

  flop = NX*NY*2.0

  return
  end subroutine accumVars


! ************************************************************
  subroutine wrt_RCSIns(nx, ny, nz, istep, time, A, U, V, HU, HV, HW, flop)

  implicit none
  real,             dimension(NX,NY,NZ)      :: HU, HV, HW
  double precision, dimension(NX,NY)         :: U, V
  double precision, dimension(2)             :: A
  integer          :: i, j, k, nx, ny, nz, istep, ii
  double precision :: flop, time
  real             :: t

  DO J=1,NY
  DO I=1,NX
    t = SNGL(U(I,J))
    HU(I,J,1)=t
    HU(I,J,2)=t
    HU(I,J,3)=t
  END DO
  END DO

  DO J=1,NY
  DO I=1,NX
    t = SNGL(V(I,J))
    HV(I,J,1)=t
    HV(I,J,2)=t
    HV(I,J,3)=t
  END DO
  END DO

!   HW(I,J,K)=0.0    initialize at init_Arrays()

  WRITE(10) ISTEP, SNGL(TIME)
  WRITE(10) (SNGL(A(II)),II=1,2)
  WRITE(10) (((HU(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
  WRITE(10) (((HV(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
  WRITE(10) (((HW(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

  flop = 0.0

  return
  end subroutine wrt_RCSIns


! *****************************************************************
!  subroutine wrt_VisIns(nx, ny, nz, istep, time, A, HU, HV, HW, flop)

!  implicit none
!  real,             dimension(NX,NY,NZ)      :: HU, HV, HW
!  double precision, dimension(2)             :: A
!  integer          :: i, j, k, istep, ii, nx, ny, nz
!  double precision :: time
!  double precision :: flop

!  WRITE(10) ISTEP, SNGL(TIME)
!  WRITE(10) (SNGL(A(II)),II=1,2)
!  WRITE(10) (((HU(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
!  WRITE(10) (((HV(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
!  WRITE(10) (((HW(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

!  flop = 0.0

!  return
!  end subroutine wrt_VisIns


! **********************************************************************************
  subroutine wrt_SPHIns(nx, ny, nz, dx, dy, istep, time, U, V, P, HUVW, HP, flop)

  implicit none
  real,             dimension(3,NX,NY,NZ)    :: HUVW
  real,             dimension(NX,NY,NZ)      :: HP
  double precision, dimension(NX,NY)         :: U, V, P
  integer          :: i, j, k, istep, nx, ny, nz, svType, dType
  double precision :: time, dx, dy
  double precision :: flop
  real             :: xorg, yorg, zorg, t, u1, u2, u3
  real             :: xpch, ypch, zpch, tm
  character(30)    :: s

  svType = 1 ! scalar
  dType  = 1 ! float 
  xorg = -5.0
  yorg = -5.0
  zorg = 0.0
  xpch = SNGL(dx)
  ypch = SNGL(dy)
  zpch = SNGL(dx)
  tm   = SNGL(time)

  DO J=1,NY
  DO I=1,NX
    t = SNGL(P(I,J))
    HP(I,J,1)=t
    HP(I,J,2)=t
    HP(I,J,3)=t
  END DO
  END DO

  write (s, '("data/prs_",I9.9,".sph")') istep

  open (unit=22,file=s,form='unformatted')
  write (22) svType, dType
  write (22) nx, ny, nz
  write (22) xorg, yorg, zorg
  write (22) xpch, ypch, zpch
  write (22) istep, tm
  write (22) HP
  close (unit=22)

  u3 = 0.0
  DO J=1,NY
  DO I=1,NX
    u1 = SNGL(U(I,J))
    u2 = SNGL(V(I,J))
    HUVW(1,I,J,1)=u1
    HUVW(2,I,J,1)=u2
    HUVW(3,I,J,1)=u3
    HUVW(1,I,J,2)=u1
    HUVW(2,I,J,2)=u2
    HUVW(3,I,J,2)=u3
    HUVW(1,I,J,3)=u1
    HUVW(2,I,J,3)=u2
    HUVW(3,I,J,3)=u3
  END DO
  END DO

  svType = 2 ! vector
  write (s, '("data/uvw_",I9.9,".sph")') istep

  open (unit=22,file=s,form='unformatted')
  write (22) svType, dType
  write (22) nx, ny, nz
  write (22) xorg, yorg, zorg
  write (22) xpch, ypch, zpch
  write (22) istep, tm
  write (22) HUVW
  close (unit=22)

  flop = 0.0

  return
  end subroutine wrt_SPHIns


! ****************************************************
  subroutine avr_Vars(nx, ny, nstep, UAVE, VAVE, flop)

  implicit none
  double precision, dimension(NX,NY)         :: UAVE, VAVE
  integer          :: i, j, nx, ny, nstep
  double precision :: dn
  double precision :: flop

  dn = 1.D0/DFLOAT(NSTEP)

  DO J=1,NY
  DO I=1,NX
    UAVE(I,J)=UAVE(I,J) * dn
    VAVE(I,J)=VAVE(I,J) * dn
!   PAVE(I,J)=PAVE(I,J) * dn
  END DO
  END DO

  flop = NX*NY*2.0

  return
  end subroutine avr_Vars


! ******************************************
  subroutine wrt_Profile(nx, ny, UAVE, flop)

  implicit none
  double precision, dimension(NX,NY)         :: UAVE
  integer          :: j, nx, ny
  double precision :: flop

  OPEN(1,FILE='5D-10D-U-profile.dat',FORM='FORMATTED')
  WRITE(1,'(2F15.7)')(UAVE(1001,j),UAVE(1501,j),j=1,NY)
  CLOSE(1)

  flop = 0.0

  return
  end subroutine wrt_Profile


! *********************************************************************
  subroutine wrt_RCSAvr(nx, ny, nz, time, nstep, UAVE, VAVE, A, HUAVE, HVAVE, HWAVE, flop)

  implicit none
  real,             dimension(NX,NY,NZ)      :: HUAVE, HVAVE, HWAVE
  double precision, dimension(NX,NY)         :: UAVE, VAVE
  double precision, dimension(2)             :: A
  integer          :: i, j, k, ii, nx, ny, nz, nstep
  double precision :: time
  double precision :: flop
  real             :: t


  DO J=1,NY
  DO I=1,NX
    t = SNGL(UAVE(I,J))
    HUAVE(I,J,1)=t
    HUAVE(I,J,2)=t
    HUAVE(I,J,3)=t
  END DO
  END DO

  DO J=1,NY
  DO I=1,NX
    t = SNGL(VAVE(I,J))
    HVAVE(I,J,1)=VAVE(I,J)
    HVAVE(I,J,2)=VAVE(I,J)
    HVAVE(I,J,3)=VAVE(I,J)
  END DO
  END DO

!   HWAVE(I,J,K)=0.0       initialize at init_Arrays()

  open(unit=11,file='3d-tave-vis.dat',form='unformatted')
  write(11) NSTEP,SNGL(TIME)
  write(11) (SNGL(A(II)),II=1,2)
  write(11) (((HUAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  write(11) (((HVAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  write(11) (((HWAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  close(11)

  flop = 0.0

  return
  end subroutine wrt_RCSAvr


! **********************************************************************************
  subroutine wrt_SPHAvr(nx, ny, nz, dx, dy, istep, time, UAVE, VAVE, HUVW, flop)

  implicit none
  real,             dimension(3,NX,NY,NZ)    :: HUVW
  double precision, dimension(NX,NY)         :: UAVE, VAVE
  double precision, dimension(NX,NY)         :: U, V
  integer          :: i, j, k, istep, nx, ny, nz, svType, dType
  double precision :: time, dx, dy
  double precision :: flop
  real             :: xorg, yorg, zorg, u1, u2, u3
  real             :: xpch, ypch, zpch, tm
  character(30)    :: s

  svType = 2 ! vector
  dType  = 1 ! float 
  xorg = -5.0
  yorg = -5.0
  zorg = 0.0
  xpch = SNGL(dx)
  ypch = SNGL(dy)
  zpch = SNGL(dx)
  tm   = SNGL(time)

  u3 = 0.0
  DO J=1,NY
  DO I=1,NX
    u1 = SNGL(UAVE(I,J))
    u2 = SNGL(VAVE(I,J))
    HUVW(1,I,J,1)=u1
    HUVW(2,I,J,1)=u2
    HUVW(3,I,J,1)=u3
    HUVW(1,I,J,2)=u1
    HUVW(2,I,J,2)=u2
    HUVW(3,I,J,2)=u3
    HUVW(1,I,J,3)=u1
    HUVW(2,I,J,3)=u2
    HUVW(3,I,J,3)=u3
  END DO
  END DO

  write (s, '("data/uvwa_",I9.9,".sph")') istep
  
  open (unit=22,file=s,form='unformatted')
  write (22) svType, dType
  write (22) nx, ny, nz
  write (22) xorg, yorg, zorg
  write (22) xpch, ypch, zpch
  write (22) istep, tm
  write (22) HUVW
  close (unit=22)

  flop = 0.0

  return
  end subroutine wrt_SPHAvr


! *********************************************************************
! sw(1,*) temporary working
! sw(2,*) accumlated time
! sw(3,*) accumlated count
! sw(4,*) flop/call

subroutine init_watch(nx, ny, N, sw, label)
implicit none
double precision, dimension(4,N)      :: sw
character(20)                         :: label(N)
integer          :: N, i, nx, ny

do i=1,N
  sw(1,i) = 0.0
  sw(2,i) = 0.0
  sw(3,i) = 0.0
  sw(4,i) = 0.0
end do

label(1) = 'meshing'
label(2) = "setup_PDM"
label(3) = "read_inflow"
label(4) = "init_arrays"
label(5) = "read_rst_file"
label(6) = "str_Nstep"
label(7) = "Uflux"
label(8) = "Vflux"
label(9) = "Poi_RHS"
label(10) = "Poi_AXB"
label(11) = "Shift_Prs"
label(12) = "BC_Prs"
label(13) = "Update_Vec"
label(14) = "BC_Vec"
label(15) = "stg2CC"
label(16) = "accumAvr"
label(17) = "display"
label(18) = "wrt_RCSIns"
label(19) = "wrt_SPHIns"
label(20) = "averaging"
label(21) = "wrt_Profile"
label(22) = "wrt_RCSAvr"
label(23) = "wrt_SPHAvr"

sw(4,1) = 6.0*NX*NY
sw(4,2) = 0.0
sw(4,3) = 0.0
sw(4,4) = 0.0
sw(4,5) = 0.0
sw(4,6) = 0.0
sw(4,7) = (NY-1)*(NX-2)*107.0
sw(4,8) = (NY-2)*(NX-1)*107.0
sw(4,9) = (NY-1)*(NX-1)*7.0
sw(4,10) = (NY-1)*(NX-1)*14.0
sw(4,11) = (NY+1)*(NX+1)*4.0
sw(4,12) = 0.0
sw(4,13) = (NY-1)*(NX-2)*4.0+(NY-2)*(NX-1)*4.0
sw(4,14) = (NY-1)*8.0
sw(4,15) = NX*NY*4.0
sw(4,16) = NX*NY*2.0
sw(4,17) = 0.0
sw(4,18) = 0.0
sw(4,19) = 0.0
sw(4,20) = NX*NY*2.0
sw(4,21) = 0.0
sw(4,22) = 0.0
sw(4,23) = 0.0

end subroutine init_watch


subroutine begin_watch(key, N, sw)
!$ use omp_lib
implicit none
double precision, dimension(4,N)      :: sw
integer          :: key, N

!$ sw(1,key) = omp_get_wtime()

end subroutine begin_watch


subroutine end_watch(key, N, sw)
!$ use omp_lib
implicit none
double precision, dimension(4,N)      :: sw
integer          :: key, N
double precision :: tm_ed, tm_st

!$ tm_ed = omp_get_wtime()
tm_st = sw(1,key)
sw(2,key) = sw(2,key) + (tm_ed - tm_st)
sw(3,key) = sw(3,key) + 1.0

end subroutine end_watch


subroutine print_watch(N, sw, label)
implicit none
double precision, dimension(4,N)          :: sw
character(20)                             :: label(N)
integer,allocatable,dimension(:)          :: idx
double precision,allocatable,dimension(:) :: val
integer          :: N, i, j
double precision :: sum

allocate( idx(N) )
allocate( val(N) )

! to avoid obtaining false Gflops
do i=1,N
  if (sw(2,i)<1.0e-3) sw(2,i)=1.0e-3
end do

! copy sort target and set index
do i=1,N
 idx(i)=i
 val(i)=sw(2,i)
end do

sum =0.0
do i=1,N
  sum = sum + val(i)
end do

! obtain sorted index
do i=1,N-1
  do j=i+1,N
    if (val(i)<val(j)) call swap_watch(i,j,N,idx,val)
  end do
end do


print "(A21,A7,A10,A15,A15,A12,A10)", 'Label', '%Time', 'Call', 'accTime', 'avrTime', 'flop/call', 'Gflops'


do j=1,N
  i = idx(j)
  print "(A20,F8.2,I10,2E15.5,E12.3,F10.3)", &
       label(i), &
       sw(2,i)/sum*100.0, &
       int(sw(3,i)), &
       sw(2,i), &
       sw(2,i)/sw(3,i), &
       sw(4,i), &
       sw(4,i)*sw(3,i)/sw(2,i)*1.0E-9
end do

end subroutine print_watch


subroutine swap_watch(i,j,N,idx,val)
implicit none
double precision, dimension(N)      :: val
integer,          dimension(N)      :: idx
integer          :: i,j,N,m
double precision :: a

a = val(i)
m = idx(i)

val(i) = val(j)
idx(i) = idx(j)

val(j) = a
idx(j) = m

end subroutine swap_watch
