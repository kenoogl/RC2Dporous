program main
!$ use omp_lib

implicit none
include 'size.fi'
integer, PARAMETER :: NSTEP=320000, nWatch=23

double precision, dimension(NX,NY)         :: RHS, U, V, P, DFUNC
!double precision, dimension(NX,NY)         :: X, Y
double precision, dimension(0:NX+2,0:NY+2) :: PP, UU, VV, UO, VO, UAVE, VAVE
real,             dimension(NX,NY,NZ)      :: HUAVE, HVAVE, HWAVE, HU, HV, HW, HP
real,             dimension(3,NX,NY,NZ)    :: HUVW
double precision, dimension(NSTEP)         :: UUIN, VVIN
double precision, dimension(NY)            :: PFUNC
double precision, dimension(2)             :: A

double precision, dimension(4,nWatch)      :: stop_watch
character(20)                              :: tm_label(nWatch)

double precision :: omega, eps, time, dt, RE, REI, ERR, CDD
double precision :: DX, DY, XLEN, YLEN
integer          :: ITER, ICON, IPROG, NSOR, ICOUNT
integer          :: NOXY, ISTEP, LOOP
integer          :: IA, IB, JA, JB

!========================================================


call init_watch(nWatch, stop_watch, tm_label)


! Parameters

! SOR
OMEGA=.9D0
EPS  =1.E-3
ITER =100

! TIME STEP ETC
ICON  = 0      ! initial start=0, restart=1
IPROG = 400000    ! Interval for vis. data output
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


!call begin_watch(1, nWatch, stop_watch)
!call gen_grid(xlen, ylen, X, Y)
!call end_watch(1, nWatch, stop_watch)


call begin_watch(2, nWatch, stop_watch)
call PDM_profile(IA, IB, JA, JB, CDD, PFUNC, DFUNC)
call end_watch(2, nWatch, stop_watch)


call begin_watch(3, nWatch, stop_watch)
call Inflow_profile(NSTEP, UUIN, VVIN)
call end_watch(3, nWatch, stop_watch)


call begin_watch(4, nWatch, stop_watch)
call init_Arrays(nstep, UU, VV, PP, UAVE, VAVE, UUIN, VVIN, HW, HWAVE)
call end_watch(4, nWatch, stop_watch)


! Restart
IF (ICON==1)THEN
  call begin_watch(5, nWatch, stop_watch)
  call read_rst_Data(time, UU, VV, PP)
  call end_watch(5, nWatch, stop_watch)
END IF


!========================================================


OPEN(10,FILE='3d-inst-vis.dat',FORM='UNFORMATTED')
ICOUNT = 0

DO ISTEP=1, NSTEP
  TIME=TIME+DT

  call begin_watch(6, nWatch, stop_watch)
  !call store_vel(UO, VO, UU, VV)
  call store_vel3(UO, VO, UU, VV)
  call end_watch(6, nWatch, stop_watch)


  call begin_watch(7, nWatch, stop_watch)
  call uflux(dt, rei, dx, dy, UO, VO, UU, DFUNC)
  call end_watch(7, nWatch, stop_watch)


  call begin_watch(8, nWatch, stop_watch)
  call vflux(dt, rei, dx, dy, UO, VO, VV, DFUNC)
  call end_watch(8, nWatch, stop_watch)


  call begin_watch(9, nWatch, stop_watch)
  call Poisson_RHS(dt, dx, dy, RHS, UU, VV)
  call end_watch(9, nWatch, stop_watch)

!--------------------------------------------------------


  DO LOOP=1,ITER
    ERR=0.D0

    call begin_watch(10, nWatch, stop_watch)
    !call Poisson_AXB(dx, dy, omega, err, RHS, PP)
    call Poisson_AXB2C(dx, dy, omega, err, RHS, PP)
    call end_watch(10, nWatch, stop_watch)

    ERR=DSQRT(ERR/DFLOAT(NOXY))
    IF (ERR < EPS) exit
  END DO
  

  call begin_watch(11, nWatch, stop_watch)
  call Shift_Prs(PP)
  call end_watch(11, nWatch, stop_watch)


  call begin_watch(12, nWatch, stop_watch)
  call BC_Prs(PP)
  call end_watch(12, nWatch, stop_watch)

  call begin_watch(13, nWatch, stop_watch)
  call Prj_Vel(dt, dx, dy, UU, VV, PP)
  call end_watch(13, nWatch, stop_watch)


  call begin_watch(14, nWatch, stop_watch)
  call BC_Vel(dt, dx, dy, nstep, istep, UU, VV, UO, UUIN, VVIN)
  call end_watch(14, nWatch, stop_watch)
  

!--------------------------------------------------------


  call begin_watch(16, nWatch, stop_watch)
  call accumVars(UU, VV, UAVE, VAVE)
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

    call begin_watch(15, nWatch, stop_watch)
    call stg2Node(U, V, UAVE, VAVE)
    call end_watch(15, nWatch, stop_watch)

    !call begin_watch(18, nWatch, stop_watch)
    !call wrt_RCSIns(istep, time, A, U, V, HU, HV, HW)
    !call end_watch(18, nWatch, stop_watch)

    call begin_watch(19, nWatch, stop_watch)
    call wrt_SPHIns(dx, dy, istep, time, U, V, PP, HUVW, HP)
    call end_watch(19, nWatch, stop_watch)

    ICOUNT=ICOUNT+1
    print "('RC-Scope date=',I10, '/', I10,'(Total)')", ICOUNT, NSTEP/IPROG
  END IF

END DO ! ISTEP

CLOSE(10)


!========================================================
      

call begin_watch(20, nWatch, stop_watch)
call avr_Vars(nstep, UAVE, VAVE)
call end_watch(20, nWatch, stop_watch)


call begin_watch(21, nWatch, stop_watch)
call wrt_Profile(UAVE)
call end_watch(21, nWatch, stop_watch)


!call begin_watch(22, nWatch, stop_watch)
!call wrt_RCSAvr(time, nstep, UAVE, VAVE, A, HUAVE, HVAVE, HWAVE)
!call end_watch(22, nWatch, stop_watch)


call begin_watch(23, nWatch, stop_watch)
call wrt_SPHAvr(dx, dy, nstep, time, UAVE, VAVE, HUVW)
call end_watch(23, nWatch, stop_watch)


!========================================================


call print_watch(nWatch, stop_watch, tm_label)

END program main


! *********************************************
  subroutine gen_grid(xlen, ylen, X, Y)

  implicit none
  include 'size.fi'
  double precision, intent(out), dimension(NX,NY)         :: X, Y
  integer          :: i, j
  double precision :: XLEN, YLEN

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP FIRSTPRIVATE(XLEN, YLEN) &
!$OMP PRIVATE(i,j)
  DO J=1,NY
  DO I=1,NX
    X(I,J)=XLEN/DFLOAT(NX-1)*DFLOAT(I-1)-5.D0
    Y(I,J)=YLEN/DFLOAT(NY-1)*DFLOAT(J-1)-5.D0
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine gen_grid


! ***********************************************************************
  subroutine PDM_profile(IA, IB, JA, JB, CDD, PFUNC, DFUNC)

  implicit none
  include 'size.fi'
  double precision, dimension(NX,NY)         :: DFUNC
  double precision, dimension(NY)            :: PFUNC
  integer          :: i, j, IA, IB, JA, JB
  double precision :: cdd

  OPEN(1,FILE='../pfunc.dat')
  READ(1,'(F15.7)')(PFUNC(J),J=JA,JB)
  CLOSE(1)

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i,j)
  DO J=1,NY
  DO I=1,NX
    DFUNC(I,J)=0.D0
  END DO
  END DO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP FIRSTPRIVATE(CDD, IA, IB, JA, JB) &
!$OMP PRIVATE(i,j)
  DO J=JA,JB
  DO I=IA,IB
    DFUNC(I,J)=PFUNC(J)*CDD
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine PDM_profile


! ********************************************
  subroutine Inflow_profile(nstep, UUIN, VVIN)

  implicit none
  double precision, dimension(NSTEP)         :: UUIN, VVIN
  integer          :: i, nstep

  OPEN(1,FILE='../time-hist-inflow_5.0deg.dat',FORM='FORMATTED')
  READ(1,'(2F20.15)') (UUIN(i),VVIN(i),i=1,NSTEP)
  CLOSE(1)


!   DO I=1,NSTEP
!     UUIN(I)=1.D0
!     VVIN(I)=0.D0
!   END DO

  return
  end subroutine Inflow_profile


! **************************************************************************************************
  subroutine init_Arrays(nstep, UU, VV, PP, UAVE, VAVE, UUIN, VVIN, HW, HWAVE)

  implicit none
  include 'size.fi'
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, PP
  double precision, dimension(0:NX+2,0:NY+2) :: UAVE, VAVE
  double precision, dimension(NSTEP)         :: UUIN, VVIN
  real,             dimension(NX,NY,NZ)      :: HW, HWAVE
  integer          :: i, j, k, nstep

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i,j)
  DO J=1,NY
  DO I=1,NX
    UU(I,J)=UUIN(1)
    VV(I,J)=VVIN(1)
    PP(I,J)=0.D0
  END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i,j)
  DO J=0,NY+2
  DO I=0,NX+2
    UAVE(I,J)=0.D0
    VAVE(I,J)=0.D0
!   PAVE(I,J)=0.D0
  END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i,j,k)
  DO K=1,NZ
  DO J=1,NY
  DO I=1,NX
    HW(I,J,K)=0.0
    HWAVE(I,J,K)=0.0
  END DO
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine init_Arrays


! **************************************************
  subroutine read_rst_Data(time, UU, VV, PP)

  implicit none
  include 'size.fi'
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, PP
  integer          :: i, j, lx, ly
  double precision :: dre, time

  OPEN(1,FILE='uvps.dat')
  READ(1,*) TIME, dre, lx, ly
  READ(1,*) ((UU(I,J),I=0,NX+2),J=0,NY+2)
  READ(1,*) ((VV(I,J),I=0,NX+2),J=0,NY+2)
  READ(1,*) ((PP(I,J),I=0,NX+2),J=0,NY+2)
  CLOSE(1)

  return
  end subroutine read_rst_Data


! ********************************************
  subroutine store_vel(UO, VO, UU, VV)

  implicit none
  include 'size.fi'
  double precision, intent(in),  dimension(0:NX+2,0:NY+2) :: UU, VV
  double precision, intent(out), dimension(0:NX+2,0:NY+2) :: UO, VO
  integer          :: i, j

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i,j)
  DO J=0,NY+2
  DO I=0,NX+2
    UO(I,J)=UU(I,J)
    VO(I,J)=VV(I,J)
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine store_vel


! *********************************************
  subroutine store_vel3(UO, VO, UU, VV)

  implicit none
  include 'size.fi'
  double precision, intent(in),  dimension((NX+3)*(NY+3)) :: UU, VV
  double precision, intent(out), dimension((NX+3)*(NY+3)) :: UO, VO
  integer          :: i

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i)
  DO I=1,(NX+3)*(NY+3)
    UO(I)=UU(I)
    VO(I)=VV(I)
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine store_vel3


! ****************************************************
  subroutine uflux(dt, rei, dx, dy, UO, VO, UU, DFUNC)

  implicit none
  include 'size.fi'
  double precision, intent(out), dimension(0:NX+2,0:NY+2) :: UU
  double precision, intent(in),  dimension(0:NX+2,0:NY+2) :: UO, VO
  double precision, intent(in),  dimension(NX,NY)         :: DFUNC
  integer          :: i, j
  double precision :: CU1, CU2, CV1, CV2
  double precision :: AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL
  double precision :: UU1, UU2, UU3, UU4, U2, U4, ADVX, VISX, SDRAGX
  double precision :: dt, rei, DX, DY, DXI, DYI, DDXI, DDYI
  double precision :: C1, C2

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  DDXI= DXI/DX
  DDYI= DYI/DY
  C1 = 1.D0/16.D0
  C2 = 1.D0/24.D0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(dt, rei, DXI, DYI, DDXI, DDYI, C1, C2) &
!$OMP PRIVATE(CU1, CU2, CV1, CV2, i, j) &
!$OMP PRIVATE(AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL) &
!$OMP PRIVATE(UU1, UU2, UU3, UU4, U2, U4, ADVX, VISX, SDRAGX)
!$OMP DO SCHEDULE(static)
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
    AX1=(CU1*C1)*(UU1*C2*DXI)
    AX2=(CU2*C1)*(UU2*C2*DXI)
    AY1=(CV1*C1)*(UU3*C2*DYI)
    AY2=(CV2*C1)*(UU4*C2*DYI)

!   4TH-ORDER NUMERICAL VISCOSITY(K-K SCHEME TYPE)
    UIJK=UO(I,J)*DXI
    VIJK=.25D0*(VO(I+1,J-1)+VO(I,J-1)+VO(I+1,J)+VO(I,J))*DYI
    U2=UO(I+2,J)-4.D0*(UO(I+1,J)+UO(I-1,J))+UO(I-2,J)+6.D0*UO(I,J)
    U4=UO(I,J+2)-4.D0*(UO(I,J+1)+UO(I,J-1))+UO(I,J-2)+6.D0*UO(I,J)
    ADVX=.5D0*(AX1+AX2+AY1+AY2)+(DABS(UIJK)*U2+DABS(VIJK)*U4)*0.25D0

!   VISCOUS TERM(2ND-ORDER CENTRAL SCHEME)
    VISX=((UO(I+1,J)-2.D0*UO(I,J)+UO(I-1,J))*DDXI &
         +(UO(I,J+1)-2.D0*UO(I,J)+UO(I,J-1))*DDYI)*REI

!   POROUS DISK MODEL
    DVEL=DSQRT(UO(I,J)*UO(I,J)+VO(I,J)*VO(I,J))
    SDRAGX=-UO(I,J)*DVEL*DFUNC(I,J)

!   INTERMEDIATE VELOCITY
    UU(I,J)=UO(I,J)+DT*(-ADVX+VISX+SDRAGX)
  END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine uflux


! ****************************************************
  subroutine vflux(dt, rei, dx, dy, UO, VO, VV, DFUNC)

  implicit none
  include 'size.fi'
  double precision, intent(out), dimension(0:NX+2,0:NY+2) :: VV
  double precision, intent(in),  dimension(0:NX+2,0:NY+2) :: UO, VO
  double precision, intent(in),  dimension(NX,NY)         :: DFUNC
  integer          :: i, j
  double precision :: CU1, CU2, CV1, CV2
  double precision :: AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL
  double precision :: VV1, VV2, VV3, VV4, V2, V4, ADVY, VISY, SDRAGY
  double precision :: dt, rei, DX, DY, DXI, DYI, DDXI, DDYI
  double precision :: C1, C2

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  DDXI= DXI/DX
  DDYI= DYI/DY
  C1 = 1.D0/16.D0
  C2 = 1.D0/24.D0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(dt, rei, DXI, DYI, DDXI, DDYI, C1, C2) &
!$OMP PRIVATE(CU1, CU2, CV1, CV2, i, j) &
!$OMP PRIVATE(AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL) &
!$OMP PRIVATE(VV1, VV2, VV3, VV4, V2, V4, ADVY, VISY, SDRAGY)
!$OMP DO SCHEDULE(static)
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
    AX1=(CU1*C1)*(VV1*C2*DXI)
    AX2=(CU2*C1)*(VV2*C2*DXI)
    AY1=(CV1*C1)*(VV3*C2*DYI)
    AY2=(CV2*C1)*(VV4*C2*DYI)

!   4TH-ORDER NUMERICAL VISCOSITY(K-K SCHEME TYPE)
    UIJK=.25D0*(UO(I-1,J+1)+UO(I-1,J)+UO(I,J+1)+UO(I,J))*DXI
    VIJK=VO(I,J)*DYI
    V2=VO(I+2,J)-4.D0*(VO(I+1,J)+VO(I-1,J))+VO(I-2,J)+6.D0*VO(I,J)
    V4=VO(I,J+2)-4.D0*(VO(I,J+1)+VO(I,J-1))+VO(I,J-2)+6.D0*VO(I,J)
    ADVY=.5D0*(AX1+AX2+AY1+AY2)+(DABS(UIJK)*V2+DABS(VIJK)*V4)*0.25D0

!   VISCOUS TERM(2ND-ORDER CENTRAL SCHEME)
    VISY=((VO(I+1,J)-2.D0*VO(I,J)+VO(I-1,J))*DDXI &
         +(VO(I,J+1)-2.D0*VO(I,J)+VO(I,J-1))*DDYI)*REI

!   POROUS DISK MODEL
    DVEL=DSQRT(UO(I,J)*UO(I,J)+VO(I,J)*VO(I,J))
    SDRAGY=-VO(I,J)*DVEL*DFUNC(I,J)

!   INTERMEDIATE VELOCITY
    VV(I,J)=VO(I,J)+DT*(-ADVY+VISY+SDRAGY)
  END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine vflux


! *******************************************************
  subroutine Poisson_RHS(dt, dx, dy, RHS, UU, VV)

  implicit none
  include 'size.fi'
  double precision, intent(out), dimension(NX,NY)         :: RHS
  double precision, intent(in),  dimension(0:NX+2,0:NY+2) :: UU, VV
  integer          :: i, j
  double precision :: dt, DX, DY, DXI, DYI, ddt

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  ddt = 1.D0/dt

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ddt, DXI, DYI) &
!$OMP PRIVATE(i, j)
!$OMP DO SCHEDULE(static)
  DO J=2,NY
  DO I=2,NX
    RHS(I,J)=-((UU(I,J)-UU(I-1,J))*DXI &
              +(VV(I,J)-VV(I,J-1))*DYI)*ddt
  END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine Poisson_RHS


! *********************************************************************
  subroutine Poisson_AXB(dx, dy, omega, err, RHS, PP)

  implicit none
  include 'size.fi'
  double precision, intent(in),    dimension(NX,NY)         :: RHS
  double precision, intent(inout), dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j
  double precision :: DX, DY, DXI, DYI, DDXI, DDYI
  double precision :: PAGE, dp, err, omega, C1

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  DDXI= DXI/DX
  DDYI= DYI/DY
  C1 = .5D0/(DDXI+DDYI)

!$OMP PARALLEL DO &
!$OMP REDUCTION(+:err) &
!$OMP FIRSTPRIVATE(OMEGA, DDXI, DDYI, C1) &
!$OMP PRIVATE(i, j, PAGE, dp)
  DO J=2,NY
  DO I=2,NX
    PAGE=(PP(I+1,J)+PP(I-1,J))*DDXI &
        +(PP(I,J+1)+PP(I,J-1))*DDYI+RHS(I,J)
    DP= PAGE*C1-PP(I,J)
    ERR=ERR+DP*DP
    PP(I,J)=PP(I,J)+OMEGA*DP
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine Poisson_AXB


! *********************************************************************
  subroutine Poisson_AXB2C(dx, dy, omega, err, RHS, PP)

  implicit none
  include 'size.fi'
  double precision, intent(in),    dimension(NX,NY)         :: RHS
  double precision, intent(inout), dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j, color
  double precision :: DX, DY, DXI, DYI, DDXI, DDYI
  double precision :: PAGE, dp, err, omega, C1

  DXI = 1.D0/DX
  DYI = 1.D0/DY
  DDXI= DXI/DX
  DDYI= DYI/DY
  C1 = .5D0/(DDXI+DDYI)

  do color=0,1

!$OMP PARALLEL DO &
!$OMP REDUCTION(+:err) &
!$OMP FIRSTPRIVATE(OMEGA, DDXI, DDYI, C1, color) &
!$OMP PRIVATE(i, j, PAGE, dp)
    DO J=2,NY
    do i=2+mod(j+color,2), nx, 2
      PAGE=(PP(I+1,J)+PP(I-1,J))*DDXI &
          +(PP(I,J+1)+PP(I,J-1))*DDYI+RHS(I,J)
      DP= PAGE*C1-PP(I,J)
      ERR=ERR+DP*DP
      PP(I,J)=PP(I,J)+OMEGA*DP
    END DO
    END DO
!$OMP END PARALLEL DO

  end do

  return
  end subroutine Poisson_AXB2C


! ********************************
  subroutine Shift_Prs(PP)

  implicit none
  include 'size.fi'
  double precision, intent(inout), dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j
  double precision :: PPP

  PPP=0.D0

!$OMP PARALLEL DO &
!$OMP REDUCTION(+:PPP) &
!$OMP PRIVATE(i, j)
  DO J=1,NY+1
  DO I=1,NX+1
    PPP=PPP+PP(I,J)
  END DO
  END DO
!$OMP END PARALLEL DO

  PPP = PPP/DFLOAT((NX+1)*(NY+1))

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(PPP) &
!$OMP PRIVATE(i, j)
  DO J=1,NY+1
  DO I=1,NX+1
    PP(I,J)=PP(I,J)-PPP
  END DO
  END DO
!$OMP END PARALLEL DO


  return
  end subroutine Shift_Prs


! *****************************
  subroutine BC_Prs(PP)

  implicit none
  include 'size.fi'
  double precision, intent(inout), dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j

!$OMP PARALLEL DO &
!$OMP PRIVATE(i)
  DO I=1,NX+1
    PP(I,1)   =PP(I,2)
    PP(I,NY+1)=PP(I,NY)
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP PRIVATE(j)
  DO J=1,NY+1
    PP(1,J)   =PP(2,J)
    PP(NX+1,J)=PP(NX,J)
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine BC_Prs


! **************************************************
  subroutine Prj_Vel(dt, dx, dy, UU, VV, PP)

  implicit none
  include 'size.fi'
  double precision, intent(inout), dimension(0:NX+2,0:NY+2) :: UU, VV
  double precision, intent(in),    dimension(0:NX+2,0:NY+2) :: PP
  integer          :: i, j
  double precision :: DX, DY, DXI, DYI, dt

  DXI = 1.D0/DX
  DYI = 1.D0/DY

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(DXI, DT) &
!$OMP PRIVATE(i, j)
  DO J=2,NY
  DO I=2,NX-1
    UU(I,J)=UU(I,J)-DT*(PP(I+1,J)-PP(I,J))*DXI
  END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(DYI, DT) &
!$OMP PRIVATE(i, j)
  DO J=2,NY-1
  DO I=2,NX
    VV(I,J)=VV(I,J)-DT*(PP(I,J+1)-PP(I,J))*DYI
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine Prj_Vel


! ***************************************************************************
  subroutine BC_Vel(dt, dx, dy, nstep, istep, UU, VV, UO, UUIN, VVIN)

  implicit none
  include 'size.fi'
  double precision, intent(inout), dimension(0:NX+2,0:NY+2) :: UU, VV
  double precision, intent(in),    dimension(0:NX+2,0:NY+2) :: UO
  double precision, intent(in),    dimension(NSTEP)         :: UUIN, VVIN
  integer          :: i, j, istep, nstep
  double precision :: DX, DY, DXI, dt

  DXI = 1.D0/DX

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(DXI, DT, ISTEP) &
!$OMP PRIVATE(j)
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
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(ISTEP) &
!$OMP PRIVATE(i)
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
!$OMP END PARALLEL DO

  return
  end subroutine BC_Vel


! *************************************
  subroutine stg2Node(U, V, UAVE, VAVE)

  implicit none
  include 'size.fi'
  double precision, intent(in),  dimension(0:NX+2,0:NY+2) :: UAVE, VAVE
  double precision, intent(out), dimension(NX,NY)         :: U, V
  integer          :: i, j

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j)
  DO J=1,NY
  DO I=1,NX
    U(I,J)=.5D0*(UAVE(I,J)+UAVE(I,J+1))
    V(I,J)=.5D0*(VAVE(I,J)+VAVE(I+1,J))
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine stg2Node


! ****************************************
  subroutine accumVars(UU, VV, UAVE, VAVE)

  implicit none
  include 'size.fi'
  double precision, intent(inout), dimension(0:NX+2,0:NY+2) :: UAVE, VAVE
  double precision, intent(in),    dimension(0:NX+2,0:NY+2) :: UU, VV
  integer          :: i, j

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j)
  DO J=1,NY+1
  DO I=1,NX+1
    UAVE(I,J)=UAVE(I,J)+UU(I,J)
    VAVE(I,J)=VAVE(I,J)+VV(I,J)
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine accumVars


! *******************************************************
  subroutine wrt_RCSIns(istep, time, A, U, V, HU, HV, HW)

  implicit none
  include 'size.fi'
  real,             dimension(NX,NY,NZ)      :: HU, HV, HW
  double precision, dimension(NX,NY)         :: U, V
  double precision, dimension(2)             :: A
  integer          :: i, j, k, istep, ii
  double precision :: time
  real             :: t

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(U(I,J))
    HU(I,J,1)=t
    HU(I,J,2)=t
    HU(I,J,3)=t
  END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(V(I,J))
    HV(I,J,1)=t
    HV(I,J,2)=t
    HV(I,J,3)=t
  END DO
  END DO
!$OMP END PARALLEL DO

!   HW(I,J,K)=0.0    initialize at init_Arrays()

  WRITE(10) ISTEP, SNGL(TIME)
  WRITE(10) (SNGL(A(II)),II=1,2)
  WRITE(10) (((HU(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
  WRITE(10) (((HV(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
  WRITE(10) (((HW(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

  return
  end subroutine wrt_RCSIns


! **************************************************************
!  subroutine wrt_VisIns(istep, time, A, HU, HV, HW)

!  implicit none
!  include 'size.fi'
!  real,             dimension(NX,NY,NZ)      :: HU, HV, HW
!  double precision, dimension(2)             :: A
!  integer          :: i, j, k, istep, ii
!  double precision :: time

!  WRITE(10) ISTEP, SNGL(TIME)
!  WRITE(10) (SNGL(A(II)),II=1,2)
!  WRITE(10) (((HU(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
!  WRITE(10) (((HV(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
!  WRITE(10) (((HW(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

!  return
!  end subroutine wrt_VisIns


! *************************************************************************
  subroutine wrt_SPHIns(dx, dy, istep, time, U, V, P, HUVW, HP)

  implicit none
  include 'size.fi'
  real,             dimension(3,NX,NY,NZ)    :: HUVW
  real,             dimension(NX,NY,NZ)      :: HP
  double precision, dimension(NX,NY)         :: U, V, P
  integer          :: i, j, k, istep, svType, dType
  double precision :: time, dx, dy
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

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(P(I,J))
    HP(I,J,1)=t
    HP(I,J,2)=t
    HP(I,J,3)=t
  END DO
  END DO
!$OMP END PARALLEL DO

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

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(u3) &
!$OMP PRIVATE(i, j, u1, u2)
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
!$OMP END PARALLEL DO

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

  return
  end subroutine wrt_SPHIns


! **********************************************
  subroutine avr_Vars(nstep, UAVE, VAVE)

  implicit none
  include 'size.fi'
  double precision, dimension(0:NX+2,0:NY+2) :: UAVE, VAVE
  integer          :: i, j, nstep
  double precision :: dn

  dn = 1.D0/DFLOAT(NSTEP)

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(dn) &
!$OMP PRIVATE(i, j)
  DO J=1,NY+1
  DO I=1,NX+1
    UAVE(I,J)=UAVE(I,J) * dn
    VAVE(I,J)=VAVE(I,J) * dn
!   PAVE(I,J)=PAVE(I,J) * dn
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine avr_Vars


! ************************************
  subroutine wrt_Profile(UAVE)

  implicit none
  include 'size.fi'
  double precision, dimension(0:NX+2,0:NY+2) :: UAVE
  integer          :: j

  OPEN(1,FILE='5D-10D-U-profile.dat',FORM='FORMATTED')
  WRITE(1,'(2F15.7)')(UAVE(1001,j),UAVE(1501,j),j=1,NY)
  CLOSE(1)

  return
  end subroutine wrt_Profile


! **********************************************************************************
  subroutine wrt_RCSAvr(time, nstep, UAVE, VAVE, A, HUAVE, HVAVE, HWAVE)

  implicit none
  include 'size.fi'
  real,             dimension(NX,NY,NZ)      :: HUAVE, HVAVE, HWAVE
  double precision, dimension(0:NX+2,0:NY+2) :: UAVE, VAVE
  double precision, dimension(2)             :: A
  integer          :: i, j, k, ii, nstep
  double precision :: time
  real             :: t

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(UAVE(I,J))
    HUAVE(I,J,1)=t
    HUAVE(I,J,2)=t
    HUAVE(I,J,3)=t
  END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(VAVE(I,J))
    HVAVE(I,J,1)=VAVE(I,J)
    HVAVE(I,J,2)=VAVE(I,J)
    HVAVE(I,J,3)=VAVE(I,J)
  END DO
  END DO
!$OMP END PARALLEL DO

!   HWAVE(I,J,K)=0.0       initialize at init_Arrays()

  open(unit=11,file='3d-tave-vis.dat',form='unformatted')
  write(11) NSTEP,SNGL(TIME)
  write(11) (SNGL(A(II)),II=1,2)
  write(11) (((HUAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  write(11) (((HVAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  write(11) (((HWAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  close(11)

  return
  end subroutine wrt_RCSAvr


! ************************************************************************
  subroutine wrt_SPHAvr(dx, dy, istep, time, UAVE, VAVE, HUVW)

  implicit none
  include 'size.fi'
  real,             dimension(3,NX,NY,NZ)    :: HUVW
  double precision, dimension(0:NX+2,0:NY+2) :: UAVE, VAVE
  double precision, dimension(NX,NY)         :: U, V
  integer          :: i, j, k, istep, svType, dType
  double precision :: time, dx, dy
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

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(u3) &
!$OMP PRIVATE(i, j, u1, u2)
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
!$OMP END PARALLEL DO

  write (s, '("data/uvwa_",I9.9,".sph")') istep
  
  open (unit=22,file=s,form='unformatted')
  write (22) svType, dType
  write (22) nx, ny, nz
  write (22) xorg, yorg, zorg
  write (22) xpch, ypch, zpch
  write (22) istep, tm
  write (22) HUVW
  close (unit=22)

  return
  end subroutine wrt_SPHAvr


! *********************************************************************
! sw(1,*) temporary working
! sw(2,*) accumlated time
! sw(3,*) accumlated count
! sw(4,*) flop/call

subroutine init_watch(N, sw, label)
implicit none
include 'size.fi'
double precision, dimension(4,N)      :: sw
character(20)                         :: label(N)
integer          :: N, i

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
label(15) = "stg2Node"
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

print "('=================================================')"
print "('Elapsed time = ', E15.3)", sum

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