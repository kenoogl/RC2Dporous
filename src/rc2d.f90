program main
!$ use omp_lib

implicit none
include 'size.fi'
integer, PARAMETER :: NSTEP=320000, nWatch=21

double precision, dimension(NX,NY)         :: RHS, P, DFUNC, UAVE, VAVE ! U, V, 
!double precision, dimension(NX,NY)         :: X, Y
double precision, dimension(0:NX+2,0:NY+2) :: PP
double precision, pointer, dimension(:,:)  :: UU, VV, UO, VO, PNT
real,             dimension(NX,NY,NZ)      :: HUAVE, HVAVE, HWAVE, HU, HV, HW, HP
real,             dimension(3,NX,NY,NZ)    :: HUVW
double precision, dimension(NSTEP)         :: UUIN, VVIN
double precision, dimension(NY)            :: PFUNC
double precision, dimension(2)             :: A

double precision, dimension(4,nWatch)      :: stop_watch
character(20)                              :: tm_label(nWatch)

double precision :: omega, eps, time, dt, RE, REI, ERR, CDD
double precision :: DX, DY, XLEN, YLEN, AVR_BEGIN, AveragedTime
integer          :: ITER, ICON, Intvl_Disp, ICOUNT
integer          :: NOXY, ISTEP, LOOP, Intvl_OutAvr, Intvl_OutIns
integer          :: IA, IB, JA, JB, AveragedStep

!========================================================

allocate(UU(0:NX+2,0:NY+2), VV(0:NX+2,0:NY+2), UO(0:NX+2,0:NY+2), VO(0:NX+2,0:NY+2))
nullify( PNT )


call init_watch(nWatch, stop_watch, tm_label)


! Parameters

! SOR
OMEGA=.9D0
EPS  =1.E-3
ITER =100

! TIME STEP ETC
ICON  = 0      ! initial start=0, restart=1
Intvl_OutIns = 4000 ! Interval for vis. data output
Intvl_OutAvr = 4000
Intvl_Disp   = 100    ! Interval for Display
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

AVR_BEGIN = 120.0 ! in case of XLEN=15
AveragedStep = 0
AveragedTime = 0.0D0


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
call init_Arrays(nstep, UU, VV, UO, VO, PP, UAVE, VAVE, UUIN, VVIN, HW, HWAVE)
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
  TIME = TIME + DT

  call begin_watch(6, nWatch, stop_watch)

  ! ポインタの付け替え, サブルーチン実装ではうまくいかない
  PNT => UU
  UU  => UO
  UO  => PNT

  PNT => VV
  VV  => VO
  VO  => PNT
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

! DISPLAY
  IF (MOD(ISTEP,Intvl_Disp)==0) THEN
    call begin_watch(15, nWatch, stop_watch)
    print "('=================================================')"
    print "(' Cd=',F5.1,' Time=',F10.5,'  Istep=',I7)", CDD, TIME, ISTEP
    print "(' Poisson equation iteration=',I7)", LOOP
    print "(' Root-Mean-Square error=',E14.7)", ERR
    call end_watch(15, nWatch, stop_watch)
  END IF
  

! FLOW VISUALIZATION
  IF (MOD(ISTEP,Intvl_OutIns)==0) THEN

    !call begin_watch(16, nWatch, stop_watch)
    !call wrt_RCSIns(istep, time, A, UU, VV, HU, HV, HW)
    !call end_watch(16, nWatch, stop_watch)

    call begin_watch(17, nWatch, stop_watch)
    call wrt_SPHIns(dx, dy, istep, time, UU, VV, PP, HUVW, HP)
    call end_watch(17, nWatch, stop_watch)

    ICOUNT=ICOUNT+1
    print "('RC-Scope date=',I10, '/', I10,'(Total)')", ICOUNT, NSTEP/Intvl_OutIns
  END IF


  if (TIME > AVR_BEGIN) then

    AveragedStep = AveragedStep + 1
    AveragedTime = AveragedTime + DT

    call begin_watch(18, nWatch, stop_watch)
    call accumVars(AveragedStep, UU, VV, UAVE, VAVE)
    call end_watch(18, nWatch, stop_watch)


    IF (MOD(AveragedStep,Intvl_OutAvr)==0) THEN

      call begin_watch(19, nWatch, stop_watch)
      call wrt_Profile(istep, UAVE)
      call end_watch(19, nWatch, stop_watch)

      !call begin_watch(20, nWatch, stop_watch)
      !call wrt_RCSAvr(AveragedTime, AveragedStep, UAVE, VAVE, A, HUAVE, HVAVE, HWAVE)
      !call end_watch(20, nWatch, stop_watch)

      call begin_watch(21, nWatch, stop_watch)
      call wrt_SPHAvr(dx, dy, AveragedStep, AveragedTime, UAVE, VAVE, HUVW)
      call end_watch(21, nWatch, stop_watch)
    end if

  end if

END DO ! ISTEP

CLOSE(10)


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
