
!================================
module array_def
implicit none
include 'size.fi'
double precision, dimension(NX,NY)         :: RHS, DFUNC, UAVE, VAVE
double precision, dimension(0:NX+2,0:NY+2) :: PP

#ifdef _WINDMILL
double precision, dimension(NSTEP)         :: UUIN, VVIN
#endif

#ifdef _SPH
real, dimension(3,NX,NY,NZ)                :: HUVW
real, dimension(NX,NY,NZ)                  :: HP
#else
real, dimension(NX,NY,NZ)                  :: HUAVE, HVAVE, HWAVE, HU, HV, HW
double precision, dimension(2)             :: A
#endif

end module array_def
!================================

program main
!$ use omp_lib
use array_def
use performance_monitor

implicit none
double precision :: time, dt, AveragedTime, err
integer          :: ICOUNT, AveragedStep
integer          :: ISTEP, LOOP
double precision, dimension(:,:), pointer :: UU, VV, UO, VO, PNT

#ifdef _WINDMILL
double precision, dimension(:),pointer    :: PFUNC
#endif

include 'param.fi'

!========================================================

allocate(UU(0:NX+2,0:NY+2), VV(0:NX+2,0:NY+2), UO(0:NX+2,0:NY+2), VO(0:NX+2,0:NY+2))
nullify( PNT )

call setParams(para)

#ifdef _WINDMILL
  call findWMindex(para)
  allocate(PFUNC(para%NOBD))
#endif

TIME = 0.D0
DT   = para%dt
ICOUNT = 0
LOOP = 0
AveragedStep = 0
AveragedTime = 0.0D0


call init_watch(nx, ny)

!call begin_watch(1)
!call gen_grid(para)
!call end_watch(1)


#ifdef _WINDMILL
call begin_watch(2)
!call PDM_profile(para)
call gen_PDMprofile(para%NOBD, PFUNC)
call set_DFUNC(para, para%NOBD, PFUNC)
call end_watch(2)


call begin_watch(3)
call Inflow_profile()
call end_watch(3)
#endif

call begin_watch(4)
call initArrays(UU, VV, UO, VO)
call end_watch(4)


! Restart
IF ( para%start_type == 1 ) THEN
  call begin_watch(5)
  call read_rst_Data(time, UU, VV)
  call end_watch(5)
END IF


!========================================================

#ifdef _RCS
OPEN(10,FILE='3d-inst-vis.dat',FORM='UNFORMATTED')
#endif

OPEN(12,FILE='history.txt',FORM='FORMATTED')


DO ISTEP=1, para%last_step
  TIME = TIME + DT

  call begin_watch(6) ! ポインタの付け替え, サブルーチン実装ではうまくいかない
  PNT => UU
  UU  => UO
  UO  => PNT

  PNT => VV
  VV  => VO
  VO  => PNT
  call end_watch(6)


  call begin_watch(7)
  call uflux(para, UU, VV, UO, VO)
  call end_watch(7)


  call begin_watch(8)
  call vflux(para, UU, VV, UO, VO)
  call end_watch(8)


  call begin_watch(9)
  call Poisson_RHS(para, UU, VV)
  call end_watch(9)

!--------------------------------------------------------


  DO LOOP=1,para%ITER_max
    ERR=0.D0

    call begin_watch(10)
    !call Poisson_AXB(para, err)
    call Poisson_AXB2C(para, err)
    call end_watch(10)

    ERR=DSQRT(ERR/DFLOAT(para%NOXY))
    IF (ERR < para%EPS) exit
  END DO
  

  call begin_watch(11)
  call Shift_Prs()
  call end_watch(11)


  call begin_watch(12)
#ifdef _WINDMILL
  call BC_Prs_WM()
#else
  call BC_Prs_DSL()
#endif
  call end_watch(12)

  call begin_watch(13)
  call Prj_Vel(para, UU, VV)
  call end_watch(13)


  call begin_watch(14)
#ifdef _WINDMILL
  call BC_Vel_WM(istep, para, UU, VV, UO)
#else
  call BC_Vel_DSL(para, UU, VV)
#endif
  call end_watch(14)
  

!--------------------------------------------------------

! DISPLAY
  IF (MOD(ISTEP, para%Intvl_Disp)==0) THEN
    call begin_watch(15)
    call log(6, istep, time, loop, err)
    call end_watch(15)
  END IF
  

! History
  IF (MOD(ISTEP, para%Intvl_OutHst)==0) THEN
    call begin_watch(22)
    call log(12, istep, time, loop, err)
    call end_watch(22)
  END IF


! FLOW VISUALIZATION
  IF (MOD(ISTEP, para%Intvl_OutIns)==0) THEN

#ifdef _SPH
    call begin_watch(17)
    call wrt_SPHIns(istep, time, para, UU, VV)
    call end_watch(17)
#else
    call begin_watch(16)
    call wrt_RCSIns(istep, time, UU, VV)
    call end_watch(16)
#endif

    ICOUNT=ICOUNT+1
    print "('RC-Scope date=',I10, '/', I10,'(Total)')", ICOUNT, para%last_step/para%Intvl_OutIns
  END IF


  if (TIME > para%AVR_BEGIN) then

    AveragedStep = AveragedStep + 1
    AveragedTime = AveragedTime + DT

    call begin_watch(18)
    call accumVars(AveragedStep, UU, VV)
    call end_watch(18)

    IF (MOD(AveragedStep,para%Intvl_OutAvr)==0) THEN

      call begin_watch(19)
      call wrt_Profile(istep)
      call end_watch(19)

#ifdef _SPH
      call begin_watch(21)
      call wrt_SPHAvr(AveragedStep, AveragedTime, para)
      call end_watch(21)
#else
      call begin_watch(20)
      call wrt_RCSAvr(AveragedStep, AveragedTime)
      call end_watch(20)
#endif
            
    end if
  end if

END DO ! ISTEP

#ifdef _RCS
CLOSE(10)
#endif

CLOSE(12)


!========================================================


call print_watch()

END program main
