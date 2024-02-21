program main
!$ use omp_lib
use performance_monitor
use param_def
use array_def

implicit none
double precision :: time, dt, AveragedTime, err
integer          :: ICOUNT, AveragedStep
integer          :: ISTEP, LOOP


!========================================================

call setParams()

TIME = 0.D0
DT   = para%dt
ICOUNT = 0
LOOP = 0
AveragedStep = 0
AveragedTime = 0.0D0

call allocArrays(nx, ny, nz, para%last_step, para%vis_mode)

call init_watch(nx, ny)

!call begin_watch(1)
!call gen_grid()
!call end_watch(1)


call begin_watch(2)
call PDM_profile()
call end_watch(2)


call begin_watch(3)
call Inflow_profile()
call end_watch(3)


call begin_watch(4)
call initArrays(nx, ny, nz, para%vis_mode)
call end_watch(4)


! Restart
IF ( para%start_type == 1 ) THEN
  call begin_watch(5)
  call read_rst_Data(time)
  call end_watch(5)
END IF


!========================================================


OPEN(10,FILE='3d-inst-vis.dat',FORM='UNFORMATTED')
OPEN(12,FILE='history.txt',FORM='FORMATTED')


DO ISTEP=1, para%last_step
  TIME = TIME + DT

  call begin_watch(6)
  call swapVelocityArray()
  call end_watch(6)


  call begin_watch(7)
  call uflux()
  call end_watch(7)


  call begin_watch(8)
  call vflux()
  call end_watch(8)


  call begin_watch(9)
  call Poisson_RHS()
  call end_watch(9)

!--------------------------------------------------------


  DO LOOP=1,para%ITER_max
    ERR=0.D0

    call begin_watch(10)
    !call Poisson_AXB(err)
    call Poisson_AXB2C(err)
    call end_watch(10)

    ERR=DSQRT(ERR/DFLOAT(para%NOXY))
    IF (ERR < para%EPS) exit
  END DO
  

  call begin_watch(11)
  call Shift_Prs()
  call end_watch(11)


  call begin_watch(12)
  call BC_Prs()
  call end_watch(12)

  call begin_watch(13)
  call Prj_Vel()
  call end_watch(13)


  call begin_watch(14)
  call BC_Vel(istep)
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

    !call begin_watch(16)
    !call wrt_RCSIns(istep, time)
    !call end_watch(16)

    call begin_watch(17)
    call wrt_SPHIns(istep, time)
    call end_watch(17)

    ICOUNT=ICOUNT+1
    print "('RC-Scope date=',I10, '/', I10,'(Total)')", ICOUNT, para%last_step/para%Intvl_OutIns
  END IF


  if (TIME > para%AVR_BEGIN) then

    AveragedStep = AveragedStep + 1
    AveragedTime = AveragedTime + DT

    call begin_watch(18)
    call accumVars(AveragedStep)
    call end_watch(18)


    IF (MOD(AveragedStep,para%Intvl_OutAvr)==0) THEN

      call begin_watch(19)
      call wrt_Profile(istep)
      call end_watch(19)

      !call begin_watch(20)
      !call wrt_RCSAvr(AveragedStep, AveragedTime)
      !call end_watch(20)

      call begin_watch(21)
      call wrt_SPHAvr(AveragedStep, AveragedTime)
      call end_watch(21)
    end if

  end if

END DO ! ISTEP

CLOSE(10)
CLOSE(12)


!========================================================


call print_watch()

END program main
