
subroutine parse_jsonfile(para, input_file)
use json_module
use json_module,  IK => json_IK, &
                  RK => json_RK, &
                  CK => json_CK, &
                  LK => json_LK

implicit none
type(json_file)          :: jsonfile
type(json_value),pointer :: p          !! a pointer for low-level manipulations
type(json_core)          :: core       !! factory for manipulating `json_value` pointers

! 読み込んだ値を格納する変数
integer(IK)                           :: ival    ! 整数型
integer(IK), allocatable              :: ia(:)   ! 整数型配列
real(RK)                              :: rval    ! 実数型
real(RK), allocatable                 :: ra(:)   ! 実数型配列
character(kind=CK,len=:), allocatable :: cval    ! 文字列
character(kind=CK,len=:), allocatable :: ca(:)   ! 文字列型配列
logical(LK)                           :: lval    ! 論理型
character(kind=CK,len=10)             :: str
integer, parameter :: io_unit = 6                ! 標準出力用の装置番号
integer            :: k, pchk, i, ecount
double precision   :: tmp, t0, refL
character(len=50)  :: char
CHARACTER(LEN=64)  :: input_file

include 'param.fi'
include 'size.fi'


! ---------------------------given values
para%dy_visc = 1.5e-5 ! air
para%last_step = nstep
para%NOXY=(NX-1)*(NY-1)
para%problem = 0

#ifdef _WINDMILL
	para%problem = 1
#elif defined _DSL
	para%problem = 2
#endif


pchk = 0 ! no check


! -------------------------
write(*,*) 'JSON-Fortran version: '//json_fortran_version()
write(*,"(' parsing file : ',A)") trim(input_file)

! ファイルの読込
call jsonfile%initialize()
call jsonfile%load(filename = input_file)

if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  stop
endif

write(io_unit,'(A)') ''

ecount = 0


! ------------------------- [dry_run]
call jsonfile%get('dry_run.char', cval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else
  ! write(io_unit,*) 'Dry run                  = '//trim(cval)
  if (trim(cval) == 'yes') pchk = 1
end if


! ------------------------- [start]
call jsonfile%get('start.char', cval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else
  ! write(io_unit,*) 'start                    = '//trim(cval)
  if (trim(cval)=='initial') then
    para%start_type = 0
  else if (trim(cval)=='restart') then
    para%start_type = 1
  else
  	write(*,*) 'keyword error for start'
  	ecount = ecount + 1
  endif
end if


! ------------------------- [Reference_Length]
call jsonfile%get('Reference_Length.real', rval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) 'Reference Length [m]     = ', rval
  para%refL = rval
  refL = rval
end if


! ------------------------- [Reference_Velocity]
call jsonfile%get('Reference_Velocity.real', rval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) 'Reference Velocity [m/s] = ', rval
  para%refV = rval
end if


t0 = refL / para%refV
para%reynolds = para%refV* refL / para%dy_visc



! ------------------------- [Computational_Region]
call jsonfile%get('Computational_Region.real array', ra)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) 'Computational Region [m] = ', ra
  para%Xlen = ra(1) / refL
  para%Ylen = ra(2) / refL
endif

para%DX  = para%XLEN/DFLOAT(NX-1)
para%DY  = para%YLEN/DFLOAT(NY-1)

if (para%DX /= para%DY) then
  write(*,*) 'dx and dy is not same'
  ecount = ecount + 1
end if


! ------------------------- [Origin_of_Region]
call jsonfile%get('Origin_of_Region.real array', ra)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) 'Origin of Region [m]     = ', ra
  para%x0 = ra(1) / refL
  para%y0 = ra(2) / refL
endif


! ------------------------- [Courant_number]
call jsonfile%get('Courant_number.real', rval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) 'Courant number           = ', rval
  para%courant = rval
end if


para%dt = para%courant * para%DX


! ------------------------- [Intervals]
call jsonfile%get('Intervals.display.int', ival)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) '  for display            = ', ival
  para%Intvl_Disp = ival
endif
    
call jsonfile%get('Intervals.history.int', ival)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) '  for history            = ', ival
  para%Intvl_OutHst = ival
endif

call jsonfile%get('Intervals.Instantaneous_file.int', ival)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) '  for Instantaneous file = ', ival
  para%Intvl_OutIns = ival
endif

call jsonfile%get('Intervals.averaged_file.int', ival)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) '  for averaged file      = ', ival
  para%Intvl_OutAvr = ival
endif


! ------------------------- [Poisson_parameter]
call jsonfile%get('Poisson_parameter.coef_acceleration.real', rval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) '  coef. acceleration     = ', rval
  para%omega = rval
endif

call jsonfile%get('Poisson_parameter.convergence_criteria.real', rval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) '  convergence_criteria   = ', rval
  para%eps = rval
endif

call jsonfile%get('Poisson_parameter.Iteration_max.int', ival)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) '  Iteration max          = ', ival
  para%iter_max = ival
endif


! ------------------------- [Start_time_for_averaging]
call jsonfile%get('Start_time_for_averaging.real', rval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) 'Start time for averaging = ', rval
  para%AVR_BEGIN = rval / t0
end if


! ------------------------- [Convection_Scheme]
call jsonfile%get('Convection_Scheme.char', cval)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else
  ! write(io_unit,*) 'Convection Scheme        = '//trim(cval)
  if (trim(cval)=='3rd_upwind') then
    para%scheme = 1
  else if (trim(cval)=='5th_wcns') then
  	para%scheme = 2
  else
    write(*,*) 'keyword error for Convection_Scheme'
  	ecount = ecount + 1
  endif
end if


! ------------------------- [turbine]
if (para%problem == 1 ) then

call jsonfile%get('Number_of_turbines.int', ival)
if (jsonfile%failed()) then
  call jsonfile%print_error_message(io_unit)
  ecount = ecount + 1
else 
  ! write(io_unit, *) 'Number of turbines       = ', ival
  para%nowm = ival
end if

do i=1,para%nowm
  write(str,fmt='(I3)') i
  str = adjustl(str)

  call jsonfile%get('turbine('//trim(str)//').PDM_profile.char', cval)
  if (jsonfile%failed()) then
    call jsonfile%print_error_message(io_unit)
    ecount = ecount + 1
  else
    ! write(io_unit,'(A)') ' PDM_profile            = ' //trim(cval)
    if (trim(cval)=='flat') then
      para%PDMprofile = 1
    else if (trim(cval)=='gaussian') then
    	para%PDMprofile = 2
    else if (trim(cval)=='double_gaussian') then
    	para%PDMprofile = 3
    else
    	write(*,*) 'keyword error for PDM_profile'
  		ecount = ecount + 1
    endif
  endif

  call jsonfile%get('turbine('//trim(str)//').center_position.real array', ra)
  if (jsonfile%failed()) then
    call jsonfile%print_error_message(io_unit)
    ecount = ecount + 1
  else
    ! write(io_unit,'(A,*(F30.16,1X))') ' center position         = ',ra
    para%pwm(1,i) = ra(1) / refL
    para%pwm(2,i) = ra(2) / refL
  end if

  call jsonfile%get('turbine('//trim(str)//').blade_diameter.real', rval)
  if (jsonfile%failed()) then
    call jsonfile%print_error_message(io_unit)
    ecount = ecount + 1
  else 
    ! write(io_unit, *) ' Blade Diameter         = ', rval
    para%pwm(3,i) = rval / refL
  end if

  call jsonfile%get('turbine('//trim(str)//').cord_length.real', rval)
  if (jsonfile%failed()) then
    call jsonfile%print_error_message(io_unit)
    ecount = ecount + 1
  else 
    ! write(io_unit, *) ' Cord length           = ', rval
    para%pwm(4,i) = rval / refL
  end if

  call jsonfile%get('turbine('//trim(str)//').pressure_loss_coef.real', rval)
  if (jsonfile%failed()) then
    call jsonfile%print_error_message(io_unit)
    ecount = ecount + 1
  else 
    ! write(io_unit, *) ' Pressure Loss coef.   = ', rval
    para%pwm(5,i) = rval
  end if

end do

end if ! WINDMILL


if (para%problem == 2 ) then
  call jsonfile%get('Thickness_of_DSL.real', rval)
  if (jsonfile%failed()) then
    call jsonfile%print_error_message(io_unit)
    ecount = ecount + 1
  else 
    ! write(io_unit, *) ' Blade Diameter         = ', rval
    para%thickness = rval
  end if
endif

call jsonfile%destroy()

!--------------------------------------

write(*,*) '========================================================='
if (para%problem == 1 ) then
  write (*,*) 'Problem                 : Wind farm'
else if (para%problem == 2 ) then
  write (*,*) 'Problem                 : Double shear layer'
endif

write(*,*) ' '

write(*,"('Grid size                : ',3I6)") nx, ny, nz
write(*,"('Reference Length         : ',E11.4,' [m]')") refL
write(*,"('Reference Velocity       : ',E11.4,' [m/s]')") para%refV
write(*,"('Origin of Region         : (',E11.4,' , ',E11.4') [m] (',E11.4,' , ',E11.4,') [-]')") &
                                      para%x0 * refL, para%y0 * refL, para%x0, para%y0
write(*,"('Calculation region for X : ',E11.4,' [m] ',E11.4,'[-]')") para%Xlen * refL, para%Xlen
write(*,"('Calculation region for Y : ',E11.4,' [m] ',E11.4,'[-]')") para%Ylen * refL, para%Ylen
write(*,"('Dynamic viscosity        : ',E11.4,' [m^2/s]')") para%dy_visc
write(*,"('Reynolds number          : ',E11.4)") para%reynolds
write(*,"('Grid spacing for X       : ',E11.4,' [m] ',E11.4,'[-]')") para%dx * refL, para%dx
write(*,"('Grid spacing for Y       : ',E11.4,' [m] ',E11.4,'[-]')") para%dy * refL, para%dy
write(*,"('Courant number           : ',E11.4)") para%courant
write(*,"('Time increment           : ',E11.4,' [sec] ',E11.4,'[-]')") para%dt * t0, para%dt

if (para%start_type==0) then
  write(*,'(A)') 'Calculation start        :  Initial'
else
  write(*,'(A)') 'Calculation start        :  Restart'
endif

tmp = para%last_step * para%dt * t0
write(*,"('Calculating steps        : ',E11.4,' [sec] ',I11,'[step]')") tmp , para%last_step
write(*,"('Start time for averaging : ',E11.4,' [sec] ',E11.4,'[-]')") para%AVR_BEGIN * t0, para%AVR_BEGIN
tmp = para%Intvl_Disp * para%dt * t0
write(*,"('Interval for printout    : ',E11.4,' [sec] ',I11,'[step]')") tmp , para%Intvl_Disp
tmp = para%Intvl_OutHst * para%dt * t0
write(*,"('Interval for history     : ',E11.4,' [sec] ',I11,'[step]')") tmp , para%Intvl_OutHst
tmp = para%Intvl_OutIns * para%dt * t0
write(*,"('Interval for Inst.files  : ',E11.4,' [sec] ',I11,'[step]')") tmp , para%Intvl_OutIns
tmp = para%Intvl_OutAvr * para%dt * t0
write(*,"('Interval for Avr.files   : ',E11.4,' [sec] ',I11,'[step]')") tmp , para%Intvl_OutAvr

write(*,"('Relax. coef. for Poisson : ',E11.4)") para%omega
write(*,"('Convergence  for Poisson : ',E11.4)") para%eps
write(*,"('Iter. max.   for Poisson : ',I11)") para%iter_max

if (para%scheme==1) then
  write(*,'(A)') 'Scheme                   : Upwind(3) / Euler Explicit(1) / STG'
else
  write(*,'(A)') 'Scheme                   : WCNS(5) / CN(2) / Collocated'
endif


#ifdef _SPH
  write(*,'(A)') 'File format for output   : sph'
#else
  write(*,'(A)') 'File format for output   : RCS'
#endif


write(*,*) ' '
write(*,'(A)') 'Problem dependent parameters'

if (para%problem == 1 ) then
  write(*,'(A)') 'Number of Windmill       :', para%NOWM

  if (para%nowm>256) then
    write(*,*)'No. of wind mill must be no more than 256'
    stop
  end if

  if (para%PDMprofile==1) then
    write(*,'(A)')'Profile of PDM          : Flat'
  else if (para%PDMprofile==2) then
    write(*,'(A)')'Profile of PDM          : Gaussian'
  else if (para%PDMprofile==3) then
    write(*,'(A)')'Profile of PDM          : Double Gaussian'
  end if

  do k=1, para%NOWM
    write(*,"('Wind turbine no.=',I3)") k
    write(*,*) '                             x[m]      y[m] bladeD[m]   cord[m]   coef[-]'
    write(*,"('     Position & coef  : ',5E10.3)") para%pwm(1,k)*refL, para%pwm(2,k)*refL, &
                                              para%pwm(3,k)*refL, para%pwm(4,k)*refL, para%pwm(5,k)
    write(*,"('     Non-dimensional  : ',5E10.3)") para%pwm(1,k), para%pwm(2,k), para%pwm(3,k), &
                                              para%pwm(4,k), para%pwm(5,k)
    ! check
    if (para%pwm(3,1) /= para%pwm(3,k)) then
      write(*,'(A)') 'Blade diameter is not same : No.', k 
    end if
  end do
  
end if ! para%problem

write(*,*) ' '


if (para%problem == 2 ) then
  write(*,*) 'Thickness of layer      :', para%thickness
endif


! -------------------------------

if (ecount>0) then
  write(*,*) 'json parameter file is something wrong.'
  stop
endif

if (pchk==1) then
	write(*,*) 'Dry run : only parameter check and stop'
  stop
endif

end subroutine parse_jsonfile
