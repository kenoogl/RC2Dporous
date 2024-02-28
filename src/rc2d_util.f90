
! **************************************
subroutine initArrays(UU, VV, UO, VO)
!$ use omp_lib
use array_def

implicit none
integer  :: i, j, k
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, UO, VO

!$OMP PARALLEL DO
DO J=1,NY
DO I=1,NX
    UU(I,J)=UUIN(1)
    VV(I,J)=VVIN(1)
    UO(I,J)=UUIN(1) ! ポインタ付け替えをするので、
    VO(I,J)=VVIN(1) ! UO, VOにも初期値を入れておく
    PP(I,J)=0.D0
END DO
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO J=1,NY
DO I=1,NX
    UAVE(I,J)=0.D0
    VAVE(I,J)=0.D0
!   PAVE(I,J)=0.D0
END DO
END DO
!$OMP END PARALLEL DO


#ifndef _SPH
    !$OMP PARALLEL DO
    DO K=1,NZ
    DO J=1,NY
    DO I=1,NX
        HW(I,J,K)=0.0
        HWAVE(I,J,K)=0.0
    END DO
    END DO
    END DO
    !$OMP END PARALLEL DO
#endif

return
end subroutine initArrays



! ************************************************
subroutine setParams(para)

implicit none
include 'size.fi'

integer            :: k, pchk
double precision   :: tmp, t0, refL
character(len=50)  :: char

include 'param.fi'


open (unit=33,file='rcin.txt')

write(*,*)'parameter check : 0-No, 1-Yes'
read (33,*) pchk, char

write(*,*)'Start type : 0-Initial 1-Restart'
read (33,*) para%start_type, char

write(*,*)'Reference Length [m]'
read (33,*) para%refL, char
refL = para%refL

write(*,*)'Reference Velocity [m/s]'
read (33,*) para%refV, char

write(*,*)'Origin of Region [m]'
read (33,*) para%x0, para%y0, char

write(*,*)'Counrant number'
read (33,*) para%courant, char

write(*,*)'Region size of X [m]'
read (33,*) para%XLEN, char

write(*,*)'Region size of Y [m]'
read (33,*) para%YLEN, char

!write(*,*)'Number of steps to calculate'
!read (33,*) para%last_step, char

write(*,*)'Interval for printout'
read (33,*) para%Intvl_Disp, char

write(*,*)'Interval for history'
read (33,*) para%Intvl_OutHst, char

write(*,*)'Interval for instantanenous outfiles'
read (33,*) para%Intvl_OutIns, char

write(*,*)'Interval for averaged outfiles'
read (33,*) para%Intvl_OutAvr, char

write(*,*)'Relaxation parameter for Poisson'
read (33,*) para%omega, char

write(*,*)'Convergence criterion for Poisson' 
read (33,*) para%eps, char

write(*,*)'Iteration max for Poisson' 
read (33,*) para%iter_max, char

write(*,*)'Start time for averaging [sec]'
read (33,*) para%AVR_BEGIN, char


write(*,*)'Scheme'
write(*,*)'1 -- Upwind(3) / Euler Explicit(1) / STG'
write(*,*)'2 -- WCNS(5) / CN(2) / Collocated'
read (33,*) para%scheme, char


#ifdef _WINDMILL
  write(*,*)'No. of wind mill'
  read (33,*) para%NOWM, char

  if (para%nowm>256) then
    write(*,*)'No. of wind mill must be no more than 256'
    stop
  end if

  write(*,*)'Profile of PDM : 0-Flat, 1-Gaussian, 2-Double Gaussian'
  read (33,*) para%PDMprofile, char

  do k=1,para%NOWM
    write(*,*)'wind mill property : xc[m], yc[m], d[m], cord[m], coef[-]'
    read (33,*) para%pwm(1,k), para%pwm(2,k), para%pwm(3,k), para%pwm(4,k), para%pwm(5,k), char
    para%pwm(1,k) = para%pwm(1,k) / refL
    para%pwm(2,k) = para%pwm(2,k) / refL
    para%pwm(3,k) = para%pwm(3,k) / refL
    para%pwm(4,k) = para%pwm(4,k) / refL
  end do
#endif

#ifdef _DSL
  write(*,*) 'Thickness of DSL layer'
  read (33,*) para%thickness
#endif
        
close (unit=33)


! given values
para%dy_visc = 1.5e-5 ! air
para%last_step = nstep
t0 = refL / para%refV

para%x0 = para%x0 / refL
para%y0 = para%y0 / refL
para%XLEN = para%XLEN / refL
para%YLEN = para%YLEN / refL

para%DX  = para%XLEN/DFLOAT(NX-1)
para%DY  = para%YLEN/DFLOAT(NY-1)
if (para%DX /= para%DY) then
  write(*,*) 'dx and dy is not same'
  stop
end if

para%reynolds = para%refV* refL / para%dy_visc
para%dt = para%courant * para%DX
para%AVR_BEGIN = para%AVR_BEGIN / t0

para%NOXY=(NX-1)*(NY-1)


! POROUS DISK MODEL
!para%IA = 500
!para%IB = 502
!para%JA = 426
!para%JB = 576
!para%CDD= 13.D0

! ===========================================

write(*,*) ' '
write(*,*) '========================================================='
#ifdef _WINDMILL
    write (*,*) 'Problem                 : Wind farm'
#endif

#ifdef _DSL
    write (*,*) 'Problem                 : Double shear layer'
#endif
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
  write(*,*) 'Calculation start       : Initial'
else
  write(*,*) 'Calculation start       : Restart'
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
  write(*,*) 'Scheme                  : Upwind(3) / Euler Explicit(1) / STG'
else
  write(*,*) 'Scheme                  : WCNS(5) / CN(2) / Collocated'
endif
    
#ifdef _SPH
  write(*,*) 'File format for output  : sph'
#else
  write(*,*) 'File format for output  : RCS'
#endif

write(*,*) ' '
write(*,*) 'Problem dependent parameters'

#ifdef _WINDMILL
  write(*,*) 'Number of Windmill      :', para%NOWM
  if (para%PDMprofile==0) then
    write(*,*)'Profile of PDM          : Flat'
  else if (para%PDMprofile==1) then
    write(*,*)'Profile of PDM          : Gaussian'
  else
    write(*,*)'Profile of PDM          : Double Gaussian'
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
      write(*,*) 'Blade diameter is not same : No.', k 
    end if
  end do
  
#endif

#ifdef _DSL
  write(*,*) 'Thickness of layer    :',para%thickness
#endif

write(*,*) ' '
if (pchk==1) then
  write(*,*) 'Parameter check mode & stop'
  stop
end if

end subroutine setParams


! **********************************************
subroutine accumVars(nadd, UU, VV)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, nadd
double precision :: c1, c2
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV


c1 = 1.0D0 / DBLE(nadd)
c2 = 1.0D0 - c1

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(c1, c2)
DO J=1,NY
DO I=1,NX
  UAVE(I,J) = c2 * UAVE(I,J) + c1 * (UU(I,J)+UU(I,J+1))*0.5D0
  VAVE(I,J) = c2 * VAVE(I,J) + c1 * (VV(I,J)+VV(I+1,J))*0.5D0
END DO
END DO
!$OMP END PARALLEL DO

return
end subroutine accumVars


! ******************************************
!  subroutine gen_grid()
!  !$ use omp_lib
!  use array_def

!  implicit none
!  include 'size.fi'
!  integer          :: i, j, nx, ny
!  double precision :: XLEN, YLEN
!include 'param.fi'

!  xlen = para%xlen
!  ylen = para%ylen

! !$OMP PARALLEL DO &
! !$OMP FIRSTPRIVATE(XLEN, YLEN)
!  DO J=1,NY
!  DO I=1,NX
!    X(I,J)=XLEN/DFLOAT(NX-1)*DFLOAT(I-1)-5.D0
!    Y(I,J)=YLEN/DFLOAT(NY-1)*DFLOAT(J-1)-5.D0
!  END DO
!  END DO
! !$OMP END PARALLEL DO

!  return
!  end subroutine gen_grid