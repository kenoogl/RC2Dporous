
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

character(len=50)  :: char

include 'param.fi'


open (unit=33,file='rcin.txt')

write(*,*)'Start type : 0-Initial 1-Restart'
read (33,*) para%start_type, char

write(*,*)'Problem type'
write(*,*)'1 - Wind farm'
write(*,*)'2 - DSL Euler'
read (33,*) para%problem, char

write(*,*)'Reference Length'
read (33,*) para%refL, char

write(*,*)'Reference Velocity'
read (33,*) para%refV, char

write(*,*)'Time increment'
read (33,*) para%dt, char

write(*,*)'Region size of X'
read (33,*) para%XLEN, char

write(*,*)'Region size of Y'
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

write(*,*)'Start step for averaging'
read (33,*) para%AVR_BEGIN, char


write(*,*)'Scheme'
write(*,*)'1 -- Upwind(3) / Euler Explicit(1) / STG'
write(*,*)'2 -- WCNS(5) / CN(2) / Colocate'
        
close (unit=33)

! given values
para%dy_visc = 1.5e-5 ! air
para%last_step = nstep

para%DX  = para%XLEN/DFLOAT(NX-1)
para%DY  = para%YLEN/DFLOAT(NY-1)

para%reynolds = para%refV* para%refL / para%dy_visc
para%cfl = para%dt / para%DX

para%NOXY=(NX-1)*(NY-1)

        
! POROUS DISK MODEL
para%IA = 500
para%IB = 502
para%JA = 426
para%JB = 576
para%CDD= 13.D0

! ===========================================

write(*,*) ' '
problem : select case (para%problem)
  case (1)
    write (*,*) 'Problem                 : Wind farm'
  case (2)
    write (*,*) 'Problem                 : Double shear layer'
  case default
    write (*,*) 'No scheme is chosen'
    stop
end select problem
write(*,*) ' '

write(*,"('Grid size                : ',3I6)") nx, ny, nz
write(*,"('Calculating steps        : ',I11)") para%last_step
write(*,"('Calculation region for X : ',E11.4)") para%Xlen
write(*,"('Calculation region for Y : ',E11.4)") para%Ylen
write(*,"('Reference Length         : ',E11.4)") para%refL
write(*,"('Reference Velocity       : ',E11.4)") para%refV
write(*,"('Dynamic viscosity        : ',E11.4)") para%dy_visc
write(*,"('Reynolds number          : ',E11.4)") para%reynolds
write(*,"('Grid spacing for X       : ',E11.4)") para%dx
write(*,"('Grid spacing for Y       : ',E11.4)") para%dy
write(*,"('Time increment           : ',E11.4)") para%dt
write(*,"('Courant number           : ',E11.4)") para%cfl

if (para%start_type==0) then
  write(*,*) 'Calculation start       : Initial'
else
  write(*,*) 'Calculation start   　  : Restart'
endif

write(*,"('Start step for averaging : ',E11.4)") para%AVR_BEGIN
    
write(*,"('Interval for printout    : ',I11)") para%Intvl_Disp
write(*,"('Interval for history     : ',I11)") para%Intvl_OutHst
write(*,"('Interval for Inst.files  : ',I11)") para%Intvl_OutIns
write(*,"('Interval for Avr.files   : ',I11)") para%Intvl_OutAvr

write(*,"('Relax. coef. for Poisson : ',E11.4)") para%omega
write(*,"('Convergence  for Poisson : ',E11.4)") para%eps
write(*,"('Iter. max.   for Poisson : ',I11)") para%iter_max
    
#ifdef _SPH
  write(*,*) 'File format for output  : sph'
#else
  write(*,*) 'File format for output  : RCS'
#endif

write(*,"('Coef. for wind turbine   : ',E11.4)") para%cdd

!write(*,*) 'Thickness of layer    :',para%thickness
write(*,*) ' '

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