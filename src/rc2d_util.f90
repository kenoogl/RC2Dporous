
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