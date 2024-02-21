
! **********************************************
  subroutine accumVars(nadd)
  !$ use omp_lib
  use param_def
  use array_def

  implicit none
  integer          :: i, j, nadd
  double precision :: c1, c2

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
!  use param_def
!  use array_def

!  implicit none
!  integer          :: i, j
!  double precision :: XLEN, YLEN

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