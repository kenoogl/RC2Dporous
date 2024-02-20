
! ************************************************************************************
  subroutine init_Arrays(nstep, UU, VV, UO, VO, PP, UAVE, VAVE, UUIN, VVIN, HW, HWAVE)

  implicit none
  include 'size.fi'
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, PP, UO, VO
  double precision, dimension(NX,NY)         :: UAVE, VAVE
  double precision, dimension(NSTEP)         :: UUIN, VVIN
  real,             dimension(NX,NY,NZ)      :: HW, HWAVE
  integer          :: i, j, k, nstep

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i,j)
  DO J=1,NY
  DO I=1,NX
    UU(I,J)=UUIN(1)
    VV(I,J)=VVIN(1)
    UO(I,J)=UUIN(1) ! ポインタ付け替えをするので、UO, VOにも初期値を入れておく
    VO(I,J)=VVIN(1)
    PP(I,J)=0.D0
  END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i,j)
  DO J=1,NY
  DO I=1,NX
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



! **********************************************
  subroutine accumVars(nadd, UU, VV, UAVE, VAVE)

  implicit none
  include 'size.fi'
  double precision, intent(inout), dimension(NX,NY)         :: UAVE, VAVE
  double precision, intent(in),    dimension(0:NX+2,0:NY+2) :: UU, VV
  integer          :: i, j, nadd
  double precision :: c1, c2

  c1 = 1.0D0 / DBLE(nadd)
  c2 = 1.0D0 - c1

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j)
  DO J=1,NY
  DO I=1,NX
    UAVE(I,J) = c1 * UAVE(I,J) + c2 * (UU(I,J)+UU(I,J+1))*0.5D0
    VAVE(I,J) = c1 * VAVE(I,J) + c2 * (VV(I,J)+VV(I+1,J))*0.5D0
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine accumVars
