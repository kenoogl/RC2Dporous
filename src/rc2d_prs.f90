
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

