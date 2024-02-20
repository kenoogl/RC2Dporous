
! ***********************************************************************
  subroutine PDM_profile(IA, IB, JA, JB, CDD, PFUNC, DFUNC)

  implicit none
  include 'size.fi'
  double precision, dimension(NX,NY)         :: DFUNC
  double precision, dimension(NY)            :: PFUNC
  integer          :: i, j, IA, IB, JA, JB
  double precision :: cdd

  OPEN(1,FILE='pfunc.dat')
  READ(1,'(F15.7)')(PFUNC(J),J=JA,JB)
  CLOSE(1)

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP PRIVATE(i,j)
  DO J=1,NY
  DO I=1,NX
    DFUNC(I,J)=0.D0
  END DO
  END DO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP FIRSTPRIVATE(CDD, IA, IB, JA, JB) &
!$OMP PRIVATE(i,j)
  DO J=JA,JB
  DO I=IA,IB
    DFUNC(I,J)=PFUNC(J)*CDD
  END DO
  END DO
!$OMP END PARALLEL DO

  return
  end subroutine PDM_profile


! ********************************************
  subroutine Inflow_profile(nstep, UUIN, VVIN)

  implicit none
  double precision, dimension(NSTEP)         :: UUIN, VVIN
  integer          :: i, nstep

  OPEN(1,FILE='time-hist-inflow_5.0deg.dat',FORM='FORMATTED')
  READ(1,'(2F20.15)') (UUIN(i),VVIN(i),i=1,NSTEP)
  CLOSE(1)


!   DO I=1,NSTEP
!     UUIN(I)=1.D0
!     VVIN(I)=0.D0
!   END DO

  return
  end subroutine Inflow_profile



! **************************************************
  subroutine read_rst_Data(time, UU, VV, PP)

  implicit none
  include 'size.fi'
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, PP
  integer          :: i, j, lx, ly
  double precision :: dre, time

  OPEN(1,FILE='uvps.dat')
  READ(1,*) TIME, dre, lx, ly
  READ(1,*) ((UU(I,J),I=0,NX+2),J=0,NY+2)
  READ(1,*) ((VV(I,J),I=0,NX+2),J=0,NY+2)
  READ(1,*) ((PP(I,J),I=0,NX+2),J=0,NY+2)
  CLOSE(1)

  return
  end subroutine read_rst_Data



! *********************************************************
  subroutine wrt_RCSIns(istep, time, A, UU, VV, HU, HV, HW)

  implicit none
  include 'size.fi'
  real,             dimension(NX,NY,NZ)      :: HU, HV, HW
  double precision, dimension(0:NX+2,0:NY+2) :: UU, VV
  double precision, dimension(2)             :: A
  integer          :: i, j, k, istep, ii
  double precision :: time
  real             :: t

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(UU(I,J)+UU(I,J+1))*0.5D0
    HU(I,J,1)=t
    HU(I,J,2)=t
    HU(I,J,3)=t
  END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(VV(I,J)+VV(I+1,J))*0.5D0
    HV(I,J,1)=t
    HV(I,J,2)=t
    HV(I,J,3)=t
  END DO
  END DO
!$OMP END PARALLEL DO

!   HW(I,J,K)=0.0    initialize at init_Arrays()

  WRITE(10) ISTEP, SNGL(TIME)
  WRITE(10) (SNGL(A(II)),II=1,2)
  WRITE(10) (((HU(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
  WRITE(10) (((HV(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
  WRITE(10) (((HW(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

  return
  end subroutine wrt_RCSIns


! **************************************************************
!  subroutine wrt_VisIns(istep, time, A, HU, HV, HW)

!  implicit none
!  include 'size.fi'
!  real,             dimension(NX,NY,NZ)      :: HU, HV, HW
!  double precision, dimension(2)             :: A
!  integer          :: i, j, k, istep, ii
!  double precision :: time

!  WRITE(10) ISTEP, SNGL(TIME)
!  WRITE(10) (SNGL(A(II)),II=1,2)
!  WRITE(10) (((HU(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
!  WRITE(10) (((HV(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
!  WRITE(10) (((HW(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

!  return
!  end subroutine wrt_VisIns


! *************************************************************************
  subroutine wrt_SPHIns(dx, dy, istep, time, UU, VV, P, HUVW, HP)

  implicit none
  include 'size.fi'
  double precision, intent(in),  dimension(0:NX+2,0:NY+2) :: UU, VV
  real,             dimension(3,NX,NY,NZ)    :: HUVW
  real,             dimension(NX,NY,NZ)      :: HP
  double precision, dimension(NX,NY)         :: P
  integer          :: i, j, k, istep, svType, dType
  double precision :: time, dx, dy
  real             :: xorg, yorg, zorg, t, u1, u2, u3
  real             :: xpch, ypch, zpch, tm
  character(30)    :: s

  svType = 1 ! scalar
  dType  = 1 ! float 
  xorg = -5.0
  yorg = -5.0
  zorg = 0.0
  xpch = SNGL(dx)
  ypch = SNGL(dy)
  zpch = SNGL(dx)
  tm   = SNGL(time)

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(P(I,J))
    HP(I,J,1)=t
    HP(I,J,2)=t
    HP(I,J,3)=t
  END DO
  END DO
!$OMP END PARALLEL DO

  write (s, '("data/prs_",I9.9,".sph")') istep

  open (unit=22,file=s,form='unformatted')
  write (22) svType, dType
  write (22) nx, ny, nz
  write (22) xorg, yorg, zorg
  write (22) xpch, ypch, zpch
  write (22) istep, tm
  write (22) HP
  close (unit=22)

  u3 = 0.0

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(u3) &
!$OMP PRIVATE(i, j, u1, u2)
  DO J=1,NY
  DO I=1,NX
    u1 = SNGL(UU(I,J)+UU(I,J+1))*0.5D0
    u2 = SNGL(VV(I,J)+VV(I+1,J))*0.5D0
    HUVW(1,I,J,1)=u1
    HUVW(2,I,J,1)=u2
    HUVW(3,I,J,1)=u3
    HUVW(1,I,J,2)=u1
    HUVW(2,I,J,2)=u2
    HUVW(3,I,J,2)=u3
    HUVW(1,I,J,3)=u1
    HUVW(2,I,J,3)=u2
    HUVW(3,I,J,3)=u3
  END DO
  END DO
!$OMP END PARALLEL DO

  svType = 2 ! vector
  write (s, '("data/uvw_",I9.9,".sph")') istep

  open (unit=22,file=s,form='unformatted')
  write (22) svType, dType
  write (22) nx, ny, nz
  write (22) xorg, yorg, zorg
  write (22) xpch, ypch, zpch
  write (22) istep, tm
  write (22) HUVW
  close (unit=22)

  return
  end subroutine wrt_SPHIns


! ****************************
  subroutine wrt_Profile(ts, UAVE)

  implicit none
  include 'size.fi'
  double precision, dimension(NX,NY) :: UAVE
  integer          :: j, ts
  character(50)    :: s

  write (s, '("data/5D-10D-U-profile_",I9.9,".dat")') ts

  OPEN(1,FILE=s,FORM='FORMATTED')
  WRITE(1,'(2F15.7)')(UAVE(1001,j),UAVE(1501,j),j=1,NY)
  CLOSE(1)

  return
  end subroutine wrt_Profile


! **********************************************************************************
  subroutine wrt_RCSAvr(time, nstep, UAVE, VAVE, A, HUAVE, HVAVE, HWAVE)

  implicit none
  include 'size.fi'
  real,             dimension(NX,NY,NZ)      :: HUAVE, HVAVE, HWAVE
  double precision, dimension(NX,NY)         :: UAVE, VAVE
  double precision, dimension(2)             :: A
  integer          :: i, j, k, ii, nstep
  double precision :: time
  real             :: t

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(UAVE(I,J))
    HUAVE(I,J,1)=t
    HUAVE(I,J,2)=t
    HUAVE(I,J,3)=t
  END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP PRIVATE(i, j, t)
  DO J=1,NY
  DO I=1,NX
    t = SNGL(VAVE(I,J))
    HVAVE(I,J,1)=VAVE(I,J)
    HVAVE(I,J,2)=VAVE(I,J)
    HVAVE(I,J,3)=VAVE(I,J)
  END DO
  END DO
!$OMP END PARALLEL DO

!   HWAVE(I,J,K)=0.0       initialize at init_Arrays()

  open(unit=11,file='3d-tave-vis.dat',form='unformatted')
  write(11) NSTEP,SNGL(TIME)
  write(11) (SNGL(A(II)),II=1,2)
  write(11) (((HUAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  write(11) (((HVAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  write(11) (((HWAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  close(11)

  return
  end subroutine wrt_RCSAvr


! ************************************************************
  subroutine wrt_SPHAvr(dx, dy, istep, time, UAVE, VAVE, HUVW)

  implicit none
  include 'size.fi'
  real,             dimension(3,NX,NY,NZ)    :: HUVW
  double precision, dimension(NX,NY)         :: UAVE, VAVE
  double precision, dimension(NX,NY)         :: U, V
  integer          :: i, j, k, istep, svType, dType
  double precision :: time, dx, dy
  real             :: xorg, yorg, zorg, u1, u2, u3
  real             :: xpch, ypch, zpch, tm
  character(30)    :: s

  svType = 2 ! vector
  dType  = 1 ! float 
  xorg = -5.0
  yorg = -5.0
  zorg = 0.0
  xpch = SNGL(dx)
  ypch = SNGL(dy)
  zpch = SNGL(dx)
  tm   = SNGL(time)

  u3 = 0.0

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(u3) &
!$OMP PRIVATE(i, j, u1, u2)
  DO J=1,NY
  DO I=1,NX
    u1 = SNGL(UAVE(I,J))
    u2 = SNGL(VAVE(I,J))
    HUVW(1,I,J,1)=u1
    HUVW(2,I,J,1)=u2
    HUVW(3,I,J,1)=u3
    HUVW(1,I,J,2)=u1
    HUVW(2,I,J,2)=u2
    HUVW(3,I,J,2)=u3
    HUVW(1,I,J,3)=u1
    HUVW(2,I,J,3)=u2
    HUVW(3,I,J,3)=u3
  END DO
  END DO
!$OMP END PARALLEL DO

  write (s, '("data/uvwa_",I9.9,".sph")') istep
  
  open (unit=22,file=s,form='unformatted')
  write (22) svType, dType
  write (22) nx, ny, nz
  write (22) xorg, yorg, zorg
  write (22) xpch, ypch, zpch
  write (22) istep, tm
  write (22) HUVW
  close (unit=22)

  return
  end subroutine wrt_SPHAvr
