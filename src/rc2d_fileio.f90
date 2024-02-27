! ******************************************
subroutine log(m, stp, tm, loop, err)
implicit none
integer          :: m, stp, loop
double precision :: tm, err

write(m,"(I10, F14.6, I7, E14.7)") STP, TM, loop, err

return
end subroutine log


#ifdef _WINDMILL
! ********************************************
subroutine Inflow_profile()
use array_def

implicit none
integer  :: i

OPEN(1,FILE='time-hist-inflow_5.0deg.dat',FORM='FORMATTED')
READ(1,'(2F20.15)') (UUIN(i),VVIN(i),i=1,NSTEP)
CLOSE(1)

!   DO I=1,NSTEP
!     UUIN(I)=1.D0
!     VVIN(I)=0.D0
!   END DO

return
end subroutine Inflow_profile
#endif


! **************************************************
subroutine read_rst_Data(time, UU, VV)
use array_def

implicit none
integer          :: i, j, lx, ly
double precision :: dre, time
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV

OPEN(1,FILE='uvps.dat')
READ(1,*) TIME, dre, lx, ly
READ(1,*) ((UU(I,J),I=0,NX+2),J=0,NY+2)
READ(1,*) ((VV(I,J),I=0,NX+2),J=0,NY+2)
READ(1,*) ((PP(I,J),I=0,NX+2),J=0,NY+2)
CLOSE(1)

return
end subroutine read_rst_Data


#ifdef _RCS
! *********************************************************
subroutine wrt_RCSIns(istep, time, UU, VV)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, k, istep, ii
double precision :: time
real             :: t
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV

!$OMP PARALLEL DO &
!$OMP PRIVATE(t)
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
!$OMP PRIVATE(t)
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
subroutine wrt_VisIns(istep, time)
use array_def

implicit none
integer          :: i, j, k, istep, ii
double precision :: time

WRITE(10) ISTEP, SNGL(TIME)
WRITE(10) (SNGL(A(II)),II=1,2)
WRITE(10) (((HU(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
WRITE(10) (((HV(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
WRITE(10) (((HW(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

return
end subroutine wrt_VisIns
#endif


#ifdef _SPH
! ****************************************************
subroutine wrt_SPHIns(istep, time, para, UU, VV)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, k, istep, svType, dType
double precision :: time, dx, dy
real             :: xorg, yorg, zorg, t, u1, u2, u3
real             :: xpch, ypch, zpch, tm
character(30)    :: s
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV

include 'param.fi'

svType = 1 ! scalar
dType  = 1 ! float 
xorg = -5.0
yorg = -5.0
zorg = 0.0
xpch = SNGL(para%dx)
ypch = SNGL(para%dy)
zpch = SNGL(para%dx)
tm   = SNGL(time)

!$OMP PARALLEL DO &
!$OMP PRIVATE(t)
DO J=1,NY
DO I=1,NX
    t = SNGL(PP(I,J))
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
!$OMP PRIVATE(u1, u2)
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
#endif


#ifdef _RCS
! ******************************************************
subroutine wrt_RCSAvr(step, time)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, k, ii, step
double precision :: time
real             :: t


!$OMP PARALLEL DO &
!$OMP PRIVATE(t)
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
!$OMP PRIVATE(t)
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
write(11) STEP,SNGL(TIME)
write(11) (SNGL(A(II)),II=1,2)
write(11) (((HUAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
write(11) (((HVAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
write(11) (((HWAVE(i,j,k),i=1,nx),j=1,ny),k=1,nz)
close(11)

return
end subroutine wrt_RCSAvr
#endif


#ifdef _SPH
! ************************************************************
subroutine wrt_SPHAvr(step, time, para)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, k, step, svType, dType
double precision :: time
real             :: xorg, yorg, zorg, u1, u2, u3
real             :: xpch, ypch, zpch, tm
character(30)    :: s

include 'param.fi'

svType = 2 ! vector
dType  = 1 ! float 
xorg = -5.0
yorg = -5.0
zorg = 0.0
xpch = SNGL(para%dx)
ypch = SNGL(para%dy)
zpch = SNGL(para%dx)
tm   = SNGL(time)

u3 = 0.0

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(u3) &
!$OMP PRIVATE(u1, u2)
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

write (s, '("data/uvwa_",I9.9,".sph")') step
  
open (unit=22,file=s,form='unformatted')
write (22) svType, dType
write (22) nx, ny, nz
write (22) xorg, yorg, zorg
write (22) xpch, ypch, zpch
write (22) step, tm
write (22) HUVW
close (unit=22)

return
end subroutine wrt_SPHAvr
#endif
