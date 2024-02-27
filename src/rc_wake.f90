#ifdef _WINDMILL

! **********************************************************
subroutine findWMindex(para)

implicit none
include 'size.fi'
integer          :: k, getindex, ia, ib, ja, jb
double precision :: p1, p2, p3, p4, d, c, x0, y0, xc, yc, dh

include 'param.fi'

x0 = para%x0
y0 = para%y0
dh = para%dx

do k=1,para%NOWM
  xc = para%pwm(1,k) ! xc
  yc = para%pwm(2,k) ! yc
  d  = para%pwm(3,k) ! diameter
  c  = para%pwm(4,k) ! cord
  p1 = xc
  p2 = xc + c
  p3 = yc - 0.5D0 * d
  p4 = yc + 0.5D0 * d
  
  ia = getIndex(x0, p1, dh)
  ib = getIndex(x0, p2, dh)
  ja = getIndex(y0, p3, dh)
  jb = getIndex(y0, p4, dh)

  para%IDXWM(1,k) = ia
  para%IDXWM(2,k) = ib
  para%IDXWM(3,k) = ja
  para%IDXWM(4,k) = jb

  write(*,"('Wind turbine no. ',I3)") k
  write(*,"('IA=',I5,' IB=',I5, ' JA=',I5,' JB=',I5)") para%IDXWM(1,k), para%IDXWM(2,k), &
                                                       para%IDXWM(3,k), para%IDXWM(4,k)
  if (ia<1 .or. NX<ia) then
    write(*,*) 'IA is out of range'
    stop
  end if
  if (ib<1 .or. NX<ib) then
    write(*,*) 'IB is out of range'
    stop
  end if
  if (ja<1 .or. NY<ja) then
    write(*,*) 'JA is out of range'
    stop
  end if
  if (jb<1 .or. NY<jb) then
    write(*,*) 'JB is out of range'
    stop
  end if

end do

para%NOBD = para%IDXWM(4,1) - para%IDXWM(3,1) + 1
write(*,*) 'Number of Grid for blade', para%NOBD
write(*,*) ' '

return
end subroutine findWMindex


! **********************************************************
integer function getIndex(x0, p, dh)
implicit none
double precision :: x0, p, dh

getIndex = int((p-x0)/dh)
end function


! **********************************************************
subroutine gen_PDMprofile(NOBD, PFUNC)
implicit none
integer          :: j, NOBD
double precision :: pi
double precision, dimension(NOBD)  :: PFUNC

pi = 2.0D0 * asin(1.0D0)

do j=1,NOBD
  PFUNC(j) = 0.5D0*( cos(2.0D0*pi*DBLE(j-1)/DBLE(NOBD-1)-pi)+1.D0 )
end do

return
end subroutine gen_PDMprofile


! **********************************************************
subroutine set_DFUNC(para, NOBD, PFUNC)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, m, IA, IB, JA, JB, NOBD
double precision :: cdd
double precision, dimension(NOBD)  :: PFUNC

include 'param.fi'

do m=1, para%NOWM

ia  = para%IDXWM(1,m)
ib  = para%IDXWM(2,m)
ja  = para%IDXWM(3,m)
jb  = para%IDXWM(4,m)
cdd = para%pwm(5,m)


!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(CDD, IA, IB, JA, JB)
DO J=JA,JB
DO I=IA,IB
  DFUNC(I,J)=PFUNC(J-JA+1)*CDD
END DO
END DO
!$OMP END PARALLEL DO

end do

return
end subroutine set_DFUNC


! ****************************
subroutine wrt_Profile(ts)
use array_def

implicit none
integer          :: j, ts
character(50)    :: s

write (s, '("data/5D-10D-U-profile_",I9.9,".dat")') ts

OPEN(1,FILE=s,FORM='FORMATTED')
WRITE(1,'(2F15.7)')(UAVE(1001,j),UAVE(1501,j),j=1,NY)
CLOSE(1)

return
end subroutine wrt_Profile

#endif // _WINDMILL


#if 0
! **********************************************************
subroutine PDM_profile(para)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, JA, JB
double precision :: cdd

include 'param.fi'

OPEN(1,FILE='pfunc.dat')
READ(1,'(F15.7)')(PFUNC(J),J=JA,JB)
CLOSE(1)

!$OMP PARALLEL DO
DO J=1,NY
DO I=1,NX
  DFUNC(I,J)=0.D0
END DO
END DO
!$OMP END PARALLEL DO

return
end subroutine PDM_profile
#endif