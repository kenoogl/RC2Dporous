
! *******************************************************
subroutine Poisson_RHS(para, UU, VV)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j
double precision :: DXI, DYI, ddt
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV

include 'param.fi'

DXI = 1.D0/para%DX
DYI = 1.D0/para%DY
ddt = 1.D0/para%dt

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(ddt, DXI, DYI)
DO J=2,NY
DO I=2,NX
    RHS(I,J)=-((UU(I,J)-UU(I-1,J))*DXI &
              +(VV(I,J)-VV(I,J-1))*DYI)*ddt
END DO
END DO
!$OMP END PARALLEL DO

return
end subroutine Poisson_RHS


! ****************************************************************
subroutine Poisson_AXB(para, err)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j
double precision :: DXI, DYI, DDXI, DDYI
double precision :: PAGE, dp, err, omega, C1

include 'param.fi'

DXI = 1.D0/para%DX
DYI = 1.D0/para%DY
DDXI= DXI/para%DX
DDYI= DYI/para%DY
C1 = .5D0/(DDXI+DDYI)
omega = para%omega

!$OMP PARALLEL DO &
!$OMP REDUCTION(+:err) &
!$OMP FIRSTPRIVATE(OMEGA, DDXI, DDYI, C1) &
!$OMP PRIVATE(PAGE, dp)
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


! **************************************************************
subroutine Poisson_AXB2C(para, err)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, color
double precision :: DXI, DYI, DDXI, DDYI
double precision :: PAGE, dp, err, omega, C1

include 'param.fi'

DXI = 1.D0/para%DX
DYI = 1.D0/para%DY
DDXI= DXI/para%DX
DDYI= DYI/para%DY
C1 = .5D0/(DDXI+DDYI)
omega = para%omega

do color=0,1

!$OMP PARALLEL DO &
!$OMP REDUCTION(+:err) &
!$OMP FIRSTPRIVATE(OMEGA, DDXI, DDYI, C1, color) &
!$OMP PRIVATE(PAGE, dp)
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
subroutine Shift_Prs()
!$ use omp_lib
use array_def

implicit none
integer          :: i, j
double precision :: PPP

PPP=0.D0

!$OMP PARALLEL DO &
!$OMP REDUCTION(+:PPP)
DO J=1,NY+1
DO I=1,NX+1
    PPP=PPP+PP(I,J)
END DO
END DO
!$OMP END PARALLEL DO

PPP = PPP/DFLOAT((NX+1)*(NY+1))

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(PPP)
DO J=1,NY+1
DO I=1,NX+1
    PP(I,J)=PP(I,J)-PPP
END DO
END DO
!$OMP END PARALLEL DO


return
end subroutine Shift_Prs


#ifdef _WINDMILL
! *****************************
subroutine BC_Prs_WM()
!$ use omp_lib
use array_def

implicit none
integer          :: i, j

!$OMP PARALLEL DO
DO I=1,NX+1
    PP(I,1)   =PP(I,2)
    PP(I,NY+1)=PP(I,NY)
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO J=1,NY+1
    PP(1,J)   =PP(2,J)
    PP(NX+1,J)=PP(NX,J)
END DO
!$OMP END PARALLEL DO

return
end subroutine BC_Prs_WM
#endif // _WINDMILL


#ifdef _DSL
! *****************************
subroutine BC_Prs_DSL()
!$ use omp_lib
use array_def

implicit none
integer          :: i, j

!$OMP PARALLEL DO
DO I=1,NX+1
    PP(I,1)   =PP(I,NY-1)
    PP(I,NY+1)=PP(I,2)
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO J=1,NY+1
    PP(1,J)   =PP(NX-1,J)
    PP(NX+1,J)=PP(2,J)
END DO
!$OMP END PARALLEL DO

return
end subroutine BC_Prs_DSL
#endif // _DSL
