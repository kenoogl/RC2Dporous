! ****************************************************
subroutine uflux(para, UU, VV, UO, VO)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j
double precision :: CU1, CU2, CV1, CV2
double precision :: AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL
double precision :: UU1, UU2, UU3, UU4, U2, U4, ADVX, VISX, SDRAGX
double precision :: dt, rei, DXI, DYI, DDXI, DDYI
double precision :: C1, C2
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, UO, VO

include 'param.fi'

dt  = para%dt
rei = 1.D0/para%Reynolds
DXI = 1.D0/para%DX
DYI = 1.D0/para%DY
DDXI= DXI/para%DX
DDYI= DYI/para%DY
C1 = 1.D0/16.D0
C2 = 1.D0/24.D0
SDRAGX = 0.D0

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(dt, rei, DXI, DYI, DDXI, DDYI, C1, C2) &
!$OMP PRIVATE(CU1, CU2, CV1, CV2) &
!$OMP PRIVATE(AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL) &
!$OMP PRIVATE(UU1, UU2, UU3, UU4, U2, U4, ADVX, VISX, SDRAGX)
DO J=2,NY
DO I=2,NX-1

    ! CONVECTIVE TERM(4TH-ORDER CENTRAL SCHEME)
    CU1=-UO(I+1,J)+9.D0*(UO(I,J)+UO(I-1,J))-UO(I-2,J)
    CU2=-UO(I+2,J)+9.D0*(UO(I+1,J)+UO(I,J))-UO(I-1,J)
    CV1=-VO(I+2,J-1)+9.D0*(VO(I+1,J-1)+VO(I,J-1))-VO(I-1,J-1)
    CV2=-VO(I+2,J)+9.D0*(VO(I+1,J)+VO(I,J))-VO(I-1,J)
    UU1=-UO(I+1,J)+27.D0*(UO(I,J)-UO(I-1,J))+UO(I-2,J)
    UU2=-UO(I+2,J)+27.D0*(UO(I+1,J)-UO(I,J))+UO(I-1,J)
    UU3=-UO(I,J+1)+27.D0*(UO(I,J)-UO(I,J-1))+UO(I,J-2)
    UU4=-UO(I,J+2)+27.D0*(UO(I,J+1)-UO(I,J))+UO(I,J-1)
    AX1=(CU1*C1)*(UU1*C2*DXI)
    AX2=(CU2*C1)*(UU2*C2*DXI)
    AY1=(CV1*C1)*(UU3*C2*DYI)
    AY2=(CV2*C1)*(UU4*C2*DYI)

!   4TH-ORDER NUMERICAL VISCOSITY(K-K SCHEME TYPE)
    UIJK=UO(I,J)*DXI
    VIJK=.25D0*(VO(I+1,J-1)+VO(I,J-1)+VO(I+1,J)+VO(I,J))*DYI
    U2=UO(I+2,J)-4.D0*(UO(I+1,J)+UO(I-1,J))+UO(I-2,J)+6.D0*UO(I,J)
    U4=UO(I,J+2)-4.D0*(UO(I,J+1)+UO(I,J-1))+UO(I,J-2)+6.D0*UO(I,J)
    ADVX=.5D0*(AX1+AX2+AY1+AY2)+(DABS(UIJK)*U2+DABS(VIJK)*U4)*0.25D0

!   VISCOUS TERM(2ND-ORDER CENTRAL SCHEME)
    VISX=((UO(I+1,J)-2.D0*UO(I,J)+UO(I-1,J))*DDXI &
         +(UO(I,J+1)-2.D0*UO(I,J)+UO(I,J-1))*DDYI)*REI

#ifdef _WINDMILL
!   POROUS DISK MODEL
    DVEL=DSQRT(UO(I,J)*UO(I,J)+VO(I,J)*VO(I,J))
    SDRAGX=-UO(I,J)*DVEL*DFUNC(I,J)
#endif

!   INTERMEDIATE VELOCITY
    UU(I,J)=UO(I,J)+DT*(-ADVX+VISX+SDRAGX)
END DO
END DO
!$OMP END PARALLEL DO

return
end subroutine uflux


! ****************************************************
subroutine vflux(para, UU, VV, UO, VO)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j
double precision :: CU1, CU2, CV1, CV2
double precision :: AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL
double precision :: VV1, VV2, VV3, VV4, V2, V4, ADVY, VISY, SDRAGY
double precision :: dt, rei, DXI, DYI, DDXI, DDYI
double precision :: C1, C2
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, UO, VO

include 'param.fi'

dt  = para%dt
rei = 1.D0/para%Reynolds
DXI = 1.D0/para%DX
DYI = 1.D0/para%DY
DDXI= DXI/para%DX
DDYI= DYI/para%DY
C1 = 1.D0/16.D0
C2 = 1.D0/24.D0
SDRAGY = 0.D0

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(dt, rei, DXI, DYI, DDXI, DDYI, C1, C2) &
!$OMP PRIVATE(CU1, CU2, CV1, CV2) &
!$OMP PRIVATE(AX1, AX2, AY1, AY2, UIJK, VIJK, DVEL) &
!$OMP PRIVATE(VV1, VV2, VV3, VV4, V2, V4, ADVY, VISY, SDRAGY)
DO J=2,NY-1
DO I=2,NX

!   CONVECTIVE TERM(4TH-ORDER CENTRAL SCHEME)
    CU1=-UO(I-1,J+2)+9.D0*(UO(I-1,J+1)+UO(I-1,J))-UO(I-1,J-1)
    CU2=-UO(I,J+2)+9.D0*(UO(I,J+1)+UO(I,J))-UO(I,J-1)
    CV1=-VO(I,J+1)+9.D0*(VO(I,J)+VO(I,J-1))-VO(I,J-2)
    CV2=-VO(I,J+2)+9.D0*(VO(I,J+1)+VO(I,J))-VO(I,J-1)
    VV1=-VO(I+1,J)+27.D0*(VO(I,J)-VO(I-1,J))+VO(I-2,J)
    VV2=-VO(I+2,J)+27.D0*(VO(I+1,J)-VO(I,J))+VO(I-1,J)
    VV3=-VO(I,J+1)+27.D0*(VO(I,J)-VO(I,J-1))+VO(I,J-2)
    VV4=-VO(I,J+2)+27.D0*(VO(I,J+1)-VO(I,J))+VO(I,J-1)
    AX1=(CU1*C1)*(VV1*C2*DXI)
    AX2=(CU2*C1)*(VV2*C2*DXI)
    AY1=(CV1*C1)*(VV3*C2*DYI)
    AY2=(CV2*C1)*(VV4*C2*DYI)

!   4TH-ORDER NUMERICAL VISCOSITY(K-K SCHEME TYPE)
    UIJK=.25D0*(UO(I-1,J+1)+UO(I-1,J)+UO(I,J+1)+UO(I,J))*DXI
    VIJK=VO(I,J)*DYI
    V2=VO(I+2,J)-4.D0*(VO(I+1,J)+VO(I-1,J))+VO(I-2,J)+6.D0*VO(I,J)
    V4=VO(I,J+2)-4.D0*(VO(I,J+1)+VO(I,J-1))+VO(I,J-2)+6.D0*VO(I,J)
    ADVY=.5D0*(AX1+AX2+AY1+AY2)+(DABS(UIJK)*V2+DABS(VIJK)*V4)*0.25D0

!   VISCOUS TERM(2ND-ORDER CENTRAL SCHEME)
    VISY=((VO(I+1,J)-2.D0*VO(I,J)+VO(I-1,J))*DDXI &
         +(VO(I,J+1)-2.D0*VO(I,J)+VO(I,J-1))*DDYI)*REI

#ifdef _WINDMILL
!   POROUS DISK MODEL
    DVEL=DSQRT(UO(I,J)*UO(I,J)+VO(I,J)*VO(I,J))
    SDRAGY=-VO(I,J)*DVEL*DFUNC(I,J)
#endif

!   INTERMEDIATE VELOCITY
    VV(I,J)=VO(I,J)+DT*(-ADVY+VISY+SDRAGY)
END DO
END DO
!$OMP END PARALLEL DO

return
end subroutine vflux



! **************************************************
subroutine Prj_Vel(para, UU, VV)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j
double precision :: DXI, DYI, dt
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV

include 'param.fi'

DXI = 1.D0/para%DX
DYI = 1.D0/para%DY
dt  = para%dt

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(DXI, DT)
DO J=2,NY
DO I=2,NX-1
    UU(I,J)=UU(I,J)-DT*(PP(I+1,J)-PP(I,J))*DXI
END DO
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(DYI, DT)
DO J=2,NY-1
DO I=2,NX
    VV(I,J)=VV(I,J)-DT*(PP(I,J+1)-PP(I,J))*DYI
END DO
END DO
!$OMP END PARALLEL DO

return
end subroutine Prj_Vel


#ifdef _WINDMILL
! ******************************************************
subroutine BC_Vel_WM(istep, para, UU, VV, UO)
!$ use omp_lib
use array_def

implicit none
integer          :: i, j, istep
double precision :: dtx, s1, s2
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV, UO

include 'param.fi'

dtx  = para%dt / para%DX

!------------ original
#if 0
!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(dtx, ISTEP)
DO J=2,NY
    UU(1,J)=UUIN(ISTEP)
    UU(0,J)=UUIN(ISTEP)
    VV(1,J)=VVIN(ISTEP)
    VV(0,J)=VVIN(ISTEP)
    UU(NX,J)  =UO(NX,J)-dtx*(UO(NX,J)-UO(NX-1,J))
    VV(NX+1,J)=VV(NX,J)
    VV(NX+2,J)=2.D0*VV(NX+1,J)-VV(NX,J)
    UU(NX+1,J)=2.D0*UU(NX,J)-UU(NX-1,J)
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(ISTEP)
DO I=0,NX+2
    UU(I,1)=UUIN(ISTEP)
    UU(I,0)=UUIN(ISTEP)
    VV(I,1)=VVIN(ISTEP)
    VV(I,0)=VVIN(ISTEP)
    UU(I,NY+1)=UUIN(ISTEP)
    UU(I,NY+2)=UUIN(ISTEP)
    VV(I,NY)  =VVIN(ISTEP)
    VV(I,NY+1)=VVIN(ISTEP)
END DO
!$OMP END PARALLEL DO
#endif
!------------


! Traction free

!$OMP PARALLEL DO FIRSTPRIVATE(dtx, ISTEP) private(s1, s2)
DO J=2,NY
  UU(1,J) = UUIN(ISTEP)
  UU(0,J) = UUIN(ISTEP)
  VV(1,J) = VVIN(ISTEP)
  VV(0,J) = VVIN(ISTEP)

  s1 = UO(NX,J)-dtx*(UO(NX,J)-UO(NX-1,J))
  UU(NX  ,J) = s1
  UU(NX+1,J) = 2.D0*s1-UU(NX-1,J)

  s2 = VV(NX,J)
  VV(NX+1,J) = s2
  VV(NX+2,J) = s2
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO private(s1, s2)
DO I=2,NX
  s1 = VV(I,2   ) + UU(I,2 ) - UU(I-1,2)
  s2 = VV(I,NY-1) - UU(i,NY) + UU(i-1,NY)
  VV(I,1   ) = s1
  VV(I,0   ) = s1
  VV(I,NY  ) = s2
  VV(I,NY+1) = s2
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO private(s1, s2)
DO I=2,NX-1
  s1 = UU(i,2 ) + VV(i+1,1 ) - VV(i,1 )
  s2 = UU(i,NY) - VV(i+1,NY) + VV(i,NY)
  UU(I,0)    = s1
  UU(I,1)    = s1
  UU(I,NY+1) = s2
  UU(I,NY+2) = s2
END DO
!$OMP END PARALLEL DO


return
end subroutine BC_Vel_WM
#endif // _WINDMILL


#ifdef _DSL
! ******************************************************
subroutine BC_Vel_DSL(para, UU, VV)
!$ use omp_lib

implicit none
integer          :: i, j
double precision :: DXI, dt
double precision, dimension(0:NX+2,0:NY+2) :: UU, VV

include 'param.fi'

DXI = 1.D0/para%DX
dt  = para%dt

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(DXI, DT, ISTEP)
DO J=2,NY
    UU(1,J)=UU(IX-2,j)
    UU(0,J)=UU(IX-3,j)
    VV(1,J)=VV(i,JX-2)
    VV(0,J)=VV(i,JX-3)
    UU(NX,J)  = UU(2,j)
    VV(NX+1,J)= VV(3,J)
    VV(NX+2,J)= VV(4,J)
    UU(NX+1,J)= UU(3,J)
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(ISTEP)
DO I=0,NX+2
    UU(I,1)=UU(i,NY-2)
    UU(I,0)=UU(i,NY-3)
    VV(I,1)=VV(i,NY-2)
    VV(I,0)=VV(i,NY-3)
    UU(I,NY+1)=UU(i,2)
    UU(I,NY+2)=UU(i,3)
    VV(I,NY)  =VV(i,2)
    VV(I,NY+1)=VV(i,3)
END DO
!$OMP END PARALLEL DO

return
end subroutine BC_Vel_DSL
#endif // _DSL
