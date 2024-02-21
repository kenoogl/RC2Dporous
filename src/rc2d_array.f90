module array_def
	implicit none
	double precision, dimension(:,:), allocatable	:: RHS, P, DFUNC, UAVE, VAVE
	double precision, dimension(:,:), allocatable	:: PP
	double precision, dimension(:,:), pointer		:: UU, VV, UO, VO, PNT
	real, dimension(:,:,:),           allocatable   :: HUAVE, HVAVE, HWAVE, HU, HV, HW, HP
	real, dimension(:,:,:,:),         allocatable   :: HUVW
	double precision, dimension(:),   allocatable 	:: UUIN, VVIN
	double precision, dimension(:),   allocatable   :: PFUNC
	double precision, dimension(2)           		:: A

contains

	! ************************************************
	subroutine allocArrays(nx, ny, nz, nstep, vmode)
		implicit none
		integer :: nx, ny, nz, nstep, vmode

		allocate(UU(0:NX+2,0:NY+2))
		allocate(VV(0:NX+2,0:NY+2))
		allocate(UO(0:NX+2,0:NY+2))
		allocate(VO(0:NX+2,0:NY+2))
		allocate(PP(0:NX+2,0:NY+2))
		
		allocate(  RHS(NX,NY))
		allocate(    P(NX,NY))
		allocate(DFUNC(NX,NY))
		allocate( UAVE(NX,NY))
		allocate( VAVE(NX,NY))

		allocate(UUIN(NSTEP))
		allocate(VVIN(NSTEP))

		allocate(PFUNC(NY))

		nullify( PNT )

		if (vmode==1) then
			allocate(HUVW(3,NX,NY,NZ))
			allocate(    HP(NX,NY,NZ))
		else
			allocate(HUAVE(NX,NY,NZ))
			allocate(HVAVE(NX,NY,NZ))
			allocate(HWAVE(NX,NY,NZ))
			allocate(   HU(NX,NY,NZ))
			allocate(   HV(NX,NY,NZ))
			allocate(   HW(NX,NY,NZ))
		end if
	end subroutine allocArrays

	! **************************************
	subroutine swapVelocityArray()
		PNT => UU
  		UU  => UO
  		UO  => PNT

  		PNT => VV
  		VV  => VO
  		VO  => PNT
	end subroutine swapVelocityArray

	! **************************************
  	subroutine initArrays(nx, ny, nz, vmode)
  	!$ use omp_lib
  	implicit none
  	integer  :: i, j, k
  	integer  :: nx, ny, nz, vmode

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
!   	PAVE(I,J)=0.D0
  	END DO
  	END DO
	!$OMP END PARALLEL DO


    if (vmode==0) then
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
	end if

  	return
  	end subroutine initArrays

end module array_def
