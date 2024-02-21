module param_def
    implicit none
    integer :: nx, ny, nz

    type parameter_type
        integer		     :: start_type, last_step
        integer		     :: scheme, problem, iter_max
        integer          :: Intvl_OutAvr, Intvl_OutIns
        integer          :: Intvl_Disp, Intvl_OutHst
        integer          :: IA, IB, JA, JB
        integer          :: NOXY, vis_mode
        double precision :: cfl, omega, thickness, eps
        double precision :: Reynolds, dt, refL, refV
        double precision :: DX, DY, XLEN, YLEN, CDD
        double precision :: AVR_BEGIN
        double precision :: dy_visc
	    character(len=50):: vfile, pfile, xfile
    end type parameter_type
    type (parameter_type):: para


contains
    subroutine setParams()

    implicit none  
    character(len=50)       ::  char

    open (unit=33,file='rcin.txt')

    write(*,*)'Input grid size X'
    read (33,*) nx, char
    write(*,*)'Input grid size Y'
    read (33,*) ny, char
    write(*,*)'Input grid size Z'
    read (33,*) nz, char

    write(*,*)'Start type : 0-Initial 1-Restart'
    read (33,*) para%start_type, char

    write(*,*)'Problem type'
    write(*,*)'1 - Wind farm'
    write(*,*)'2 - DSL Euler'
    read (33,*) para%problem, char

    write(*,*)'Reference Length'
    read (33,*) para%refL, char

    write(*,*)'Reference Velocity'
    read (33,*) para%refV, char

    write(*,*)'Time increment'
    read (33,*) para%dt, char

    write(*,*)'Region size of X'
    read (33,*) para%XLEN, char

    write(*,*)'Region size of Y'
    read (33,*) para%YLEN, char

    write(*,*)'Number of steps to calculate'
    read (33,*) para%last_step, char

    write(*,*)'Interval for printout'
    read (33,*) para%Intvl_Disp, char

    write(*,*)'Interval for history'
    read (33,*) para%Intvl_OutHst, char

    write(*,*)'Interval for instantanenous outfiles'
    read (33,*) para%Intvl_OutIns, char

    write(*,*)'Interval for averaged outfiles'
    read (33,*) para%Intvl_OutAvr, char

    write(*,*)'Relaxation parameter for Poisson'
    read (33,*) para%omega, char

    write(*,*)'Convergence criterion for Poisson' 
    read (33,*) para%eps, char

    write(*,*)'Iteration max for Poisson' 
    read (33,*) para%iter_max, char

    write(*,*)'Start step for averaging'
    read (33,*) para%AVR_BEGIN, char

    write(*,*)'Output file format : 0-RCS, 1-sph'
    read (33,*) para%vis_mode, char

    write(*,*)'Scheme'
    write(*,*)'1 -- Upwind(3) / Euler Explicit(1) / STG'
    write(*,*)'2 -- WCNS(5) / CN(2) / Colocate'
        
    close (unit=33)

    ! given values
    para%dy_visc = 1.5e-5 ! air


    para%DX  = para%XLEN/DFLOAT(NX-1)
    para%DY  = para%YLEN/DFLOAT(NY-1)

    para%reynolds = para%refV* para%refL / para%dy_visc
    para%cfl = para%dt / para%DX

    para%NOXY=(NX-1)*(NY-1)

        
    ! POROUS DISK MODEL
    para%IA = 500
    para%IB = 502
    para%JA = 426
    para%JB = 576
    para%CDD= 13.D0

    ! ===========================================

    write(*,*) ' '
    problem : select case (para%problem)
        case (1)
            write (*,*) 'Problem                 : Wind farm'
        case (2)
            write (*,*) 'Problem                 : Double shear layer'
        case default
            write (*,*) 'No scheme is chosen'
            stop
    end select problem
    write(*,*) ' '

    write(*,"('Calculation region for X : ',E11.4)") para%Xlen
    write(*,"('Calculation region for Y : ',E11.4)") para%Ylen
    write(*,"('Reference Length         : ',E11.4)") para%refL
    write(*,"('Reference Velocity       : ',E11.4)") para%refV
    write(*,"('Dynamic viscosity        : ',E11.4)") para%dy_visc
    write(*,"('Reynolds number          : ',E11.4)") para%reynolds
    write(*,"('Grid spacing for X       : ',E11.4)") para%dx
    write(*,"('Grid spacing for Y       : ',E11.4)") para%dy
    write(*,"('Time increment           : ',E11.4)") para%dt
    write(*,"('Courant number           : ',E11.4)") para%cfl

    if (para%start_type==0) then
        write(*,*) 'Calculation start       : Initial'
    else
        write(*,*) 'Calculation start   　  : Restart'
    endif

    write(*,"('Calculating steps        : ',I11)") para%last_step
    write(*,"('Start step for averaging : ',E11.4)") para%AVR_BEGIN
    
    write(*,"('Interval for printout    : ',I11)") para%Intvl_Disp
    write(*,"('Interval for history     : ',I11)") para%Intvl_OutHst
    write(*,"('Interval for Inst.files  : ',I11)") para%Intvl_OutIns
    write(*,"('Interval for Avr.files   : ',I11)") para%Intvl_OutAvr

    write(*,"('Relax. coef. for Poisson : ',E11.4)") para%omega
    write(*,"('Convergence  for Poisson : ',E11.4)") para%eps
    write(*,"('Iter. max.   for Poisson : ',I11)") para%iter_max
    
    if (para%vis_mode==0) then
        write(*,*) 'File format for output  : RCS'
    else
        write(*,*) 'File format for output  : sph'
    end if

    write(*,"('Coef. for wind turbine   : ',E11.4)") para%cdd

    !write(*,*) 'Thickness of layer    :',para%thickness

    end subroutine setParams

end module param_def
