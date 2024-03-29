
module performance_monitor
	implicit none
	integer, PARAMETER :: nWatch=22

	double precision, dimension(4,nWatch), private	:: sw
	integer, dimension(nWatch)      	 , private  :: idx
	double precision, dimension(nWatch)  , private  :: val
	character(30)                    	 , private  :: label(nWatch)

	! sw(1,*) temporary working
	! sw(2,*) accumlated time
	! sw(3,*) accumlated count
	! sw(4,*) flop/call

	public init_watch
	public begin_watch
	public end_watch
	public print_watch
	private swap_watch

contains

	! *********************************
	subroutine init_watch(nx, ny)
	implicit none
	integer          :: i, nx, ny

	do i=1,nWatch
  		sw(1,i) = 0.0
  		sw(2,i) = 0.0
  		sw(3,i) = 0.0
  		sw(4,i) = 0.0
	end do

	label(1) = 'meshing'
	label(2) = "setup_PDM"
	label(3) = "read_inflow"
	label(4) = "init_arrays"
	label(5) = "read_rst_file"
	label(6) = "str_Nstep"
	label(7) = "Uflux"
	label(8) = "Vflux"
	label(9) = "Poi_RHS"
	label(10) = "Poi_AXB"
	label(11) = "Shift_Prs"
	label(12) = "BC_Prs"
	label(13) = "Update_Vec"
	label(14) = "BC_Vec"
	label(15) = "display"
	label(16) = "wrt_RCSIns"
	label(17) = "wrt_SPHIns"
	label(18) = "accumAvr"
	label(19) = "wrt_Profile"
	label(20) = "wrt_RCSAvr"
	label(21) = "wrt_SPHAvr"
	label(22) = "History"

	sw(4,1) = 6.0*NX*NY
	sw(4,2) = 0.0
	sw(4,3) = 0.0
	sw(4,4) = 0.0
	sw(4,5) = 0.0
	sw(4,6) = 0.0
	sw(4,7) = (NY-1)*(NX-2)*107.0
	sw(4,8) = (NY-2)*(NX-1)*107.0
	sw(4,9) = (NY-1)*(NX-1)*7.0
	sw(4,10) = (NY-1)*(NX-1)*14.0
	sw(4,11) = (NY+1)*(NX+1)*4.0
	sw(4,12) = 0.0
	sw(4,13) = (NY-1)*(NX-2)*4.0+(NY-2)*(NX-1)*4.0
	sw(4,14) = (NY-1)*8.0
	sw(4,15) = 0.0
	sw(4,16) = NX*NY*4.0
	sw(4,17) = NX*NY*4.0
	sw(4,18) = NX*NY*10.0
	sw(4,19) = 0.0
	sw(4,20) = 0.0
	sw(4,21) = 0.0
	sw(4,22) = 0.0

	end subroutine init_watch


	! *********************************
	subroutine begin_watch(key)
	!$ use omp_lib
	implicit none

	integer          :: key

	!$ sw(1,key) = omp_get_wtime()

	end subroutine begin_watch


	! *********************************
	subroutine end_watch(key)
	!$ use omp_lib
	implicit none

	integer          :: key
	double precision :: tm_ed, tm_st

	!$ tm_ed = omp_get_wtime()
	tm_st = sw(1,key)
	sw(2,key) = sw(2,key) + (tm_ed - tm_st)
	sw(3,key) = sw(3,key) + 1.0

	end subroutine end_watch


	! *********************************
	subroutine print_watch()
	implicit none

	integer          :: i, j
	double precision :: sum

	! to avoid obtaining false Gflops
	do i=1,nWatch
  		if (sw(2,i)<1.0e-3) sw(2,i)=1.0e-3
	end do

	! copy sort target and set index
	do i=1,nWatch
 		idx(i)=i
 		val(i)=sw(2,i)
	end do

	sum =0.0
	do i=1,nWatch
  		sum = sum + val(i)
	end do

	! obtain sorted index
	do i=1,nWatch-1
  		do j=i+1,nWatch
    		if ( val(i)<val(j) ) call swap_watch(i,j)
  		end do
	end do


	print "(A21,A7,A10,A15,A15,A12,A10)", 'Label', '%Time', 'Call', 'accTime', 'avrTime', 'flop/call', 'Gflops'


	do j=1,nWatch
  		i = idx(j)
  		print "(A20,F8.2,I10,2E15.5,E12.3,F10.3)", &
       	label(i), &
       	sw(2,i)/sum*100.0, &
       	int(sw(3,i)), &
       	sw(2,i), &
       	sw(2,i)/sw(3,i), &
       	sw(4,i), &
       	sw(4,i)*sw(3,i)/sw(2,i)*1.0E-9
	end do

	print "('=================================================')"
	print "('Elapsed time = ', E15.5)", sum

	end subroutine print_watch


	! *********************************
	subroutine swap_watch(i, j)
	implicit none

	integer          :: i,j,m
	double precision :: a

	a = val(i)
	m = idx(i)

	val(i) = val(j)
	idx(i) = idx(j)

	val(j) = a
	idx(j) = m

	end subroutine swap_watch

end module performance_monitor














