!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                  matshiftapplyfunc.f90              !!!
!!!                                                     !!!
!!!  This file is part of the R package 'matrixutils'.  !!!
!!!  This code can be used to calculate predefined      !!!
!!!  functions on the overlap of two shifted matrices.  !!!
!!!                                                     !!!
!!!  written by: Yann Buechau <nobodyinperson@gmx.de>   !!!
!!!            last update: 30.01.2016                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! module for variables !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vars
	implicit none
	
	! copy of passed variables
	integer :: NROWBASE,NCOLBASE,NROWSHIFT,NCOLSHIFT,NROWRES,NCOLRES ! dimensions
	
	! non-passed variables
	logical :: VERBOSE = .FALSE. ! Be verbose?
	logical :: DEBUG   = .FALSE. ! Perform several checks to debug
	integer :: j,k   = 0 ! loop indices
	
	! variables for shiftbounds
	integer :: rowshift,colshift ! current shifts
	integer :: mirosh,marosh,micosh,macosh,miroba,maroba,micoba,macoba ! shift bounds
	integer :: rowdiff,coldiff,rowdiffceil,rowdifffloor,coldiffceil,coldifffloor ! row/col diff variables
	integer :: noverlap ! number of overlapping points
	
	! other global variables
	integer :: NAVAL  = -99999 ! Caution: predefined!

end module vars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine to calculate function by shifts !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine matshiftapplyfunc(MBASE,MBROW,MBCOL,MSHIFT,MSROW,MSCOL,MRES,MRROW,MRCOL,ROWSHIFTS,COLSHIFTS,FUNC,RADIAL,NA)
	use vars
	implicit none
	
	! passed arguments
	logical          :: RADIAL ! Only shift to circular shape?
	integer          :: FUNC ! function nr. to apply to overlap
	integer          :: NA   ! Value of missing values
	integer          :: MBROW,MBCOL,MSROW,MSCOL,MRROW,MRCOL ! dimensions
	integer          :: ROWSHIFTS(MSROW),COLSHIFTS(MSCOL)   ! vectors of row- and colshifts 
	double precision :: MBASE (MBROW,MBCOL)     ! base matrix
	double precision :: MSHIFT(MSROW,MSCOL)     ! shift matrix
	double precision :: MRES  (MRROW,MRCOL)     ! resulting matrix
	
	
	! functions
	double precision :: rmse ! RMSE function
	double precision :: multav ! multiplication and averaging
	
	! test output of passed variables
	if( VERBOSE ) then
		print *,"mbase:" 
		do j=1,MBROW
			print *, MBASE(j,1:MBCOL)
		end do
		print *,"mshift:" 
		do j=1,MSROW
			print *, MSHIFT(j,1:MSCOL)
		end do
		print *,
		print ('(6A8)'),"MBROW","MBCOL","MSROW","MSCOL","MRROW","MRCOL"
		print ('(6I8)'),MBROW,MBCOL,MSROW,MSCOL,MRROW,MRCOL
	end if
	
	! copy passed variables to make them available to the module
	NROWBASE  = MBROW
	NCOLBASE  = MBCOL
	NROWSHIFT = MSROW
	NCOLSHIFT = MSCOL
	NAVAL     = NA
	

	! fill the resulting matrix
	do k=1,MRCOL
		colshift     = COLSHIFTS(k)    
		do j=1,MRROW
			rowshift  = ROWSHIFTS(j)

			if(RADIAL) then ! if RADIAL, check if now outside of circular region
				if( (2.*real(rowshift)/real(MBROW))**2 + (2.*real(colshift)/real(MBCOL))**2 .gt. 1. ) then
					MRES(j,k) = NAVAL ! outside areas filled with NAs
					cycle
				end if
			end if
	
			if(VERBOSE) then ! verbose output
				print *,"--------------------------------------"
				print '(2A4)',"j","k"
				print '(2I4)',j,k
			end if
	
			call shiftbounds ! calculate shiftbounds
			
			! Function to be applied to the overlap
			select case (FUNC)
			case(1)
				MRES(j,k) = rmse  ( MBASE(miroba:maroba,micoba:macoba), MSHIFT(mirosh:marosh,micosh:macosh), noverlap )
			case(2)
				MRES(j,k) = multav( MBASE(miroba:maroba,micoba:macoba), MSHIFT(mirosh:marosh,micosh:macosh), noverlap )
			end select
	
			if(VERBOSE) print *,"MRES(j,k):",MRES(j,k)
		end do
	end do
	
	if( VERBOSE ) then ! verbose output matrix
		print *,"mres:" 
		do j=1,MRROW
				print *, MRES(j,1:MRCOL)
		end do
	end if
end subroutine matshiftapplyfunc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine to calculate shift bounds !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine shiftbounds
	use vars
	implicit none
	
	! helper variables
	rowdiff = NROWBASE - NROWSHIFT
	coldiff = NCOLBASE - NCOLSHIFT
	rowdiffceil  = ceiling ( real(rowdiff) / 2. )
	rowdifffloor = floor   ( real(rowdiff) / 2. )
	coldiffceil  = ceiling ( real(coldiff) / 2. )
	coldifffloor = floor   ( real(coldiff) / 2. )
	
	! shift bounds
	mirosh = max(         1,         1 - ( rowshift  + rowdifffloor ) )
	marosh = min( NROWSHIFT, NROWSHIFT - ( rowshift  - rowdiffceil  ) )
	micosh = max(         1,         1 - ( colshift  + coldifffloor ) ) 
	macosh = min( NCOLSHIFT, NCOLSHIFT - ( colshift  - coldiffceil  ) ) 
	miroba = max(         1,         1 + ( rowshift  + rowdifffloor ) ) 
	maroba = min(  NROWBASE,  NROWBASE + ( rowshift  - rowdiffceil  ) ) 
	micoba = max(         1,         1 + ( colshift  + coldifffloor ) ) 
	macoba = min(  NCOLBASE,  NCOLBASE + ( colshift  - coldiffceil  ) ) 
	
	! number of remaining points in the overlap
	noverlap = ( maroba - miroba + 1 ) * ( macoba - micoba + 1 )
	if(.NOT. DEBUG) then ! Only check this at debugging
		if(noverlap .ne. (marosh-mirosh+1)*(macosh-micosh+1)) print *,"CAUTION: Something with the indices is wrong!!!"
	end if
	
	if(VERBOSE) then ! verbose output
		print ('(2A10)'),"rowshift","colshift"
		print ('(2I10)'),rowshift,colshift
		print ('(2A10)'),"rowdiff","coldiff"
		print ('(2I10)'),rowdiff,coldiff
		print ('(4A14)'),"rowdiffceil","rowdifffloor","coldiffceil","coldifffloor"
		print ('(4I14)'),rowdiffceil,rowdifffloor,coldiffceil,coldifffloor
		print ('(8A8)'),"mirosh","marosh","micosh","macosh","miroba","maroba","micoba","macoba"
		print ('(8I8)'),mirosh,marosh,micosh,macosh,miroba,maroba,micoba,macoba
		print ('(A10)'),"noverlap"
		print ('(I10)'),noverlap
	end if
	
end subroutine shiftbounds


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Root Mean Square Error !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! x,y : values
!!! n   : number of values
double precision function rmse(x,y,n) 
	use vars
	implicit none

	integer         , intent(in) :: n ! number of values
	double precision, intent(in) :: x(n),y(n) ! input
	integer                      :: i = 0 ! loop variabe
	integer                      :: skip = 0  ! count of skipped values

	skip = 0   ! initialize (yes, double initialization)
	rmse = 0.0 ! initialize

	do i=1,n ! Loop over all given values
		if(int(x(i)) .eq. NAVAL .OR. int(y(i)) .eq. NAVAL) then
			skip = skip + 1 ! Count up number of missing values
			cycle ! skip missing values
		end if
		rmse = rmse + ( x(i) - y(i) ) ** 2
	end do

	if(skip .ne. n) then ! Return value
		rmse = sqrt( rmse / dble(n-skip) )
	else
		rmse = NAVAL
	end if
end function rmse


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Simple multiplication and averaging !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! x,y : values
!!! n   : number of values
double precision function multav(x,y,n) 
	use vars
	implicit none

	integer         , intent(in) :: n ! number of values
	double precision, intent(in) :: x(n),y(n) ! input
	integer                      :: i = 0 ! loop variabe
	integer                      :: skip = 0  ! count of skipped values

	skip = 0   ! initialize (yes, double initialization)
	multav = 0.0 ! initialize

	do i=1,n ! Loop over all given values
		if(int(x(i)) .eq. NAVAL .OR. int(y(i)) .eq. NAVAL) then
			skip = skip + 1 ! count up number of missing values
			cycle ! skip missing values
		end if
		multav = multav + x(i) * y(i) 
	end do

	if(skip .ne. n) then ! Return value
		multav = multav / dble(n-skip)
	else
		multav = NAVAL
	end if
end function multav


