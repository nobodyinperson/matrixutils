#### Advanced matrix functions ####

#' Remaining bounds after shifting two matrices against each other
#' @description Determine the indices that limit the remainig 
#'  matrix after shifting two matrices against each other.
#' @param dimshift Dimension of shifted matrix, as returned by \code{\link{dim}}.
#' @param dimbase Dimension of static matrix, as returned by \code{\link{dim}}.
#' @param rowshift integer. column shift. See details.
#' @param colshift integer. row shift. See details.
#' @details
#'  Imagine two matrices \code{mshift} and \code{mbase}. You shift matrix \code{mshift} against \code{mbase} by 
#'  a certain amount of rows and columns. Now you have a remaining 
#'  matrix that is possibly smaller than \code{mshift} or \code{mbase}. The bounds of this 
#'  remaining matrix are in both matrices.\cr\cr
#'  This function returns only the bounds of the matrices \code{mshift} and \code{mbase}
#'  that limit the remainig matrix.
#'  \cr\cr
#'  The \code{rowshift} and \code{colshift} parameters describe the shifts in units of rows and columns respectively. 
#'  A positive value indicates an \strong{downwards} (!) / rightwards shift of matrix \code{m1} against \code{m2}.
#'  A negative value indicates a \strong{upwards} (!) / leftwards shift. Values greater than 
#'  the respective dimensions \code{dim} lead to indefinite bounds because there 
#'  is no remaining matrix after shifting. In this case, all bounds are set to \code{NA}.
#'  \cr\cr
#' @return 
#'  vector of bounds in the following order:
#'  \itemize{
#'  \item{\code{mirosh}: lower rowbound of \code{mshift}}
#'  \item{\code{marosh}: upper rowbound of \code{mshift}}
#'  \item{\code{micosh}: lower columnbound of \code{mshift}}
#'  \item{\code{macosh}: upper columnbound of \code{mshift}}
#'  \item{\code{miroba}: lower rowbound of \code{mbase}}
#'  \item{\code{maroba}: upper rowbound of \code{mbase}}
#'  \item{\code{micoba}: lower columnbound of \code{mbase}}
#'  \item{\code{macoba}: upper columnbound of \code{mbase}}
#'  }
#' @examples
#' # Two equal 3x3-matrices
#' (m1 <- m2 <- matrix(1:9,3,3))
#' #      [,1] [,2] [,3]
#' # [1,]    1    4    7
#' # [2,]    2    5    8
#' # [3,]    3    6    9
#'  
#' # Shift by 0 rows and 0 columns
#' (bounds = matrix.shift.bounds(dim(m1),dim(m2),0,0))
#' # [1] 1 3 1 3 1 3 1 3
#' # Remaining part of m1:
#'  m1[bounds[1]:bounds[2],bounds[3]:bounds[4]]
#' #      [,1] [,2] [,3]
#' # [1,]    1    4    7
#' # [2,]    2    5    8
#' # [3,]    3    6    9
#' # Everything! Because nothing was shifted!

#' @export
matrix.shift.bounds = function(dimbase,dimshift,rowshift,colshift) {
	NROWSHIFT = dimshift[1]
	NROWBASE  = dimbase[1]
	NCOLSHIFT = dimshift[2]
	NCOLBASE  = dimbase[2]
	ROWDIFF = (NROWBASE - NROWSHIFT) / 2
	COLDIFF = (NCOLBASE - NCOLSHIFT) / 2
	MAXROW = (NROWBASE + NROWSHIFT) / 2
	MAXCOL = (NCOLBASE + NCOLSHIFT) / 2
	
	# Gucken, ob Verschiebung zu groß ist
	if(rowshift <= -floor(MAXROW) | rowshift >= ceiling(MAXROW) | colshift <= -floor(MAXCOL) | colshift >= ceiling(MAXCOL) )
		return(rep(NA,8)) # Nur NA zurückgeben
	
	c(max(         1,         1 - ( rowshift + floor   (ROWDIFF) ) ),
		min( NROWSHIFT, NROWSHIFT - ( rowshift - ceiling (ROWDIFF) ) ),
		max(         1,         1 - ( colshift + floor   (COLDIFF) ) ),
		min( NCOLSHIFT, NCOLSHIFT - ( colshift - ceiling (COLDIFF) ) ),
		max(         1,         1 + ( rowshift + floor   (ROWDIFF) ) ),
		min(  NROWBASE,  NROWBASE + ( rowshift - ceiling (ROWDIFF) ) ),
		max(         1,         1 + ( colshift + floor   (COLDIFF) ) ),
		min(  NCOLBASE,  NCOLBASE + ( colshift - ceiling (COLDIFF) ) )
	)
}
# Kompilieren, wenn möglich
if(requireNamespace("compiler",quietly=T))
	matrix.shift.bounds <- compiler::cmpfun(matrix.shift.bounds)



#' Shift two matrices against each other and apply function to the overlap
#' @description Shift two matrices against each other and apply a function to the overlap
#' @param mshift,mbase matrices of to be shifted against each other. \code{NA}s are allowed, see argument \code{na.value}.
#' @param func character. One of "rmse" or "multav". See details.
#' @param keepbasedim logical. Only shift \code{mshift} as much that the center stays inside \code{mbase}.
#' 	This results in an output matrix of equal shape like \code{mbase}. Handy for sliding means. Defaults to FALSE.
#' @param radial logical. Only shift \code{mshift} up to a circle-like shape against \code{mbase}. 
#' 	That is, The center of \code{mshift} won't exceed the biggest ellipsis fitting into the center of \code{mbase}. 
#' 	This results in an output matrix of equal shape like \code{mbase} with \code{NA}s in the edges. Handy for circular images.
#' @param na.value integer. The value used internally in FORTRAN for \code{NA}s. Make sure this value does not appear in neither \code{mbase} nor \code{mshift}.
#' 	Defaults to \code{-99999}.
#' @details
#'  Matrix \code{mshift} is shifted against \code{mbase} and the overlapping matrices passed to \code{fun}.
#'  The return value of selected function \code{func} is put into the resulting matrix at the position of the shift.
#'  
#'  This function uses a FORTRAN implementation for efficiency.
#'  
#'  The \code{func} argument specifies which type of function is to be 
#'  applied to the overlap. It can be:
#'  \itemize{
#'  	\item \code{rmse}: The Root Mean Square Error RMSE is calculated on the overlap.
#'  				Handy for finding the shift of the maximum resemblance of \code{mbase} and \code{mshift}.
#'  	\item \code{multav}: The overlapping matrix parts are multiplicated elementwise and then averaged.
#'  	      Handy for sliding means.
#'  }
#' @return a matrix filled with the results of \code{func} and row- and colnames that indicate the shift of \code{mshift} against \code{mbase}.
#' @useDynLib matrixutils
#' @export
matrices.shift.apply.function <- function(mbase,mshift,func,keepbasedim=FALSE,radial=FALSE,na.value=-99999) {
	if(!is.matrix(mbase) | !is.matrix(mshift))
		stop("arguments 'mshift' and 'mbase' have to be matrices!")
	
	funcavail <- c("rmse","multav") # available functions in the order of the FORTRAN code
	funcnr <- pmatch(func,funcavail)
	if(is.na(funcnr))
		stop("argument 'func' must be one of \"",paste(funcavail,collapse='","'),'"!',sep="")
	
	# prepare variables for FORTRAN code
	if(keepbasedim) {
		mres      <- array(0,dim = dim(mbase))
		colshifts <- -floor((ncol(mbase)-1)/2):ceiling((ncol(mbase)-1)/2) 
		rowshifts <- -floor((nrow(mbase)-1)/2):ceiling((nrow(mbase)-1)/2)
	} else {
		mres      <- array(0,dim = dim(mbase) + dim(mshift) - 1)
		rowshifts <- -(floor((nrow(mshift)+nrow(mbase))/2)-1):(ceiling((nrow(mshift)+nrow(mbase))/2)-1)
		colshifts <- -(floor((ncol(mshift)+ncol(mbase))/2)-1):(ceiling((ncol(mshift)+ncol(mbase))/2)-1)
	}
	
	na.value <- -99999 # value for NAs 
	mbase  [is.na(mbase)] <- na.value # Replace NAs for FORTRAN
	mshift[is.na(mshift)] <- na.value # Replace NAs for FORTRAN
	
	# Call FORTRAN
	res <- .Fortran("matshiftapplyfunc"
									,mbase     = as.double  (mbase)
									,mbrow     = as.integer (nrow(mbase))
									,mbcol     = as.integer (ncol(mbase))
									,mshift    = as.double  (mshift)
									,msrow     = as.integer (nrow(mshift))
									,mscol     = as.integer (ncol(mshift))
									,mres      = as.double  (mres)
									,mrrow     = as.integer (nrow(mres))
									,mrcol     = as.integer (ncol(mres))
									,rowshifts = as.integer (rowshifts)
									,colshifts = as.integer (colshifts)
									,func      = as.integer (funcnr)
									,radial    = as.logical (radial)
									,na        = as.integer (na.value)
	)
	
	# Reshape to array and add row and colnames
	res <- array(res$mres, dim = dim(mres))
	res[res==na.value] <- NA # Re-replace missing values
	row.names(res) <- rowshifts
	colnames(res)  <- colshifts
	
	return(res)
}
