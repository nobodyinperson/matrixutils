#### Matrix-Funktionen ####
#' Flip matrix vertically
#' 
#' @description Flip given matrix on a horizontal axis, thus returning it upside down.
#' @param x a 2d-matrix
#' @return the given matrix flipped upside down
#' @details \code{apply(x,2,rev)}
#' @examples
#' a<-matrix(1:9,3,3)
#' a
#' #     [,1] [,2] [,3]
#' #[1,]    1    4    7
#' #[2,]    2    5    8
#' #[3,]    3    6    9
#' flip.vert(a)
#' #     [,1] [,2] [,3]
#' #[1,]    3    6    9
#' #[2,]    2    5    8
#' #[3,]    1    4    7
#' @export
flip.vert = function(x) {
	if(length(dim(x))!=2) stop("Argument 'x' has to be a 2d-matrix!")
	apply(x,2,rev)
}

#' Flip matrix horizontally
#' 
#' @description Flip given matrix on a vertical axis, returning it mirrored from left to right.
#' @param x a 2d-matrix
#' @return the given matrix flipped from left to right
#' @details \code{t(apply(x,1,rev))}
#' @examples
#' a<-matrix(1:9,3,3)
#' a
#' #     [,1] [,2] [,3]
#' #[1,]    1    4    7
#' #[2,]    2    5    8
#' #[3,]    3    6    9
#' flip.horz(a)
#' #     [,1] [,2] [,3]
#' #[1,]    7    4    1
#' #[2,]    8    5    2
#' #[3,]    9    6    3
#' @export
flip.horz = function(x) {
	if(length(dim(x))!=2) stop("Argument 'x' has to be a 2d-matrix!")
	t(apply(x,1,rev))
}


#' Rotate matrix 90 degrees clockwise
#' 
#' @param x a 2d-matrix
#' @return the given matrix rotated 90 degrees clockwise
#' @details \code{t(apply(x,2,rev))}
#' @examples
#' a<-matrix(1:9,3,3)
#' a
#' #     [,1] [,2] [,3]
#' #[1,]    1    4    7
#' #[2,]    2    5    8
#' #[3,]    3    6    9
#' rotate.right(a)
#' #     [,1] [,2] [,3]
#' #[1,]    3    2    1
#' #[2,]    6    5    4
#' #[3,]    9    8    7
#' @export
rotate.right = function(x) {
	if(length(dim(x))!=2) stop("Argument 'x' has to be a 2d-matrix!")
	t(apply(x,2,rev))
}

#' Rotate matrix 90 degrees counter-clockwise
#' 
#' @param x a 2d-matrix
#' @return the given matrix rotated 90 degrees counter-clockwise
#' @details \code{t(apply(x,2,rev))}
#' @examples
#' a<-matrix(1:9,3,3)
#' a
#' #     [,1] [,2] [,3]
#' #[1,]    1    4    7
#' #[2,]    2    5    8
#' #[3,]    3    6    9
#' rotate.left(a)
#' #     [,1] [,2] [,3]
#' #[1,]    7    8    9
#' #[2,]    4    5    6
#' #[3,]    1    2    3
#' @export
rotate.left  = function(x) {
	if(length(dim(x))!=2) stop("Argument 'x' has to be a 2d-matrix!")
	apply(x,1,rev)
}

