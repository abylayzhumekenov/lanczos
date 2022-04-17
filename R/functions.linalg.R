#' Condition number
#'
#' Computes the condition number of the given matrix.
#' @param A A non-singular matrix. No default value.
#' @keywords condition number
#' @export
#' @examples
#' A = diag(c(1,2,3))
#' k = cond(A)
cond = function(A){
    eig = eigen(A)$values
    return(max(abs(eig))/min(abs(eig)))
}

#' Tridiagonal matrix
#'
#' Constructs a tridiagonal matrix from three vectors supplied. The matrix size is adjusted according
#' to the first vector (n). The function fails, if the other two dimensions are less than n-1.
#' However, if the dimensions are larger, the function will truncate the vectors b and c.
#' @param a A vector of size n, whose entries will form a diagonal. No default value.
#' @param b A vector of size n-1, whose entries will form a sub-diagonal. No default value.
#' @param c A vector of size n-1, whose entries will form a super-diagonal. Defaults is c = b.
#' @keywords tridiagonal
#' @export
#' @examples
#' tridiag(c(1,2,3), c(4,5,6))
#' tridiag(c(1,2,3), c(4,5,6), c(7,8,9))
tridiag = function(a, b, c = b){
    if(length(a)==1 & (missing(b) | missing(c))){
        return(matrix(a))
    } else if(length(a)==2){
        return(matrix(c(a[1],b[1],c[1],a[2]), 2))
    } else {
        H = diag(a)
        diag(H[-1,]) = b[1:(length(a)-1)]
        diag(H[,-1]) = c[1:(length(a)-1)]
        return(H)
    }
}

#' Determinant of a tridiagonal matrix
#'
#' Computes the determinant of a tridiagonal matrix using continuants three term relation.
#' @param H A tridiagonal matrix whose determinant is to be computed. No default value.
#' @keywords determinant
#' @export
#' @examples
#' H = tridiag(c(1,2,3), c(4,5,6), c(7,8,9))
#' det.tridiag(H)
det.tridiag = function(H){
    m = dim(H)[1]
    val = numeric(m)
    val[1] = H[1,1]
    val[2] = H[2,2]*val[1] - H[2,1]^2
    for(i in 3:m){
        val[i] = H[i,i]*val[i-1] - H[i,i-1]^2*val[i-2]
    }
    return(val[m])
}
