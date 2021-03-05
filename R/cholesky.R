#' A complete Cholesky factorization
#'
#' The function computes Cholesky factor as a lower triangular matrix.
#' @param A A symmetric positive definite matrix to be factored. No default value.
#' @keywords Cholesky
#' @export
#' @examples
#' A = matrix(c(1,-1,-1,1), ncol=2)
#' L = chol(A)
chol = function(A){
    return(t(chol.default(A)))
}

#' An Incomplete Cholesky factorization
#'
#' The function compues IC(0) factorization as a lower triangular matrix.
#' @param A A symmetric positive definite matrix to be factored. No default value.
#' @param diag.comp A type of diagonal compensation. Either "aj" (default) for Ajiz-Jennings,
#' "as" for absolute shift (diag(A) = diag(A) + alpha), "rs" for relative shift
#' (diag(A) = diag(A) + alpha * diag(A)).
#' @param alpha Value for diagonal compensation. Default is 1. For matrices not diagonally dominant,
#' small values of the parameter could lead to negative pivots. When such are encountered,
#' the pivot is replaced by 1.
#' @keywords Cholesky, Incomplete Cholesky, Incomplete LU
#' @export
#' @examples
#' A = matrix(c(1,-1,-1,1), ncol=2)
#' L = ichol(A, diag.comp="as", alpha=1)
ichol = function(A, diag.comp = "aj", alpha = 1){
    n = dim(A)[1]
    L = matrix(0, n, n)
    if(diag.comp=="rs"){
        diag(A) = diag(A) + alpha * diag(A) # relative shift
    } else if(diag.comp=="as"){
        diag(A) = diag(A) + alpha # absolute shift
    }
    for(i in 1:n){
        for(j in i:n){
            if(A[i,j]!=0){
                A[i,j] = A[i,j] - sum(L[1:i,i]*L[1:i,j])
                # Above, it is actually L[1:(i-1),j], but L[1:i,j] is zero anyway at this point.
            } else if(diag.comp=="aj"){
                c = abs(sum(L[1:i,i]*L[1:i,j]))
                a = sqrt(A[i,i]/A[j,j])
                A[i,i] = A[i,i] + alpha * c * a # Ajiz-Jennings
                A[j,j] = A[j,j] + alpha * c / a # compensation
            }
        }
        if(A[i,i]>0){
            L[i,i] = sqrt(A[i,i])
        } else {
            L[i,i] = 1
            cat(paste("Incomplete Cholesky: non-positive pivot at [", i, ",", i, "].\n", sep=""))
            cat(paste("Incomplete Cholesky: consider increasing the alpha parameter.\n"))
        }
        if(i<n){
            L[i,(i+1):n] = A[i,(i+1):n] / L[i,i]
        }
    }
    return(t(L))
}
