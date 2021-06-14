#' Conjugate Gradient algorithm
#'
#' The function solves a symmetric positive definite system by constructing conjugate directions.
#' Similar to Lanczos, R is matrix containing Krylov orthogonal vectors (non-normalized), 
#' P is a matrix of conjugate directions, and H is a tridiagonal matrix with Lanczos coefficients.
#' @param A A symmetric positive definite matrix for the system Ax = b. No default value.
#' @param b A right hand side vector in Ax = b. No default value.
#' @param x An initial guess which produces the first residual r = b - Ax. Default is a random vector.
#' @param m Maximum number of iterations to be run. Default is length(b).
#' @param tol The relative tolerance for residual norms. Default is 1e-7.
#' @keywords Conjugate Gradient, Lanczos, Krylov
#' @export
#' @examples
#' A = matrix(rnorm(100^2), 100)
#' A = t(A)%*%A
#' b = rnorm(100)
#' res = cg(A, b, 80)
#' plot(solve(A, b), t="l")
#' lines(res$x, col="red")
cg = function(A, b, x = rnorm(length(b)), m = length(b), tol = 1e-7){
    # initializing vectors and matrices
    n = length(b)
    R = matrix(NA, n, m+1)
    P = matrix(NA, n, m+1)
    alpha = numeric(m)
    beta = numeric(m)
    d = numeric(m)
    e = numeric(m)
    cont = numeric(m)
    # the first residual
    R[,1] = b - drop(A%*%x)
    P[,1] = R[,1]
    for(j in 1:m){
        if(j==1 | j%%50==0){
            cat(paste("CG: iteration #", j, "..\n", sep=""))
        }
        # generating vectors
        w = drop(A%*%P[,j])
        alpha[j] = drop(t(R[,j])%*%R[,j]) / drop(t(w)%*%P[,j])
        x = x + alpha[j] * P[,j]
        R[,j+1] = R[,j] - alpha[j] * w
        # computing entries of tridiagonal Lanczos matrix
        if(j==1){
            d[j] = 1/alpha[j]
        } else {
            d[j] = 1/alpha[j] + beta[j-1]/alpha[j-1]
            e[j] = sqrt(beta[j-1])/alpha[j-1]
        }
        # computing continuant (i.e. determinant)
        if(j == 1){
            cont[j] = d[j]
        } else if(j == 2){
            cont[j] = d[j] * cont[j-1] - e[j]^2
        } else {
            cont[j] = d[j] * cont[j-1] - e[j]^2 * cont[j-2]
        }
        # checking early convergence
        if(drop(t(R[,j+1])%*%R[,j+1]) / drop(t(R[,1])%*%R[,1]) < tol){
            m = j
            cat(paste("CG: converged at iteration #", j, ".\n", sep=""))
            break
        } else if(j==m){
            cat(paste("CG: reached max iteration #", j, ".\n", sep=""))
        }
        beta[j] = drop(t(R[,j+1])%*%R[,j+1]) / drop(t(R[,j])%*%R[,j])
        P[,j+1] = R[,j+1] + beta[j] * P[,j]
    }
    # forming matrices
    R = R[,1:m]
    P = P[,1:m]
    H = tridiag(d[1:m], e[2:m])
    cont = cont[1:m]
    return(list(x = x, R = R, P = P, H = H, det = log(cont)))
}