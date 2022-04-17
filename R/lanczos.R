#' Lanczos algorithm
#'
#' The function solves a symmetric positive definite system by factoring the matrix into A=VHV^T form.
#' Here, V is matrix containing Krylov orthogonal vectors, and H is a tridiagonal matrix with Lanczos
#' coefficients.
#' @param A A symmetric positive definite matrix for the system Ax = b. No default value.
#' @param b A right hand side vector in Ax = b. No default value.
#' @param x An initial guess which produces the first residual r = b - Ax. Default is a random vector.
#' @param m Maximum number of iterations to be run. Default is length(b).
#' @param tol The relative tolerance for residual norms. Default is 1e-7.
#' @keywords Lanczos, Krylov, tridiagonalization
#' @export
#' @examples
#' A = matrix(rnorm(100^2), 100)
#' A = t(A)%*%A
#' b = rnorm(100)
#' res = lanczos(A, b, 80)
#' plot(solve(A, b), t="l")
#' lines(res$x, col="red")
lanczos = function(A, b, x = rnorm(length(b)), m = length(b), tol = 1e-7){
    # initializing vectors and matrices
    n = length(b)
    V = matrix(NA, n, m+1)
    alpha = numeric(m)
    beta = numeric(m+1)
    det = numeric(m)
    # cont = numeric(m)
    # the first residual
    w = b - A%*%x
    beta[1] = norm(w, "2")
    V[,1] = w/beta[1]
    for(j in 1:m){
        if(j==1 | j%%50==0){
            cat(paste("Lanczos: iteration #", j, "..\n", sep=""))
        }
        # three term reccurence relation
        w = A%*%V[,j]
        alpha[j] = drop(t(V[,j])%*%w)
        # # Iterated MGS is costly and slow
        # for(k in 1:j){
        #     w = w - V[,k] * drop(t(V[,k])%*%w)
        # }
        # for(k in 1:j){
        #     w = w - V[,k] * drop(t(V[,k])%*%w)
        # }
        # # Use iterated CGS instead, cheaper but of the same quality
        w = w - V[,1:j] %*% (t(V[,1:j])%*%w)
        w = w - V[,1:j] %*% (t(V[,1:j])%*%w)
        beta[j+1] = norm(w, "2")
        # computing continuant (i.e. determinant)
        det[j] = determinant(as.matrix(tridiag(alpha[1:j], beta[2:j])))$modulus[1]
        # # Overflows...
        # if(j == 1){
        #     cont[j] = alpha[j]
        # } else if(j == 2){
        #     cont[j] = alpha[j] * cont[j-1] - beta[j]^2
        # } else {
        #     cont[j] = alpha[j] * cont[j-1] - beta[j]^2 * cont[j-2]
        # }
        # check for early convergence or append Krylov vector
        if(beta[j+1]/beta[1] < tol){
            cat(paste("Lanczos: converged at iteration #", j, ".\n", sep=""))
            m = j
            break
        } else if(j==m){
            cat(paste("Lanczos: reached max iteration #", j, ".\n", sep=""))
        } else {
            V[,(j+1)] = w/beta[j+1]
        }
    }
    # constructing smaller system and solving it
    cat(paste("Lanczos: back solving...\n"))
    det = det[1:m]
    # cont = cont[1:m]
    H = tridiag(alpha[1:m], beta[2:m])
    V = V[,1:m]
    y = as.matrix(solve(H, c(beta[1], numeric(m-1)), tol = 1e-16))
    x = x + drop(V%*%y)
    return(list(x = x, V = V, H = H, det = det))
}
