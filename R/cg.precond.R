#' Preconditioned Conjugate Gradient algorithm
#'
#' The function solve a symmetric positive definite system by left preconditioning the system Ax = b.
#' The Preconditoined CG algorithm is implemented using the (preconditioner) M inner product.
#' @param A A symmetric positive definite matrix for the system Ax = b. No default value.
#' @param b A right hand side vector in Ax = b. No default value.
#' @param x An initial guess which is produces the first residual r. Default is a random vector.
#' @param m Maximum number of iterations to be run. Default is length(b).
#' @param tol The relative tolerance for residual norms. Default is 1e-7.
#' @param diag.comp A type of diagonal compensation for IC(0) preconditioning. Default is "aj".
#' @param alpha The scale of diagonal compensation. Default is 1.
#' @param L The preconditioner to supply. If missing, the IC(0) is above parameters is computed.
#' @keywords Conjugate Gradient, Lanczos, Krylov
#' @export
#' @examples
#' A = matrix(rnorm(100^2), 100)
#' A = t(A)%*%A
#' b = rnorm(100)
#' res = pcg(A, b)
#' plot(solve(A, b), t="l")
#' lines(res$x, col="red")
pcg = function(A, b, x = rnorm(length(b)), m = length(b), tol = 1e-7, diag.comp = "aj", alpha = 1, L){
    # preconditioning
    cat(paste("Preconditioning: started...\n"))
    time.0 = as.numeric(Sys.time())
    if(missing(L)){
        L = ichol(A, diag.comp, alpha)
    } else {
        cat(paste("Preconditioning: using supplied preconditioner.\n"))
    }
    time.1 = as.numeric(Sys.time())
    time.precond = time.1 - time.0
    cat(paste("Preconditioning: finished (", round(time.precond, 5), " sec).\n", sep=""))
    
    cat(paste("CG: started...\n"))
    # initializing vectors and matrices
    n = length(b)
    R = matrix(NA, n, m+1)
    Z = matrix(NA, n, m+1)
    P = matrix(NA, n, m+1)
    alpha = numeric(m)
    beta = numeric(m)
    d = numeric(m)
    e = numeric(m)
    cont = numeric(m)
    # the first residual
    R[,1] = b - drop(A%*%x)
    Z[,1] = solve(L, R[,1])
    Z[,1] = solve(t(L), Z[,1])
    P[,1] = Z[,1]
    for(j in 1:m){
        if(j==1 | j%%50==0){
            cat(paste("CG: iteration #", j, "..\n", sep=""))
        }
        # generating vectors
        w = drop(A%*%P[,j])
        alpha[j] = drop(t(R[,j])%*%Z[,j]) / drop(t(w)%*%P[,j])
        x = x + alpha[j] * P[,j]
        R[,j+1] = R[,j] - alpha[j] * w
        Z[,j+1] = solve(L, R[,j+1])
        Z[,j+1] = solve(t(L), Z[,j+1])
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
        if(drop(t(R[,j+1])%*%Z[,j+1]) / drop(t(R[,1])%*%Z[,1]) < tol){
            m = j
            cat(paste("CG: converged at iteration #", j, ".\n", sep=""))
            break
        } else if(j==m){
            cat(paste("CG: reached max iteration #", j, ".\n", sep=""))
        }
        beta[j] = drop(t(R[,j+1])%*%Z[,j+1]) / drop(t(R[,j])%*%Z[,j])
        P[,j+1] = Z[,j+1] + beta[j] * P[,j]
    }
    # forming matrices
    R = R[,1:m]
    Z = Z[,1:m]
    P = P[,1:m]
    H = tridiag(d[1:m], e[2:m])
    cont = cont[1:m]
    time.2 = as.numeric(Sys.time())
    time.solve = time.2 - time.1
    cat(paste("CG: finished (", round(time.solve, 5), " sec).\n", sep=""))
    
    time.total = time.2 - time.0
    cat(paste("Total run time: ", round(time.total, 5), " sec.\n\n", sep=""))
    return(list(x = x, L = L, R = R, Z = Z, P = P, H = H, det = log(cont),
                time.precond = time.precond, time.solve = time.solve))
}