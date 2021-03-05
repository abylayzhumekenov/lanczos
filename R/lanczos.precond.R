#' Preconditioned Lanczos algorithm
#'
#' The function solve a symmetric positive definite system by split preconditioning the system Ax = b.
#' Then it is handled by a plain Lanczos algorithm.
#' @param A A symmetric positive definite matrix for the system Ax = b. No default value.
#' @param b A right hand side vector in Ax = b. No default value.
#' @param m Maximum number of iterations to be run. Default is length(b).
#' @param tol The relative tolerance for residual norms. Default is 1e-7.
#' @param diag.comp A type of diagonal compensation for IC(0) preconditioning. Default is "aj".
#' @param alpha The scale of diagonal compensation. Default is 1.
#' @param L The preconditioner to supply. If missing, the IC(0) is above parameters is computed.
#' @keywords Lanczos, Krylov, tridiagonalization
#' @export
#' @examples
#' A = matrix(rnorm(100^2), 100)
#' A = t(A)%*%A
#' b = rnorm(100)
#' res = planczos(A, b)
#' plot(solve(A, b), t="l")
#' lines(res$x, col="red")
planczos = function(A, b, m = length(b), tol = 1e-7, diag.comp = "aj", alpha = 1, L){
    # preconditioning
    cat(paste("Preconditioning: started...\n"))
    time.0 = as.numeric(Sys.time())
    if(missing(L)){
        L = ichol(A, diag.comp, alpha)
    } else {
        cat(paste("Preconditioning: using supplied preconditioner.\n"))
    }
    Q = solve(L, A)%*%solve(t(L))
    c = solve(L, b)
    time.1 = as.numeric(Sys.time())
    time.precond = time.1 - time.0
    cat(paste("Preconditioning: finished (", round(time.precond, 5), " sec).\n", sep=""))

    # solving
    cat(paste("Lanczos: started...\n"))
    out.lanc = lanczos(Q, c, m, tol)
    x = drop(solve(t(L), out.lanc$x))
    time.2 = as.numeric(Sys.time())
    time.solve = time.2 - time.1
    cat(paste("Lanczos: finished (", round(time.solve, 5), " sec).\n", sep=""))

    time.total = time.2 - time.0
    cat(paste("Total run time: ", round(time.total, 5), " sec.\n\n", sep=""))
    return(list(x = x, L = L, V = out.lanc$V, H = out.lanc$H, det = out.lanc$det,
                time.precond = time.precond, time.solve = time.solve))
}
