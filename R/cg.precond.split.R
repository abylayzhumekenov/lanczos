#' Split Preconditioned Conjugate Gradient algorithm
#'
#' The function solve a symmetric positive definite system by split preconditioning the system Ax = b.
#' Then it is handled by a plain CG algorithm.
#' @param A A symmetric positive definite matrix for the system Ax = b. No default value.
#' @param b A right hand side vector in Ax = b. No default value.
#' @param x An initial guess which is produces the first residual r. Default is a random vector.
#' @param m Maximum number of iterations to be run. Default is length(b).
#' @param tol The relative tolerance for residual norms. Default is 1e-7.
#' @param diag.comp A type of diagonal compensation for IC(0) preconditioning. Default is "aj".
#' @param alpha The scale of diagonal compensation. Default is 1.
#' @param L The preconditioner to supply. If missing, the IC(0) with parameters above is computed.
#' @keywords Conjugate Gradient, Lanczos, Krylov
#' @export
#' @examples
#' A = matrix(rnorm(100^2), 100)
#' A = t(A)%*%A
#' b = rnorm(100)
#' res = pcgs(A, b)
#' plot(solve(A, b), t="l")
#' lines(res$x, col="red")
pcgs = function(A, b, x = rnorm(length(b)), m = length(b), tol = 1e-7, diag.comp = "aj", alpha = 1, L){
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
    x = t(L)%*%x
    time.1 = as.numeric(Sys.time())
    time.precond = time.1 - time.0
    cat(paste("Preconditioning: finished (", round(time.precond, 5), " sec).\n", sep=""))
    
    # solving
    cat(paste("CG: started...\n"))
    out.cg = cg(Q, c, x, m, tol)
    x = drop(solve(t(L), out.cg$x))
    time.2 = as.numeric(Sys.time())
    time.solve = time.2 - time.1
    cat(paste("CG: finished (", round(time.solve, 5), " sec).\n", sep=""))
    
    time.total = time.2 - time.0
    cat(paste("Total run time: ", round(time.total, 5), " sec.\n\n", sep=""))
    return(list(x = x, L = L, R = out.cg$R, P = out.cg$P, H = out.cg$H, det = out.cg$det,
                time.precond = time.precond, time.solve = time.solve))
}
