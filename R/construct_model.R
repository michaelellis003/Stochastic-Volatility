
construct_model <- function(M) {
    
    mod <- list()
    mod$FF <- matrix(rep(1, M), ncol = 1)
    mod$GG <- diag(M)
    mod$V <- diag(M)
    mod$W <- diag(M)
    mod$m0 <- rep(0, M)
    mod$C0 <- 10^7 * diag(M)
    mod$n0 <- M+1
    mod$D0 <- 10^7 * diag(M)
    
    return(mod)
}