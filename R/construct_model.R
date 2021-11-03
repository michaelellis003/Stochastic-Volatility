
construct_model <- function(M) {
    
    mod <- list()
    mod$FF <- diag(M)
    mod$GG <- diag(M)
    
    return(mod)
}