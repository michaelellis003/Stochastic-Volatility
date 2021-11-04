
inverse_smw <- function(A, U, C, V) {
    
    A_inv <- diag(1/diag(A))
    C_inv <- solve(C)
    
    return(
        A_inv - A_inv %*% U %*% solve(C_inv + V %*% A_inv %*% U) %*% V %*% A_inv
    )
    
}