# DLM Kalman Filter ------------------------------------------------------------
kalman_filter <- function(y, mod) {
    # y - N x M matrix. N is the # of observations and M is the # of time series
    
    # unpack model
    FF <- mod$FF
    GG <- mod$GG
    V <- mod$V
    W <- mod$W
    m0 <- mod$m0
    C0 <- mod$C0
    N = nrow(y)
    M = ncol(y)
    P = nrow(GG)
    
    # mean and variance-covariance matix for filtering distribution of
    #state t given all observations through time t
    m <- matrix(0, nrow = P, ncol = N + 1)
    C <- array(0, dim = c(P, P, N + 1))
    U.C <- matrix(0, nrow = P, ncol = P)
    D.C <- matrix(0, nrow = P, ncol = P)
    
    # mean and variance-covariance matix for one-step-ahead (Gaussian) predictive
    # distribution of state t given all observations through time t-1
    a <- matrix(0, nrow = P, ncol = N)
    R <- array(0, dim = c(P, P, N))
    U.R <- matrix(0, nrow = P, ncol = P)
    D.R <- matrix(0, nrow = P, ncol = P)
    
    # mean and variance-covariance matix for one-step-ahead (Gaussian) predictive
    # distribution of observation t given all observations through time t-1
    f <-  matrix(0, nrow = M, ncol = N)
    Q <- array(0, dim = c(M, M, N))
    
    # get sqrt of V
    Vsvd <- La.svd(V, nu = 0)
    U.V <- t(Vsvd$vt)
    D.V <- sqrt(Vsvd$d)
    Dv.inv <- 1/D.V
    Dv.inv[abs(Dv.inv) == Inf] <- 0
    sqrtVinv <- Dv.inv * t(U.V)
    sqrtV <- D.V * U.V
    
    # get sqrt of W
    Wsvd <- La.svd(W , nu = 0)
    sqrtW <- sqrt(Wsvd$d) * Wsvd$vt
    
    # Kalman filter
    # Intialize
    m[, 1] <- m0
    C[, , 1] <- C0
    C_svd <- La.svd(C0)
    D.C <- diag(C_svd$d, nrow = P, ncol = P)
    U.C = t(C_svd$vt)
    
    for(n in 1:N){
        
        a[, n] <- GG %*% matrix(m[, n], nrow = P, ncol = 1)
        MR <- rbind(sqrt(D.C) %*% t(U.C) %*% t(GG),
                    sqrtW)
        R[, , n] <- t(MR) %*% MR
        MRsvd <- La.svd(MR)
        U.R <- t(MRsvd$vt)
        D.R <- t(diag(MRsvd$d, P)) %*% diag(MRsvd$d, P)
        D.Rinv <- 1/D.R
        D.Rinv[abs(D.Rinv) == Inf] <- 0
        
        f[, n] = FF %*% a[, n]
        Q[, , n] = FF %*% R[, , n] %*% t(FF) + V
        
        # No observations in a row are missing
        Qinv = solve(Q[, , n])
        m[, n+1] = a[, n] + R[, , n] %*% t(FF) %*% Qinv %*% (y[n] - f[, n]);
        
        MC <- rbind(sqrtVinv %*% FF %*% U.R,
                    sqrt(D.Rinv))
        MCsvd <- La.svd(MC)
        D.MCinv = diag(1/MCsvd$d, M+P, P)
        
        U.C = U.R %*% t(MCsvd$vt)
        D.C = t(D.MCinv) %*% D.MCinv
        C[, , n+1] = U.C %*% D.C %*% t(U.C)
    }
    
    return(list(
        a = a,
        R = R,
        f = f,
        Q = Q,
        m = m,
        C = C
    ))
    
}