# Kalman filter using discount factors for stochastic volatility models --------
filter <- function(y, mod, delta) {
    # y - N x M matrix. N is the # of observations and M is the # of time series
    
    # unpack model
    FF <- mod$FF
    GG <- mod$GG
    V <- mod$V
    m0 <- mod$m0
    C0 <- mod$C0
    N = nrow(y)
    M = ncol(y)
    P = nrow(GG)
    
    # mean and variance-covariance matrix for filtering distribution of
    # state t given all observations through time t
    m <- matrix(0, nrow = P, ncol = N + 1)
    C <- array(0, dim = c(P, P, N + 1))
    
    # mean and variance-covariance matrix for one-step-ahead (Gaussian) predictive
    # distribution of state t given all observations through time t-1
    a <- matrix(0, nrow = P, ncol = N)
    R <- array(0, dim = c(P, P, N))
    
    # mean and variance-covariance matrix for one-step-ahead (Gaussian) predictive
    # distribution of observation t given all observations through time t-1
    f <-  matrix(0, nrow = M, ncol = N)
    Q <- array(0, dim = c(M, M, N))
    
    # Kalman filter
    # Intialize
    m[, 1] <- m0
    C[, , 1] <- C0
    
    for(n in 1:N){
        
        # prior at time t
        a[, n] <- GG %*% matrix(m[, n], nrow = P, ncol = 1)
        R[, , n] <- (GG %*% C[, , n] %*% t(GG))/delta
        
        # one-step forecast
        f[, n] <- FF %*% a[, n]
        Q[, , n] = FF %*% R[, , n] %*% t(FF) + V
        
        # posterior at t
        At <- R[, , n] %*% t(FF) %*% solve(Q[, , n])
        
        m[, n+1] <- a[, n] + At %*% (y[n, ] - f[, n])
        C[, , n+1] <- R[, , n] - At %*% Q[, , n] %*% t(At)

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