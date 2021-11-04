source("R/construct_model.R")
source("R/get_data.R")

dat <- get_data()
y <- matrix(dat[, 1:2], ncol = 1)

mod <- construct_model(M = 1)
output <- kalman_filter(y, mod)

# dlm --------------------------------------------------------------------------
library(dlm)

y <- dat$AAPL[-1]
mod <- dlmModPoly(order = 1)
dlm_output <- dlmFilter(y, mod)
