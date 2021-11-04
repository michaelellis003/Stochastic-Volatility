# Clean and save raw data ------------------------------------------------------
get_data <- function() {
    # library(tidyquant)
    # 
    # ## get all symbols from the S&P 500
    # symbols <- tq_index("SP500")
    # symbols <- symbols$symbol[1:10]
    # symbols[which(symbols == "BRB.B")] <- "BRB-B"
    # 
    # prices <- tq_get(symbols,
    #                  from = "2017-01-01",
    #                  to = "2017-03-01",
    #                  get = "stock.prices")
    # 
    # unique(prices$symbol)
    
    dat <- read.csv(file = "data/logreturns.csv")
    
    return(dat)
}