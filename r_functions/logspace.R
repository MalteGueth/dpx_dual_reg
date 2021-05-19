logspace <- function(d1, d2, n) {
  
  exp(log(10)*seq(d1, d2, length.out=n))
  
}

nextpow2 <- function(x) {
  if (!(is.numeric(x) || is.complex(x))) {
    stop(sprintf("argument %s must be numeric or complex",
                 sQuote('x')))
  }
  
  if (length(x) == 0) {
    return(numeric(0))
  }
  
  x[x == 0] <- 1
  return(ceiling(log2(abs(x))))
}