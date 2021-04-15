# Thompson, J. Bootstrapping a Confidence Interval for a Difference in Means.
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.579.1330&rep=rep1&type=pdf
bootdif <- function(y, g) {
  # Ensure treatment group is a factor
  g <- as.factor(g)

  # Use the smean.cl.boot function to bootstrap means for
  # variable y for each treatment group (A and B); this code
  # uses 5000 samples, but can easily be changed
  require(Hmisc)
  a <- attr(smean.cl.boot(y[g == levels(g)[1]], B = 5000, reps=TRUE),
            'reps')
  b <- attr(smean.cl.boot(y[g == levels(g)[2]], B = 5000, reps=TRUE),
            'reps')

  # Calculate the observed mean difference between groups
  meandif <- diff(tapply(y, g, mean, na.rm = TRUE))

  # Calculate the 2.5 and 97.5 percentiles of the differences
  # in bootstrapped means
  # (can easily be changed for 90% CI, 99% CI, etc.)
  a.b <- quantile(b - a, c(.025, .975))

  # Prepare object to return
  res <- c(meandif, a.b)
  names(res) <- c('Mean', 'Lower', 'Upper')

  res
}
