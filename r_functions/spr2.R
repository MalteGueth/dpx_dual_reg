# Semi-partion R2
spR2 <- function(model) {
  anov_mod <- as.data.frame(model)

  # remove intercept
  anov_mod <- anov_mod[-1, ]

  rr <- ((anov_mod$Df / anov_mod$Df.res) *
    anov_mod[, grepl(names(anov_mod), pattern = '^F')] ) /
    (1 +
      ((anov_mod$Df / anov_mod$Df.res) *
        anov_mod[, grepl(names(anov_mod), pattern = '^F')] ))

  names(rr) <- row.names(anov_mod)

  print(rr)

}
