library(MASS)

omnibus <- function(Y, E = NULL, M, X = NULL) {
  Y = as.matrix(Y)
  E = as.matrix(E)
  M = as.matrix(M)
  X = as.matrix(X)
  p = ncol(E)
  q = ncol(M)
  r = ncol(X)

  # split samples
  comp = complete.cases(Y, E, M, X)
  incomp = !comp
  n = length(Y)
  mis = sum(incomp) / n

  igwas <- function(Y, E = NULL, M, X = NULL)
  {
    Y = as.matrix(Y)
    E = as.matrix(E)
    M = as.matrix(M)
    X = as.matrix(X)
    p = ncol(E)
    q = ncol(M)
    r = ncol(X)
    # fit IG model
    fit1 = lm(Y ~ E + X)
    ystar = E %*% coef(fit1)[2:(1 + p)]
    fit2 = lm(ystar ~ M + X)
    return(list(
      beta = coef(fit2)[2],
      Tstat = summary(fit2)$coefficients[2, 3],
      pvalue = summary(fit2)$coefficients[2, 4]
    ))

  }

  if (q > 1) {
    fit_gwas = apply(M[comp, ], 2, function(s)
    {
      fit = igwas(Y[comp, ], E[comp, ], s, X[comp, ])

      return(unlist(fit))

    })
    p_ig = fit_gwas[3, ]
    mis_lm = apply(M[incomp, ], 2, function(s) {
      fit = lm(Y[incomp, ] ~ s + X[incomp, ])

      return(unlist(summary(fit)$coefficients[2, 3:4]))

    })
    p_lm = mis_lm[2, ]
  } else if (q == 1) {
    fit_gwas = unlist(igwas(Y[comp, ], E[comp, ], M[comp, ], X[comp, ]))
    p_ig = fit_gwas[3]
    fit = lm(Y[incomp, ] ~ M[incomp, ] + X[incomp, ])

    mis_lm = unlist(summary(fit)$coefficients[2, 3:4])
    p_lm = mis_lm[2]
  }


  # 1. general weight
  weight_gen = ifelse(p_lm < 0.05, sqrt(-log10(p_lm)), 1)
  weight_lin = weight_gen / mean(weight_gen)
  p_gen = p_ig / weight_lin

  # 2. Reverse weight
  weight_gwas_general = ifelse(p_ig < 0.05, sqrt(-log10(p_ig)), 1)
  weight_gwas_general = weight_gwas_general / mean(weight_gwas_general)
  p_rev = p_lm / weight_gwas_general

  # omnibus method
  omn = cbind(p_gen, p_rev)
  p_omn = 0.5 - atan(apply(omn, 1, function(x) {
    (1 - mis) * tan((0.5 - x[1]) * pi) + mis * tan((0.5 - x[2]) * pi)
  })) / pi

  return(list(
    general = p_gen,
    reverse = p_rev,
    omnibus = p_omn
  ))
}


