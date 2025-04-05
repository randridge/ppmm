#' Calculate Bayes estimates for continuous outcomes
#'
#' @param data_YZ_s list for selected sample (Microdata or Summary statistics). See details below for more details.
#' @param data_Z_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param phi value of phi, fixed to a value between 0 and 1 (defaults to 0)
#' @param phi_character draw phi from specific distribution (i.e."`runif(1)`", defaults to `NULL`)
#' * If `phi_character=NULL`, fix at value set to `phi`;
#' * If `phi_character=` specific distribution (i.e. `runif(1)` or `rbeta(1,1,1)`), then the value of `phi` input to the function is ignored
#' @param ndraws number of draws
#'
#' @returns A `ndraws*10` matrix of Bayes estimates with each row containing:
#' \itemize{
#' \item `phi` : drawn phi in this draw
#' \item `muY_s` : drawn mean of outcome variable for selected sample
#' \item `muY_ns` : drawn mean of outcome variable for non-selected sample
#' \item `sigmaYY_s` : drawn variance of outcome variable for selected sample
#' \item `sigmaYY_ns` : drawn variance of outcome variable for non-selected sample
#' \item `rho_XY_s` : drawn correlation between the proxy and the outcome variable for selected sample
#' \item `muY` : drawn grand mean of outcome variable
#' \item `smub` : drawn SMUB with phi = phi
#' \item `smub0` : drawn SMUB with phi = 0
#' \item `smab` : drawn SMAB
#' }
#' @import stats
#' @import MCMCpack
#' @import mnormt
#' @import MASS
#' @import mvtnorm
#' @importFrom stringr str_sub
#' @export
#' @details
#' \describe{
#'  \item{`data_YZ_s=list(Y,Z)`(Microdata)}{
#'    \itemize{
#'      \item `Y` : `n*1` matrix of outcome variable `Y`
#'      \item `Z` : `n*nZparams` matrix of auxiliary variables `Z`
#'    }
#'  }
#'  \item{`data_YZ_s=list(mean_YZ,var_YZ,n_YZ)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_YZ` : `(1+nZparams)*1` vector of mean (in order) of outcome variable `Y` and auxiliary variables `Z`
#'      \item `var_YZ` : `(1+nZparams)*(1+nZparams)` matrix of variance (in order) of outcome variable `Y` and auxiliary variables `Z`
#'      \item `n_YZ` : number of selected sample
#'    }
#'  }
#'  \item{`data_Z_ns=list(Z)`(Microdata)}{
#'    \itemize{
#'      \item `Z` : `(N-n)*nZparams` matrix of auxiliary variables `Z`
#'    }
#'  }
#'  \item{`data_Z_ns=list(mean_Z, var_Z, n_Z)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_Z` : `nZparams*1` vector of mean of auxiliary variables `Z`
#'      \item `var_Z` : `nZparams*nZparams` matrix of variance of auxiliary variables `Z`
#'      \item `n_Z` : number of non-selected sample
#'    }
#'  }
#' }
#'
#' @seealso [means_mle()]
#' @seealso [prop_bayes()]
#' @seealso [prop_mle()]
#'
#' @examples
#' set.seed(123)
#' # num of total sample
#' N <- 10000
#' # generate auxiliary variables Z
#' Z <- rnorm(N, 0, 1)
#' # define the biserial correlation rho between outcome variable Y and the proxy X
#' rho <- 0.5
#' # define the mean of the outcome variable Y
#' mu_y <- 10
#' a1 <- rho / sqrt(1 - rho^2)
#' # generate the observed variable Y
#' Y <- rnorm(Z, 10 + a1 * Z, 1)
#' # generate the sample selection indicator S with 0.05 of the population selected
#' prob <- plogis(-4 + 0.3 * Z + 0.1 * Y)
#' S <- rbinom(prob, 1, prob)
#' index_s <- S == 1
#' index_ns <- S == 0
#' # generate the list of microdata for selected (data_YZ_s) and non-selected sample (data_Z_ns)
#' data_YZ_s <- list(Y = Y[index_s], Z = Z[index_s])
#' data_Z_ns <- list(Z = Z[index_ns])
#'
#' # draw bayesian result from ppmm model
#' result <- means_bayes(data_YZ_s, data_Z_ns, phi_character = "runif(1)")
#'
#' # Use summary statistics for both selected and non-selected sample ------------------------
#' # generate the list of summary statistics for selected sample (data_YZ_sumry)
#' data_YZ_sumry <- list(
#'   mean_YZ = c(mean(Y[index_s]), mean(Z[index_s])),
#'   var_YZ = var(cbind(Y[index_s], Z[index_s])),
#'   n_YZ = length(Y[index_s])
#' )
#' # generate the list of summary statistics for non-selected sample (data_Z_sumry)
#' data_Z_sumry <- list(
#'   mean_Z = mean(Z[index_ns]),
#'   var_Z = as.matrix(var(Z[index_ns])),
#'   n_Z = length(Z[index_ns])
#' )
#' colnames(data_YZ_sumry$var_YZ) <- c("Y", "Z")
#' colnames(data_Z_sumry$var_Z) <- c("Z")
#' result_sumry <- means_bayes(data_YZ_sumry, data_Z_sumry, phi_character = "runif(1)")
means_bayes <- function(data_YZ_s,
                        data_Z_ns,
                        phi = 0,
                        phi_character = NULL,
                        ndraws = 1500) {
  sumry_list <- ErrorCheck(data_YZ_s, data_Z_ns, prop_check = FALSE)
  data_YZ_s <- sumry_list$data_YZ_s
  sumry_YZ_s <- sumry_list$sumry_YZ_s
  sumry_Z_ns <- sumry_list$sumry_Z_ns

  # Selected sample size
  n_s <- sumry_YZ_s$n_YZ
  # Non-Selected sample size
  n_ns <- sumry_Z_ns$n_Z

  var_YZ_s <- sumry_YZ_s$var_YZ
  mean_YZ_s <- sumry_YZ_s$mean_YZ
  var_Z_s <- var_YZ_s[-1, -1]
  mean_Z_s <- mean_YZ_s[-1]
  zparams <- length(mean_YZ_s) - 1

  mean_Z_ns <- sumry_Z_ns$mean_Z
  var_Z_ns <- sumry_Z_ns$var_Z

  D.s <- solve(var_Z_s * (n_s - 1))
  C.s <- drop(crossprod(mean_Z_s, D.s))
  F.s <- drop(C.s %*% mean_Z_s)

  # Calculate mult.mat.s as (Z'Z)^{-1}
  mult.mat.s <- matrix(nrow = 1 + zparams, ncol = 1 + zparams)
  mult.mat.s[1, 1] <- F.s + 1 / n_s
  mult.mat.s[1, -1] <- mult.mat.s[-1, 1] <- -C.s
  mult.mat.s[-1, -1] <- D.s

  # MLEs of regression coefficients for Y|Z,S=1
  beta_YZ.Z_s <- drop(var_YZ_s[1, -1] %*% solve(var_YZ_s[-1, -1])) # slopes
  beta_YZ.Z_s <- c(mean_YZ_s[1] - beta_YZ.Z_s %*% mean_YZ_s[-1], beta_YZ.Z_s) # intercept
  # Residual variance
  var_Y.Z_s <- ((n_s - 1) / (n_s - (zparams + 1))) * drop(var_YZ_s[1, 1] - var_YZ_s[1, -1] %*% solve(var_YZ_s[-1, -1]) %*% var_YZ_s[-1, 1])

  draws <- matrix(nrow = ndraws, ncol = 10)

  # posterior draws for sigma2
  DRAWS_var_Y.Z_s <- (n_s - (zparams + 1)) * var_Y.Z_s / rchisq(ndraws, n_s - (zparams + 1))
  # posterior draws for beta | sigma2
  DRAWS_beta_YZ.Z_s <- t(apply(matrix(DRAWS_var_Y.Z_s), 1, function(s) rmnorm(mean = beta_YZ.Z_s, varcov = s * mult.mat.s)))


  # Compute predicted values of X for all cases in subpopulation
  DRAWS_var_X_s <- rowSums(tcrossprod(DRAWS_beta_YZ.Z_s[, -1], var_YZ_s[-1, -1]) * DRAWS_beta_YZ.Z_s[, -1])
  DRAWS_var_X_ns <- double(ndraws)
  DRAWS_var_XY_s <- matrix(0, nrow = ndraws, ncol = 3)
  DRAWS_S_X_ns <- double(ndraws)
  DRAWS_mean_X_s <- double(ndraws)
  DRAWS_mean_XY_s <- matrix(0, nrow = ndraws, ncol = 2)
  DRAWS_mean_X_ns <- double(ndraws)

  DRAWS_SMUB0_muY_ns <- double(ndraws)
  DRAWS_SMUB0_sigmaYY_ns <- double(ndraws)
  DRAWS_SMUB0_cov_XY_ns <- double(ndraws)
  DRAWS_SMUB0_muY <- double(ndraws)

  for (d in 1:ndraws)
  {
    pd <- 0
    while (pd == 0) {
      # Calculate the observed sample variance-covariance matrix for X and Y (selected sample)
      var_XY_s_d <- matrix(0, 2, 2)
      var_XY_s_d[1, 1] <- DRAWS_var_X_s[d]
      var_XY_s_d[1, 2] <- var_XY_s_d[2, 1] <- DRAWS_beta_YZ.Z_s[d, -1] %*% var_YZ_s[1, -1]
      var_XY_s_d[2, 2] <- var_YZ_s[1, 1]

      # Calculate S_xxns^(0)(d)
      S_X_ns <- sum(tcrossprod(DRAWS_beta_YZ.Z_s[d, -1], var_Z_ns) * DRAWS_beta_YZ.Z_s[d, -1])
      # Draw cov(x,y) matrix for selected from inverse wishart given the observed sample variance-covariance matrix
      DRAW_var_XY_s <- riwish(n_s - 1, var_XY_s_d) * (n_s - 1)
      # Draw var_X_ns for non-selected from inverse chi-square given the observed sample variance of proxy X
      DRAW_var_X_ns <- (n_ns - 1) * S_X_ns / rchisq(1, n_ns - 1)
      # Calculate mean_X_s = beta_y0.z^(d)+beta_yz.z^(d)*mean_Z_s
      mean_X_s <- DRAWS_beta_YZ.Z_s[d, ] %*% c(1, mean_YZ_s[-1])
      # Draw mean_XY_s for selected sample
      DRAW_mean_XY_s <- rmnorm(mean = c(mean_X_s, mean_YZ_s[1]), varcov = DRAW_var_XY_s / n_s)
      # Draw mean_X_ns for non-selected sample
      DRAW_mean_X_ns <- rnorm(1, DRAWS_beta_YZ.Z_s[d, ] %*% c(1, mean_Z_ns), sqrt(DRAW_var_X_ns / n_ns))

      # Draw phi from prior specified by phi_character
      if (!is.null(phi_character)) {
        # Draw phi using provided string
        phi <- eval(parse(text = phi_character))
        phi <- pmax(.Machine$double.eps, phi)
      }

      # Calculate rho_XY_s
      DRAW_rho_XY_s <- DRAW_var_XY_s[1, 2] / sqrt(DRAW_var_XY_s[1, 1] * DRAW_var_XY_s[2, 2])
      # Calculate g(phi) multiplier
      g_d <- ((phi + (1 - phi) * DRAW_rho_XY_s) / ((1 - phi) + phi * DRAW_rho_XY_s))
      # Calculate ratio of variance components for Y|Z, X|Z from sigmaYY_s/sigmaXX_s
      vratio_d <- DRAW_var_XY_s[2, 2] / DRAW_var_XY_s[1, 1]

      # Calculate Regression coefficients for non-selected cases, Y|Z,S=0
      DRAW_mean_Y_ns <- DRAW_mean_XY_s[2] + g_d * sqrt(vratio_d) * (DRAW_mean_X_ns - DRAW_mean_XY_s[1])
      DRAW_var_Y_ns <- DRAW_var_XY_s[2, 2] + g_d^2 * vratio_d * (DRAW_var_X_ns - DRAW_var_XY_s[1, 1])
      DRAW_cov_XY_ns <- DRAW_var_XY_s[1, 2] + g_d * sqrt(vratio_d) * (DRAW_var_X_ns - DRAW_var_XY_s[1, 1])

      if (DRAW_var_X_ns * DRAW_var_Y_ns - DRAW_cov_XY_ns^2 > 0) { # check pd for non-select cov XY
        pd <- 1 # pd check pass! GOOD!!
      } else {
        DRAWS_beta_YZ.Z_s[d, ] <- rmnorm(mean = beta_YZ.Z_s, varcov = DRAWS_var_Y.Z_s[d] * mult.mat.s)
        DRAWS_var_X_s[d] <- sum(tcrossprod(DRAWS_beta_YZ.Z_s[d, -1], var_YZ_s[-1, -1]) * DRAWS_beta_YZ.Z_s[d, -1])
      }
    }


    # save into the DRAWS matrix
    DRAWS_S_X_ns[d] <- S_X_ns
    DRAWS_var_XY_s[d, ] <- DRAW_var_XY_s[c(1, 2, 4)]
    DRAWS_var_X_ns[d] <- DRAW_var_X_ns
    DRAWS_mean_X_s[d] <- mean_X_s
    DRAWS_mean_XY_s[d, ] <- DRAW_mean_XY_s
    DRAWS_mean_X_ns[d] <- DRAW_mean_X_ns

    DRAW_Y_all <- (n_s * DRAW_mean_XY_s[2] + n_ns * DRAW_mean_Y_ns) / (n_ns + n_s)
    DRAW_SMUB <- (mean_YZ_s[1] - DRAW_Y_all) / sqrt(DRAW_var_XY_s[2, 2])

    # Calculation of SMUB(0)
    DRAW_SMUB0_scale <- DRAW_rho_XY_s * sqrt(vratio_d)
    DRAW_SMUB0_mean_Y_ns <- DRAW_mean_XY_s[2] + DRAW_SMUB0_scale * (DRAW_mean_X_ns - DRAW_mean_XY_s[1])
    DRAW_SMUB0_var_Y_ns <- DRAW_var_XY_s[2, 2] + DRAW_SMUB0_scale^2 * (DRAW_var_X_ns - DRAW_var_XY_s[1, 1])
    DRAW_SMUB0_cov_XY_ns <- DRAW_var_XY_s[1, 2] + DRAW_SMUB0_scale * (DRAW_var_X_ns - DRAW_var_XY_s[1, 1])
    DRAW_SMUB0_Y_all <- (n_s * DRAW_mean_XY_s[2] + n_ns * DRAW_SMUB0_mean_Y_ns) / (n_ns + n_s)
    DRAW_SMUB0 <- (mean_YZ_s[1] - DRAW_SMUB0_Y_all) / sqrt(DRAW_var_XY_s[2, 2])

    DRAWS_SMUB0_muY_ns[d] <- DRAW_SMUB0_mean_Y_ns
    DRAWS_SMUB0_sigmaYY_ns[d] <- DRAW_SMUB0_var_Y_ns
    DRAWS_SMUB0_cov_XY_ns[d] <- DRAW_SMUB0_cov_XY_ns
    DRAWS_SMUB0_muY[d] <- DRAW_SMUB0_Y_all

    drawsPPM <- list(
      phi = phi,
      muY_s = DRAW_mean_XY_s[2],
      muY_ns = DRAW_mean_Y_ns,
      sigmaYY_s = DRAW_var_XY_s[2, 2],
      sigmaYY_ns = DRAW_var_Y_ns,
      rho_XY_s = DRAW_rho_XY_s,
      muY = DRAW_Y_all,
      smub = DRAW_SMUB,
      smub0 = DRAW_SMUB0,
      smab = DRAW_SMUB - DRAW_SMUB0
    )
    draws[d, ] <- unlist(drawsPPM)
  }
  draws <- as.data.frame(draws)
  names(draws) <- names(drawsPPM)

  return(draws)
}


#' Calculate maximum likelihood estimates for continuous outcomes
#'
#' @param data_YZ_s list for selected sample (Microdata or Summary statistics). See details below for more details.
#' @param data_Z_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param phi scalar or vector of phi values at which the (S)MUB is calculated, should be in `[0,1]` but values from `[-Inf,1]` are allowed (defaults to `c(0, 0.5, 1)`)
#'
#' @returns A `length(phi)*10` matrix of Bayes estimates with each row containing:
#' \itemize{
#' \item `phi` : corresponding phi
#' \item `muY_s` : ML estimate of the mean of outcome variable for selected sample
#' \item `muY_ns` : ML estimate of the mean of outcome variable for non-selected sample
#' \item `sigmaYY_s` : ML estimate of the variance of outcome variable for selected sample
#' \item `sigmaYY_ns` : ML estimate of the variance of outcome variable for non-selected sample
#' \item `rho_XY_s` : ML estimate of the correlation between the proxy and the outcome variable for selected sample
#' \item `muY` : ML estimate of the population mean of outcome variable
#' \item `smub` : ML estimate of SMUB with phi = phi
#' \item `smub0` : ML estimate of SMUB with phi = 0
#' \item `smab` :  ML estimate of SMAB
#' }
#' @import stats
#' @import MCMCpack
#' @import mnormt
#' @import MASS
#' @import mvtnorm
#' @importFrom stringr str_sub
#' @export
#' @details
#' \describe{
#'  \item{`data_YZ_s=list(Y,Z)`(Microdata)}{
#'    \itemize{
#'      \item `Y` : `n*1` matrix of outcome variable `Y`
#'      \item `Z` : `n*nZparams` matrix of auxiliary variables `Z`
#'    }
#'  }
#'  \item{`data_YZ_s=list(mean_YZ,var_YZ,n_YZ)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_YZ` : `(1+nZparams)*1` vector of mean (in order) of outcome variable `Y` and auxiliary variables `Z`
#'      \item `var_YZ` : `(1+nZparams)*(1+nZparams)` matrix of variance (in order) of outcome variable `Y` and auxiliary variables `Z`
#'      \item `n_YZ` : number of selected sample
#'    }
#'  }
#'  \item{`data_Z_ns=list(Z)`(Microdata)}{
#'    \itemize{
#'      \item `Z` : `(N-n)*nZparams` matrix of auxiliary variables `Z`
#'    }
#'  }
#'  \item{`data_Z_ns=list(mean_Z, var_Z, n_Z)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_Z` : `nZparams*1` vector of mean of auxiliary variables `Z`
#'      \item `var_Z` : `nZparams*nZparams` matrix of variance of auxiliary variables `Z`
#'      \item `n_Z` : number of non-selected sample
#'    }
#'  }
#' }
#'
#' @seealso [means_bayes()]
#' @seealso [prop_bayes()]
#' @seealso [prop_mle()]
#'
#' @examples
#' set.seed(123)
#' # num of total sample
#' N <- 10000
#' # generate auxiliary variables Z
#' Z <- rnorm(N, 0, 1)
#' # define the biserial correlation rho between outcome variable Y and the proxy X
#' rho <- 0.5
#' # define the mean of the outcome variable Y
#' mu_y <- 10
#' a1 <- rho / sqrt(1 - rho^2)
#' # generate the observed variable Y
#' Y <- rnorm(Z, 10 + a1 * Z, 1)
#' # generate the sample selection indicator S with 0.05 of the population selected
#' prob <- plogis(-4 + 0.3 * Z + 0.1 * Y)
#' S <- rbinom(prob, 1, prob)
#' index_s <- S == 1
#' index_ns <- S == 0
#' # generate the list of microdata for selected (data_YZ_s) and non-selected sample (data_Z_ns)
#' data_YZ_s <- list(Y = Y[index_s], Z = Z[index_s])
#' data_Z_ns <- list(Z = Z[index_ns])
#'
#' # calculate result from ppmm model with maximum likelihood estimation
#' result <- means_mle(data_YZ_s, data_Z_ns)
#'
#' # Use summary statistics for both selected and non-selected sample ------------------------
#' # generate the list of summary statistics for selected sample (data_YZ_sumry)
#' data_YZ_sumry <- list(
#'   mean_YZ = c(mean(Y[index_s]), mean(Z[index_s])),
#'   var_YZ = var(cbind(Y[index_s], Z[index_s])),
#'   n_YZ = length(Y[index_s])
#' )
#' # generate the list of summary statistics for non-selected sample (data_Z_sumry)
#' data_Z_sumry <- list(
#'   mean_Z = mean(Z[index_ns]),
#'   var_Z = as.matrix(var(Z[index_ns])),
#'   n_Z = length(Z[index_ns])
#' )
#' colnames(data_YZ_sumry$var_YZ) <- c("Y", "Z")
#' colnames(data_Z_sumry$var_Z) <- c("Z")
#' result_sumry <- means_mle(data_YZ_sumry, data_Z_sumry)
#'
means_mle <- function(data_YZ_s,
                      data_Z_ns,
                      phi = c(0, 0.5, 1)) {
  sumry_list <- ErrorCheck(data_YZ_s, data_Z_ns)
  sumry_YZ_s <- sumry_list$sumry_YZ_s
  sumry_Z_ns <- sumry_list$sumry_Z_ns

  # Sufficient statistics: non-selected/missing population
  mean_Z_ns <- sumry_Z_ns$mean_Z
  var_Z_ns <- sumry_Z_ns$var_Z
  n0 <- sumry_Z_ns$n_Z

  # Sufficient statistics: selected/non-missing population
  mean_YZ_s <- sumry_YZ_s$mean_YZ
  var_YZ_s <- sumry_YZ_s$var_YZ
  n1 <- sumry_YZ_s$n_YZ
  # First step calculates the slopes
  beta_YZ.Z_s <- drop(var_YZ_s[1, -1] %*% solve(var_YZ_s[-1, -1]))
  # Second step calculates the intercept
  beta_YZ.Z_s <- c(mean_YZ_s[1] - beta_YZ.Z_s %*% mean_YZ_s[-1], beta_YZ.Z_s)
  npreds_as_dummy <- length(mean_YZ_s)
  var_Y.Z_s <- (n1 / (n1 - npreds_as_dummy + 1)) * drop(var_YZ_s[1, 1] - var_YZ_s[1, -1] %*% solve(var_YZ_s[-1, -1]) %*% var_YZ_s[-1, 1])
  ZtZinv_s <- solve(rbind(c(n1, n1 * mean_YZ_s[-1]), cbind(n1 * mean_YZ_s[-1], (var_YZ_s[-1, -1] * (n1 - 1) + n1 * tcrossprod(mean_YZ_s[-1])))))

  ## NonBayes-calculation: point estimates of SMUB
  mean_X_s <- drop(beta_YZ.Z_s %*% c(1, mean_YZ_s[-1]))
  mean_X_ns <- drop(beta_YZ.Z_s %*% c(1, mean_Z_ns))
  var_X_s <- drop(tcrossprod(beta_YZ.Z_s[-1], var_YZ_s[-1, -1]) %*% beta_YZ.Z_s[-1])
  var_X_ns <- drop(tcrossprod(beta_YZ.Z_s[-1], var_Z_ns) %*% beta_YZ.Z_s[-1])
  mean_X_pop <- (mean_X_s * n1 + mean_X_ns * n0) / (n0 + n1)
  cor_XY_s <- drop(beta_YZ.Z_s[-1] %*% var_YZ_s[1, -1]) / sqrt(var_X_s * var_YZ_s[1, 1])

  g.phi <- (phi + (1 - phi) * cor_XY_s) / (phi * cor_XY_s + (1 - phi))
  vratio <- var_YZ_s[1, 1] / var_X_s

  mean_Y_s <- mean_YZ_s[1]

  ## NonBayes-calculation: point estimates of Y|Z,S=0 (mean and variance)
  mean_Y_ns <- mean_Y_s + g.phi * sqrt(vratio) * (mean_X_ns - mean_X_s)
  sigmaYY_ns <- var_YZ_s[1, 1] + g.phi^2 * vratio * (var_X_ns - var_X_s)

  ## NonBayes-calculation: point estimates of SMUB(0)
  smub0_point_est <- cor_XY_s * (mean_X_s - mean_X_pop) / sqrt(var_X_s)

  ## NonBayes-calculation: point estimates of SMAB
  smub_point_est <- g.phi * (mean_X_s - mean_X_pop) / sqrt(var_X_s)

  ## NonBayes-calculation: point estimates of MUB
  mub_point_est <- g.phi * sqrt(vratio) * (mean_X_s - mean_X_pop)
  smab_point_est <- smub_point_est - smub0_point_est

  ## population Y
  mean_Y_pop <- mean_Y_s - mub_point_est

  ## final result
  sumry_table <- data.frame(
    phi = phi,
    muY_s = mean_Y_s,
    muY_ns = mean_Y_ns,
    sigmaYY_s = var_YZ_s[1, 1],
    sigmaYY_ns = sigmaYY_ns,
    rho_XY_s = cor_XY_s,
    muY = mean_Y_pop,
    smub = smub_point_est,
    smub0 = smub0_point_est,
    smab = smab_point_est,
    row.names = NULL
  )

  ## return final result
  return(sumry_table)
}
