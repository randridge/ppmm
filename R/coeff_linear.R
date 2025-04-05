#' Function to calculate MUB for non-binary outcomes
#'
#' `mub_coeff_linear_mle` calculate maximum likelihood estimates of MUB for non-binary outcomes
#'
#' @param data_YZA_s list for selected sample (Microdata or Summary statistics). See details below for more details.
#' @param data_ZA_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param nZvars dimension of predictor variables `Z`
#' @param phi scalar or vector of phi values at which the (S)MUB is calculated, should be in `[0,1]` but values from `[-Inf,1]` are allowed
#'
#' @returns A dataframe with dimension of `ndraws*nZvars` with each row the values of MUB for each Z variable with the given phi (Column 1).
#' @import stats
#' @import MCMCpack
#' @import mnormt
#' @import MASS
#' @import mvtnorm
#' @importFrom stringr str_sub
#' @export
#' @details
#' \describe{
#'  \item{`data_YZA_s=list(Y,Z,A)`(Microdata)}{
#'    \itemize{
#'      \item `Y` : `n*1` matrix of outcome variable `Y`
#'      \item `Z` : `n*nZvars` matrix of predictor variables `Z`
#'      \item `A` : `n*nAvars`matrix of auxiliary variables `A`
#'    }
#'  }
#'  \item{`data_YZA_s=list(mean_YZA, var_YZA, n_YZA)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_YZA` : `(1+nZvars+nAvars)*1` vector of mean (in order) of outcome variable `Y`, predictor variables `Z` and auxiliary variables `A`
#'      \item `var_YZA` : `(1+nZvars+nAvars)*(1+nZvars+nAvars)` matrix of variance (in order) of outcome variable `Y`, predictor variables `Z` and auxiliary variables `A`
#'      \item `n_YZA` : number of selected sample
#'    }
#'  }
#'
#'  \item{`data_ZA_ns=list(Z,A)`(Microdata)}{
#'    \itemize{
#'      \item `Z` : `(N-n)*nZvars` matrix of predictor variables `Z`
#'      \item `A` : `(N-n)*nAvars`matrix of auxiliary variables `A`
#'    }
#'  }
#'  \item{`data_ZA_ns=list(mean_ZA, var_ZA, n_ZA)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_ZA` : `(nZvars+nAvars)*1` vector of mean (in order) of predictor variables `Z` and auxiliary variables `A`
#'      \item `var_ZA` : `(nZvars+nAvars)*(nZvars+nAvars)` matrix of variance (in order) of predictor variables `Z` and auxiliary variables `A`
#'      \item `n_ZA` : number of non-selected sample
#'    }
#'  }
#'
#' }
#'
#' @references
#' West, B. T., Little, R. J., Andridge, R. R., Boonstra, P. S., Ware, E. B., Pandit, A., & Alvarado-Leiton, F. (2021). Assessing selection bias in regression coefficients estimated from nonprobability samples with applications to genetics and demographic surveys. \emph{The annals of applied statistics}, 15(3), 1556–1581.
#'
#' @examples
#' require(mvtnorm)
#' set.seed(123)
#' # num of total sample
#' N <- 10000
#' nZvars <- 2
#' mu <- c(10, 0, 0, 0)
#' rhoy1 <- 0.4
#' rhoy2 <- 0.6
#' rho1a <- 0.2
#' sigmaya <- 0.3
#' Sigma <- matrix(c(
#'   4, 2 * rhoy1, 2 * rhoy2, sigmaya,
#'   2 * rhoy1, 1, 0, rho1a,
#'   2 * rhoy2, 0, 1, 0,
#'   sigmaya, rho1a, 0, 1
#' ), 4, 4)
#' raw_data <- rmvnorm(N, mu, Sigma)
#' colnames(raw_data) <- c("Y", "Z1", "Z2", "A")
#' Y <- raw_data[, 1]
#' Z <- raw_data[, 2:3]
#' A <- raw_data[, 4]
#' prob <- plogis(-4 + Z %*% c(log(1.1), log(1.1)) + log(1.1) * Y + log(2) * A)
#' S <- rbinom(prob, 1, prob)
#' index_s <- S == 1
#' index_ns <- S == 0
#' # generate the list of microdata for selected (data_YZA_s) and non-selected sample (data_ZA_ns)
#' data_YZA_s <- list(Y = Y[index_s], Z = Z[index_s, ], A = A[index_s])
#' data_ZA_ns <- list(Z = Z[index_ns, ], A = A[index_ns])
#'
#' # draw mle result from ppmm model
#' result <- mub_coeff_linear_mle(data_YZA_s, data_ZA_ns, nZvars)
#'
#' # Use summary statistics for both selected and non-selected sample ------------------------
#' # generate the list of summary statistics for selected sample (data_YZA_sumry)
#' data_YZA_sumry <- list(
#'   mean_YZA = colMeans(raw_data[index_s, ]),
#'   var_YZA = var(raw_data[index_s, ]),
#'   n_YZA = sum(index_s)
#' )
#' # generate the list of summary statistics for non-selected sample (data_ZA_sumry)
#' data_ZA_sumry <- list(
#'   mean_ZA = colMeans(raw_data[index_ns, -1]),
#'   var_ZA = var(raw_data[index_ns, -1]),
#'   n_ZA = sum(index_ns)
#' )
#' result_sumry <- mub_coeff_linear_mle(data_YZA_sumry, data_ZA_sumry, nZvars)
#'
mub_coeff_linear_mle <- function(data_YZA_s,
                                 data_ZA_ns,
                                 nZvars,
                                 phi = c(0, 0.5, 1)) {
  sumry_list <- Coeff_ErrorCheck(data_YZA_s, data_ZA_ns, prop_check = FALSE)
  sumry_YZA_s <- sumry_list$sumry_YZA_s
  sumry_ZA_ns <- sumry_list$sumry_ZA_ns

  # Pieces needed for calculations below
  n_s <- sumry_YZA_s$n_YZA
  n_ns <- sumry_ZA_ns$n_ZA
  # (Y,Z,A|S=1)
  var_YZA_s <- sumry_YZA_s$var_YZA
  # (Y,Z|S=1)
  mean_YZ_s <- sumry_YZA_s$mean_YZA[1:(1 + nZvars)]
  var_YZ_s <- sumry_YZA_s$var_YZA[1:(1 + nZvars), 1:(1 + nZvars)]
  # (Y,A|S=1)
  mean_YA_s <- sumry_YZA_s$mean_YZA[-c(2:(nZvars + 1))]
  var_YA_s <- sumry_YZA_s$var_YZA[-c(2:(nZvars + 1)), -c(2:(nZvars + 1))]
  # (Z,A|S=1)
  var_ZA_s <- sumry_YZA_s$var_YZA[-1, -1]
  # (A|S=0)
  mean_A_ns <- sumry_ZA_ns$mean_ZA[-c(1:nZvars)]
  var_A_ns <- sumry_ZA_ns$var_ZA[-c(1:nZvars), -c(1:nZvars)]
  # (Z|S=0)
  mean_Z_ns <- sumry_ZA_ns$mean_ZA[1:nZvars]
  # (Z,A|S=0)
  var_ZA_ns <- sumry_ZA_ns$var_ZA

  #### Regression of Y|Z,S=1 --> regression of interest
  # First step calculates the slopes
  beta_YZ.Z_s <- drop(var_YZ_s[1, -1] %*% solve(var_YZ_s[-1, -1]))
  # Second step calculates the intercept
  beta_YZ.Z_s <- c(mean_YZ_s[1] - beta_YZ.Z_s %*% mean_YZ_s[-1], beta_YZ.Z_s)
  # Residual variance
  npreds_as_dummy <- nZvars + 1
  var_Y.Z_s <- (n_s / (n_s - npreds_as_dummy)) * drop(var_YZ_s[1, 1] - var_YZ_s[1, -1] %*% solve(var_YZ_s[-1, -1]) %*% var_YZ_s[-1, 1])

  #### Regression of Y|A,Z,S=1 --> to calculate proxy X
  # First step calculates the slopes (Don't need the intercept)
  beta_YZA.ZA_s <- drop(var_YZA_s[1, -1] %*% solve(var_YZA_s[-1, -1]))
  # Subset to only the A terms
  beta_YA.ZA_s <- beta_YZA.ZA_s[-c(1:nZvars)]

  #### Means and variances for (X,Z|S)
  # Selected
  mean_XZ_s <- c(mean_YA_s[-1] %*% beta_YA.ZA_s, mean_YZ_s[-1])
  var_XZ_s <- var_YZ_s # start with this and replace the Y with X
  var_XZ_s[1, 1] <- drop(crossprod(beta_YA.ZA_s, var_YA_s[-1, -1]) %*% beta_YA.ZA_s)
  var_XZ_s[1, -1] <- var_XZ_s[-1, 1] <- drop(crossprod(beta_YA.ZA_s, var_ZA_s[-c(1:nZvars), 1:nZvars]))

  # Non-selected
  mean_XZ_ns <- c(mean_A_ns %*% beta_YA.ZA_s, mean_Z_ns)
  var_XZ_ns <- matrix(nrow = 1 + nZvars, ncol = 1 + nZvars)
  var_XZ_ns[1, 1] <- drop(crossprod(beta_YA.ZA_s, var_A_ns) %*% beta_YA.ZA_s)
  var_XZ_ns[1, -1] <- var_XZ_ns[-1, 1] <- drop(crossprod(beta_YA.ZA_s, var_ZA_ns[-c(1:nZvars), 1:nZvars]))
  var_XZ_ns[-1, -1] <- var_ZA_ns[1:nZvars, 1:nZvars]

  #### Regression of X|Z,S=1
  # First step calculates the slopes
  beta_XZ.Z_s <- drop(var_XZ_s[1, -1] %*% solve(var_XZ_s[-1, -1]))
  # Second step calculates the intercept
  beta_XZ.Z_s <- c(mean_XZ_s[1] - beta_XZ.Z_s %*% mean_XZ_s[-1], beta_XZ.Z_s)
  # Residual variance
  npreds_as_dummy <- nZvars + 1
  var_X.Z_s <- (n_s / (n_s - npreds_as_dummy)) * drop(var_XZ_s[1, 1] - var_XZ_s[1, -1] %*% solve(var_XZ_s[-1, -1]) %*% var_XZ_s[-1, 1])

  #### Regression of X|Z,S=0
  # First step calculates the slopes
  beta_XZ.Z_ns <- drop(var_XZ_ns[1, -1] %*% solve(var_XZ_ns[-1, -1]))
  # Second step calculates the intercept
  beta_XZ.Z_ns <- c(mean_XZ_ns[1] - beta_XZ.Z_ns %*% mean_XZ_ns[-1], beta_XZ.Z_ns)

  #### Cor(X,Y|S==1) *not conditional on Z
  cor_XY_s <- drop(beta_YA.ZA_s %*% var_YA_s[1, -1]) / sqrt(var_XZ_s[1, 1] * var_YZ_s[1, 1])

  #### Cor(X,Y|Z,S==1) *conditional on Z
  cor_XY.Z_s <- drop(beta_YA.ZA_s %*% var_YA_s[1, -1] - var_XZ_s[1, -1] %*% solve(var_XZ_s[-1, -1]) %*% var_YZ_s[1, -1]) / sqrt(var_X.Z_s * var_Y.Z_s)

  #### MUB
  g <- (phi + (1 - phi) * cor_XY.Z_s) / (phi * cor_XY.Z_s + (1 - phi))
  g[phi == -Inf] <- -1 # replaces the NaN that results from the calculation involving -Inf
  mub_point_est <- as.matrix(cbind(phi = phi, outer(g, (beta_XZ.Z_s - beta_XZ.Z_ns) * sqrt(var_Y.Z_s / var_X.Z_s))))

  # add column names
  for (j in 0:nZvars)
  {
    colnames(mub_point_est)[2 + j] <- paste("index.B", j, sep = "")
  }

  return(mub_point_est)
}

#' Function to calculate Bayes estimates for non-binary outcomes
#'
#' `mub_coeff_linear_bayes` calculate Bayes estimates of MUB for non-binary outcomes
#'
#' @param data_YZA_s list for selected sample (Microdata or Summary statistics). See details below for more details.
#' @param data_ZA_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param nZvars dimension of predictor variables `Z`
#' @param phi value of phi, fixed to a value between 0 and 1 (defaults to 0)
#' @param phi_character draw phi from specific distribution (i.e."`runif(1)`", defaults to `NULL`)
#' * If `phi_character=NULL`, fix at value set to `phi`;
#' * If `phi_character=` specific distribution (i.e. `runif(1)` or `rbeta(1,1,1)`), then the value of `phi` input to the function is ignored
#' @param ndraws number of draws
#'
#' @returns A dataframe with dimension of `ndraws*nZvars` with each row a draw of MUB for each Z variable.
#' @import stats
#' @import MCMCpack
#' @import mnormt
#' @import MASS
#' @import mvtnorm
#' @importFrom stringr str_sub
#' @export
#' @details
#' \describe{
#'  \item{`data_YZA_s=list(Y,Z,A)`(Microdata)}{
#'    \itemize{
#'      \item `Y` : `n*1` matrix of outcome variable `Y`
#'      \item `Z` : `n*nZvars` matrix of predictor variables `Z`
#'      \item `A` : `n*nAvars`matrix of auxiliary variables `A`
#'    }
#'  }
#'  \item{`data_YZA_s=list(mean_YZA, var_YZA, n_YZA)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_YZA` : `(1+nZvars+nAvars)*1` vector of mean (in order) of outcome variable `Y`, predictor variables `Z` and auxiliary variables `A`
#'      \item `var_YZA` : `(1+nZvars+nAvars)*(1+nZvars+nAvars)` matrix of variance (in order) of outcome variable `Y`, predictor variables `Z` and auxiliary variables `A`
#'      \item `n_YZA` : number of selected sample
#'    }
#'  }
#'
#'  \item{`data_ZA_ns=list(Z,A)`(Microdata)}{
#'    \itemize{
#'      \item `Z` : `(N-n)*nZvars` matrix of predictor variables `Z`
#'      \item `A` : `(N-n)*nAvars`matrix of auxiliary variables `A`
#'    }
#'  }
#'  \item{`data_ZA_ns=list(mean_ZA, var_ZA, n_ZA)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_ZA` : `(nZvars+nAvars)*1` vector of mean (in order) of predictor variables `Z` and auxiliary variables `A`
#'      \item `var_ZA` : `(nZvars+nAvars)*(nZvars+nAvars)` matrix of variance (in order) of predictor variables `Z` and auxiliary variables `A`
#'      \item `n_ZA` : number of non-selected sample
#'    }
#'  }
#'
#' }
#' @seealso [mub_coeff_binary_bayes()]
#'
#' @references
#' West, B. T., Little, R. J., Andridge, R. R., Boonstra, P. S., Ware, E. B., Pandit, A., & Alvarado-Leiton, F. (2021). Assessing selection bias in regression coefficients estimated from nonprobability samples with applications to genetics and demographic surveys. \emph{The annals of applied statistics}, 15(3), 1556–1581.
#'
#' @examples
#' require(mvtnorm)
#' set.seed(123)
#' # num of total sample
#' N <- 10000
#' nZvars <- 2
#' mu <- c(10, 0, 0, 0)
#' rhoy1 <- 0.4
#' rhoy2 <- 0.6
#' rho1a <- 0.2
#' sigmaya <- 0.3
#' Sigma <- matrix(c(
#'   4, 2 * rhoy1, 2 * rhoy2, sigmaya,
#'   2 * rhoy1, 1, 0, rho1a,
#'   2 * rhoy2, 0, 1, 0,
#'   sigmaya, rho1a, 0, 1
#' ), 4, 4)
#' raw_data <- rmvnorm(N, mu, Sigma)
#' colnames(raw_data) <- c("Y", "Z1", "Z2", "A")
#' Y <- raw_data[, 1]
#' Z <- raw_data[, 2:3]
#' A <- raw_data[, 4]
#' prob <- plogis(-4 + Z %*% c(log(1.1), log(1.1)) + log(1.1) * Y + log(2) * A)
#' S <- rbinom(prob, 1, prob)
#' index_s <- S == 1
#' index_ns <- S == 0
#' # generate the list of microdata for selected (data_YZA_s) and non-selected sample (data_ZA_ns)
#' data_YZA_s <- list(Y = Y[index_s], Z = Z[index_s, ], A = A[index_s])
#' data_ZA_ns <- list(Z = Z[index_ns, ], A = A[index_ns])
#'
#' # draw bayesian result from ppmm model
#' result <- mub_coeff_linear_bayes(data_YZA_s, data_ZA_ns, nZvars)
#'
#' # Use summary statistics for both selected and non-selected sample ------------------------
#' # generate the list of summary statistics for selected sample (data_YZA_sumry)
#' data_YZA_sumry <- list(
#'   mean_YZA = colMeans(raw_data[index_s, ]),
#'   var_YZA = var(raw_data[index_s, ]),
#'   n_YZA = sum(index_s)
#' )
#' # generate the list of summary statistics for non-selected sample (data_ZA_sumry)
#' data_ZA_sumry <- list(
#'   mean_ZA = colMeans(raw_data[index_ns, -1]),
#'   var_ZA = var(raw_data[index_ns, -1]),
#'   n_ZA = sum(index_ns)
#' )
#' result_sumry <- mub_coeff_linear_bayes(data_YZA_sumry, data_ZA_sumry, nZvars)
mub_coeff_linear_bayes <- function(data_YZA_s,
                                   data_ZA_ns,
                                   nZvars,
                                   phi = 0,
                                   phi_character = NULL,
                                   ndraws = 1500) {
  sumry_list <- Coeff_ErrorCheck(data_YZA_s, data_ZA_ns, prop_check = FALSE)
  sumry_YZA_s <- sumry_list$sumry_YZA_s
  sumry_ZA_ns <- sumry_list$sumry_ZA_ns

  # Pieces needed for calculations
  n_s <- sumry_YZA_s$n_YZA
  n_ns <- sumry_ZA_ns$n_ZA
  aparams <- length(sumry_YZA_s$mean_YZA) - 1 - nZvars
  # Y,Z,A|S=1
  mean_YZA_s <- sumry_YZA_s$mean_YZA
  var_YZA_s <- sumry_YZA_s$var_YZA
  # Z,A|S=0
  mean_ZA_ns <- sumry_ZA_ns$mean_ZA
  var_ZA_ns <- sumry_ZA_ns$var_ZA
  # Z,A|S=1
  mean_ZA_s <- mean_YZA_s[-1]
  var_ZA_s <- var_YZA_s[-1, -1]
  # Z|S=1
  mean_Z_s <- mean_YZA_s[2:(nZvars + 1)]
  var_Z_s <- var_YZA_s[2:(nZvars + 1), 2:(nZvars + 1)]
  # Z|S=0
  mean_Z_ns <- mean_ZA_ns[1:nZvars]
  var_Z_ns <- var_ZA_ns[1:nZvars, 1:nZvars]
  # Y,Z|S=1
  mean_YZ_s <- mean_YZA_s[1:(nZvars + 1)]
  var_YZ_s <- var_YZA_s[1:(nZvars + 1), 1:(nZvars + 1)]
  # Y,A|S=1
  mean_YA_s <- mean_YZA_s[-c(2:(nZvars + 1))]
  var_YA_s <- var_YZA_s[-c(2:(nZvars + 1)), -c(2:(nZvars + 1))]
  # A|S=0
  mean_A_ns <- mean_ZA_ns[-c(1:nZvars)]
  var_A_ns <- var_ZA_ns[-c(1:nZvars), -c(1:nZvars)]

  #### Regression of Y|Z,A,S=1 --> to calculate proxy X
  D.za.s <- solve(var_ZA_s * (n_s - 1)) # resulting matrix is (nZvars x nZvars)
  C.za.s <- drop(crossprod(mean_ZA_s, D.za.s)) # resulting matrix is (1 x nZvars)
  F.za.s <- drop(C.za.s %*% mean_ZA_s) # result is scalar
  mult.mat.za.s <- matrix(nrow = 1 + nZvars + aparams, ncol = 1 + nZvars + aparams)
  mult.mat.za.s[1, 1] <- F.za.s + 1 / n_s
  mult.mat.za.s[1, -1] <- mult.mat.za.s[-1, 1] <- -C.za.s
  mult.mat.za.s[-1, -1] <- D.za.s
  # MLE of slopes
  beta_YZA.ZA_s <- drop(var_YZA_s[1, -1] %*% solve(var_YZA_s[-1, -1])) # slopes
  beta_YZA.ZA_s <- c(mean_YZA_s[1] - beta_YZA.ZA_s %*% mean_ZA_s, beta_YZA.ZA_s) # intercept
  # MLE of Residual variance
  var_Y.ZA_s <- (n_s / (n_s - (nZvars + aparams))) * drop(var_YZA_s[1, 1] - var_YZA_s[1, -1] %*% solve(var_YZA_s[-1, -1]) %*% var_YZA_s[-1, 1])

  #### Regression of Y|Z,S=1 (target regresion of interest)
  # MLEs of regression coefficients for Y|Z,S=1
  beta_YZ.Z_s <- drop(var_YZ_s[1, -1] %*% solve(var_YZ_s[-1, -1])) # slopes
  beta_YZ.Z_s <- c(mean_YZ_s[1] - beta_YZ.Z_s %*% mean_YZ_s[-1], beta_YZ.Z_s) # intercept
  # Residual variance
  var_Y.Z_s <- (n_s / (n_s - nZvars)) * drop(var_YZ_s[1, 1] - var_YZ_s[1, -1] %*% solve(var_YZ_s[-1, -1]) %*% var_YZ_s[-1, 1])

  ##################
  # pieces involving Z needed for draws
  # selected cases
  D.s <- solve(var_Z_s * (n_s - 1)) # resulting matrix is (nZvars x nZvars)
  C.s <- drop(crossprod(mean_Z_s, D.s)) # resulting matrix is (1 x nZvars)
  F.s <- drop(C.s %*% mean_Z_s) # result is scalar
  mult.mat.s <- matrix(nrow = nZvars + 1, ncol = nZvars + 1)
  mult.mat.s[1, 1] <- F.s + 1 / n_s
  mult.mat.s[1, 2:(nZvars + 1)] <- mult.mat.s[2:(nZvars + 1), 1] <- -C.s
  mult.mat.s[-1, -1] <- D.s
  # non-selected cases
  D.ns <- solve(var_Z_ns * (n_ns - 1)) # resulting matrix is (nZvars x nZvars)
  C.ns <- drop(crossprod(mean_Z_ns, D.ns)) # resulting matrix is (1 x nZvars)
  F.ns <- drop(C.ns %*% mean_Z_ns) # result is scalar
  mult.mat.ns <- matrix(nrow = nZvars + 1, ncol = nZvars + 1)
  mult.mat.ns[1, 1] <- F.ns + 1 / n_ns
  mult.mat.ns[1, 2:(nZvars + 1)] <- mult.mat.ns[2:(nZvars + 1), 1] <- -C.ns
  mult.mat.ns[-1, -1] <- D.ns

  ###############
  # DRAWS STEP 1a: Draw residual Var(Y|Z,A,S=1) from inverse-ChiSq
  ###############
  DRAWS_var_Y.ZA_s <- (n_s - (nZvars + aparams + 1)) * var_Y.ZA_s / rchisq(ndraws, n_s - (nZvars + aparams + 1))

  ###############
  # DRAWS STEP 1b: Draw regression coefs from Y|Z,A,S=1 conditional on draws of resid var (columns are replicates)
  ###############
  DRAWS_beta_YZA.ZA_s <- apply(matrix(DRAWS_var_Y.ZA_s), 1, function(s) rmnorm(mean = beta_YZA.ZA_s, varcov = s * mult.mat.za.s))

  ###############
  # Loop over draws for remaining parameters
  ###############
  # Matrix to hold draws of MUB(phi)
  mubdraws <- matrix(0, nrow = ndraws, ncol = (nZvars + 1))
  for (d in 1:ndraws)
  {
    pd <- 0
    while (pd == 0) # Checking draws to ensure positive-definite covariance matrix
    {
      smat_invertible <- 0
      while (smat_invertible == 0) # Checking to ensure residual covariance matrix (X,Y|Z,S=1) is invertible
      {
        #### Means and variances for (X,Z|S)
        # Selected
        mean_XZ_s_d <- c(mean_YA_s[-1] %*% DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d], mean_YZ_s[-1])
        var_XZ_s_d <- var_YZ_s # start with this and replace the Y with X
        var_XZ_s_d[1, 1] <- drop(crossprod(DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d], var_YA_s[-1, -1]) %*% DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d])
        var_XZ_s_d[1, -1] <- var_XZ_s_d[-1, 1] <- drop(crossprod(DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d], var_ZA_s[-c(1:nZvars), 1:nZvars]))
        # Non-selected
        mean_XZ_ns_d <- c(mean_A_ns %*% DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d], mean_Z_ns)
        var_XZ_ns_d <- matrix(nrow = 1 + nZvars, ncol = 1 + nZvars)
        var_XZ_ns_d[1, 1] <- drop(crossprod(DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d], var_A_ns) %*% DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d])
        var_XZ_ns_d[1, -1] <- var_XZ_ns_d[-1, 1] <- drop(crossprod(DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d], var_ZA_ns[-c(1:nZvars), 1:nZvars]))
        var_XZ_ns_d[-1, -1] <- var_ZA_ns[1:nZvars, 1:nZvars]

        #### Regression of X|Z,S=1
        beta_XZ.Z_s_d <- drop(var_XZ_s_d[1, -1] %*% solve(var_XZ_s_d[-1, -1])) # slopes
        beta_XZ.Z_s_d <- c(mean_XZ_s_d[1] - beta_XZ.Z_s_d %*% mean_XZ_s_d[-1], beta_XZ.Z_s_d) # intercept
        # Residual variance
        var_X.Z_s_d <- (n_s / (n_s - nZvars)) * drop(var_XZ_s_d[1, 1] - var_XZ_s_d[1, -1] %*% solve(var_XZ_s_d[-1, -1]) %*% var_XZ_s_d[-1, 1])

        #### Residual covariance of X,Y|Z,S=1
        var_XY.Z_s_d <- drop(DRAWS_beta_YZA.ZA_s[-c(1:(nZvars + 1)), d] %*% var_YA_s[1, -1] - var_XZ_s_d[1, -1] %*% solve(var_XZ_s_d[-1, -1]) %*% var_YZ_s[1, -1])

        #### Regression of X|Z,S=0
        beta_XZ.Z_ns_d <- drop(var_XZ_ns_d[1, -1] %*% solve(var_XZ_ns_d[-1, -1])) # slopes
        beta_XZ.Z_ns_d <- c(mean_XZ_ns_d[1] - beta_XZ.Z_ns_d %*% mean_XZ_ns_d[-1], beta_XZ.Z_ns_d) # intercept
        # Residual variance
        var_X.Z_ns_d <- (n_ns / (n_ns - nZvars)) * drop(var_XZ_ns_d[1, 1] - var_XZ_ns_d[1, -1] %*% solve(var_XZ_ns_d[-1, -1]) %*% var_XZ_ns_d[-1, 1])

        # Selected cases residual covariance matrix: Cov(X_d,Y|Z,S=1)
        Smat_s_d <- matrix(c(var_X.Z_s_d, var_XY.Z_s_d, var_XY.Z_s_d, var_Y.Z_s), nrow = 2, ncol = 2, byrow = TRUE) # 2x2 matrix

        # Check to see if residual covariance matrix for S=1 will be invertible
        # If not, re-draw the regression coefs from Y|Z,A,S=1 and go back to top of loop (recreate means/variances based on new proxy)
        check <- tryCatch(solve(Smat_s_d), silent = T)
        if (class(check)[1] == "try-error") {
          print(paste("Residual covariance matrix non-invertible (draw = ", d, ")", sep = ""))
          print("Redrawing regression coefficients that create the proxy")
          # Redraw betas in Y|Z,A,S==1 (coefficients that create the proxy)
          DRAWS_beta_YZA.ZA_s[, d] <- rmnorm(mean = beta_YZA.ZA_s, varcov = DRAWS_var_Y.ZA_s[d] * mult.mat.za.s)
        } else if (class(check)[1] == "matrix") {
          smat_invertible <- 1
        }
      }

      ###############
      # DRAWS STEP 2a: Draw residual covariance matrices for (X_d,Y|Z,S)
      ###############
      DRAW_var_XY.Z_s <- riwish(n_s - nZvars - 1, Smat_s_d) * (n_s - nZvars - 1)
      # non-selected: Var(X_d|Z,S=0)
      DRAW_var_X.Z_ns <- (n_ns - nZvars - 1) * var_X.Z_ns_d / rchisq(1, n_ns - nZvars - 1) # scalar

      ###############
      # DRAWS STEP 2b: Draw regression coefficients for X_d,Y|Z,S=1
      ###############
      # Covariance matrix for betas
      coef.vc.s_d <- matrix(0, nrow = (nZvars + 1) * 2, ncol = (nZvars + 1) * 2)
      # upper left (X_d|Z,S=1)
      coef.vc.s_d[1:(nZvars + 1), 1:(nZvars + 1)] <- DRAW_var_XY.Z_s[1, 1] * mult.mat.s
      # lower right (Y|Z,S=1)
      coef.vc.s_d[-c(1:(nZvars + 1)), -c(1:(nZvars + 1))] <- DRAW_var_XY.Z_s[2, 2] * mult.mat.s
      # upper right and lower left
      coef.vc.s_d[1:(nZvars + 1), -c(1:(nZvars + 1))] <- coef.vc.s_d[-c(1:(nZvars + 1)), 1:(nZvars + 1)] <- DRAW_var_XY.Z_s[1, 2] * mult.mat.s
      # Draw of the betas for selected: X_d|Z,S=1, Y|Z,S=1
      DRAW_coef.s <- rmnorm(mean = c(beta_XZ.Z_s_d, beta_YZ.Z_s), varcov = coef.vc.s_d)
      DRAW_beta_X0.Z_s <- DRAW_coef.s[1]
      DRAW_beta_XZ.Z_s <- DRAW_coef.s[2:(nZvars + 1)]
      DRAW_beta_Y0.Z_s <- DRAW_coef.s[nZvars + 2]
      DRAW_beta_YZ.Z_s <- DRAW_coef.s[-c(1:(nZvars + 2))]

      ###############
      # DRAWS STEP 2c: Draw regression coefficients for X_d|Z,S=0
      ###############
      # Covariance matrix for betas
      coef.vc.ns_d <- DRAW_var_X.Z_ns * mult.mat.ns
      # Draw of the betas for non-selected: X|Z,S=0
      DRAW_coef.ns <- rmnorm(mean = beta_XZ.Z_ns_d, varcov = coef.vc.ns_d)
      DRAW_beta_X0.Z_ns <- DRAW_coef.ns[1]
      DRAW_beta_XZ.Z_ns <- DRAW_coef.ns[-1]

      ###############
      # DRAWS STEP 3: Compute draws of other parameters for non-selected (S=1)
      ###############

      # Draw phi from prior specified by phi_character
      if (!is.null(phi_character)) {
        # Draw phi using provided string
        phi <- eval(parse(text = phi_character))
        phi <- pmax(.Machine$double.eps, phi)
      }

      # Cor(X_d,Y|Z,S=1)
      DRAW_rho_XY.Z_s <- DRAW_var_XY.Z_s[1, 2] / sqrt(DRAW_var_XY.Z_s[1, 1] * DRAW_var_XY.Z_s[2, 2])
      # g(phi) multiplier
      g_d <- ((phi + (1 - phi) * DRAW_rho_XY.Z_s) / ((1 - phi) + phi * DRAW_rho_XY.Z_s))
      # ratio of variance components for Y|Z, X|Z
      vratio_d <- DRAW_var_XY.Z_s[2, 2] / DRAW_var_XY.Z_s[1, 1]
      # Regression coefficients for non-selected cases, Y|Z,S=0
      DRAW_beta_Y0.Z_ns <- DRAW_beta_Y0.Z_s + g_d * sqrt(vratio_d) * (DRAW_beta_X0.Z_ns - DRAW_beta_X0.Z_s)
      DRAW_beta_YZ.Z_ns <- DRAW_beta_YZ.Z_s + g_d * sqrt(vratio_d) * (DRAW_beta_XZ.Z_ns - DRAW_beta_XZ.Z_s)
      DRAW_var_Y.Z_ns <- DRAW_var_XY.Z_s[2, 2] + g_d^2 * vratio_d * (DRAW_var_X.Z_ns - DRAW_var_XY.Z_s[1, 1])
      DRAW_covar_XY.Z_ns <- DRAW_var_XY.Z_s[1, 2] + g_d * sqrt(vratio_d) * (DRAW_var_X.Z_ns - DRAW_var_XY.Z_s[1, 1])
      DRAW_beta_YX.XZ_ns <- DRAW_covar_XY.Z_ns / DRAW_var_X.Z_ns
      DRAW_beta_Y0.XZ_ns <- DRAW_beta_Y0.Z_ns - DRAW_beta_YX.XZ_ns * DRAW_beta_X0.Z_ns
      if (DRAW_var_Y.Z_ns - DRAW_beta_YX.XZ_ns^2 * DRAW_var_X.Z_ns > 0) pd <- 1
      #  DRAW_var_Y.Z_ns - DRAW_covar_XY.Z_ns^2/DRAW_var_X.Z_ns   #equivalent check
    }
    ###############
    # DRAWS STEP 4: Compute draws of MUB(phi)
    ###############
    mubdraws[d, 1] <- g_d * sqrt(vratio_d) * (DRAW_beta_X0.Z_s - DRAW_beta_X0.Z_ns)
    mubdraws[d, 2:(nZvars + 1)] <- g_d * sqrt(vratio_d) * (DRAW_beta_XZ.Z_s - DRAW_beta_XZ.Z_ns)
  }

  mubdraws <- as.data.frame(mubdraws)
  # add column names
  for (j in 0:nZvars)
  {
    colnames(mubdraws)[1 + j] <- paste("index.B", j, sep = "")
  }
  return(mubdraws)
}
