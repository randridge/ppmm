#' Calculate Bayes estimates of MUB for binary outcomes
#'
#' @param data_YZA_s list of Y, Z, A for selected sample (Microdata only)
#' @param data_ZA_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param nZvars dimension of predictor variables `Z`
#' @param phi value of phi, fixed to a value between 0 and 1 (defaults to 0)
#' @param phi_character draw phi from specific distribution (i.e."`runif(1)`", defaults to `NULL`)
#' * If `phi_character=NULL`, fix at value set to `phi`;
#' * If `phi_character="discrete"`, draws from `{0,0.5,1}` with 1/3 prob each;
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
#'
#' @seealso [mub_coeff_linear_bayes()]
#' @references
#' West, B. T., Little, R. J., Andridge, R. R., Boonstra, P. S., Ware, E. B., Pandit, A., & Alvarado-Leiton, F. (2021). Assessing selection bias in regression coefficients estimated from nonprobability samples with applications to genetics and demographic surveys. \emph{The annals of applied statistics}, 15(3), 1556–1581.
#'
#' @examples
#' require(mvtnorm)
#' set.seed(123)
#' # num of total sample
#' N <- 10000
#' nZvars <- 2
#' mu <- c(0, 0, 0, 0)
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
#' U <- raw_data[, 1]
#' Z <- raw_data[, 2:3]
#' A <- raw_data[, 4]
#' Y <- ifelse(U > 0, 1, 0)
#' prob <- plogis(-4 + Z %*% c(log(1.1), log(1.1)) + log(1.1) * U + log(2) * A)
#' S <- rbinom(prob, 1, prob)
#' index_s <- S == 1
#' index_ns <- S == 0
#' # generate the list of microdata for selected (data_YZA_s) and non-selected sample (data_ZA_ns)
#' data_YZA_s <- list(Y = Y[index_s], Z = Z[index_s, ], A = A[index_s])
#' data_ZA_ns <- list(Z = Z[index_ns, ], A = A[index_ns])
#'
#' # draw mle result from ppmm model
#' result <- mub_coeff_binary_bayes(data_YZA_s, data_ZA_ns, nZvars)
#'
#' # Use summary statistics for non-selected sample ------------------------
#' # generate the list of summary statistics for non-selected sample (data_ZA_sumry)
#' data_ZA_sumry <- list(
#'   mean_ZA = colMeans(raw_data[index_ns, -1]),
#'   var_ZA = var(raw_data[index_ns, -1]),
#'   n_ZA = sum(index_ns)
#' )
#' result_sumry <- mub_coeff_binary_bayes(data_YZA_s, data_ZA_sumry, nZvars)
#'
mub_coeff_binary_bayes <- function(data_YZA_s,
                                   data_ZA_ns,
                                   nZvars,
                                   phi = 0,
                                   phi_character = NULL,
                                   ndraws = 1500) {
  ### Some functions to do truncated normals in R
  rnorm.lt <- function(n, lv = rep(0, n), mv = rep(0.5, n), sv = rep(1, n)) {
    lstd <- (lv - mv) / sv
    sv * qnorm(runif(n, pnorm(lstd), rep(1, n))) + mv
  }
  rnorm.rt <- function(n, rv = rep(0, n), mv = rep(0.5, n), sv = rep(1, n)) {
    rstd <- (rv - mv) / sv
    sv * qnorm(runif(n, rep(0, n), pnorm(rstd))) + mv
  }

  sumry_list <- Coeff_ErrorCheck(data_YZA_s, data_ZA_ns, prop_check = TRUE)
  data_YZA_s <- sumry_list$data_YZA_s
  sumry_YZA_s <- sumry_list$sumry_YZA_s
  sumry_ZA_ns <- sumry_list$sumry_ZA_ns

  # Some scalars and vectors
  n_s <- length(data_YZA_s$Y)
  Y_s <- data_YZA_s$Y
  A_s <- data_YZA_s$A
  Z_s <- data_YZA_s$Z
  intZA_s <- cbind(rep(1, n_s), Z_s, A_s) # attach intercept for convenience
  intZ_s <- cbind(rep(1, n_s), Z_s)
  ZA_sTZA_sInv <- solve(t(intZA_s) %*% intZA_s) # so only have to calculate once
  mean_A_ns <- sumry_ZA_ns$mean_ZA[-c(1:nZvars)]
  var_A_ns <- sumry_ZA_ns$var_ZA[-c(1:nZvars), -c(1:nZvars)]
  var_ZA_ns <- sumry_ZA_ns$var_ZA

  ##################
  # pieces involving Z needed for draws
  var_Z_s <- var(Z_s)
  # mean_Z_s <- colMeans(Z_s)
  if (nZvars == 1) {
    mean_Z_s <- mean(Z_s)
  } else {
    mean_Z_s <- colMeans(Z_s)
  }
  mean_Z_ns <- sumry_ZA_ns$mean_ZA[1:nZvars]
  var_Z_ns <- sumry_ZA_ns$var_ZA[1:nZvars, 1:nZvars]
  n_ns <- sumry_ZA_ns$n_ZA
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

  #### Regression of Y|Z,A,S=1 --> to create latent U and to create proxy X
  fit.y.za_s <- glm(data_YZA_s$Y ~ data_YZA_s$Z + data_YZA_s$A, family = binomial(link = "probit"))
  B <- as.matrix(fit.y.za_s$coef)
  rownames(B) <- NULL

  # Starting proxy value, X = A*Ba -- selected sample
  X_s <- as.matrix(A_s) %*% B[-c(1:(1 + nZvars))]

  # Matrix to hold draws of MUB(phi)
  mubdraws <- matrix(0, nrow = ndraws, ncol = (nZvars + 1))

  # Intialize latent U vector -- selected sample
  U_s <- rep(NA, n_s)

  # Looping to obtain draws
  for (d in 1:ndraws)
  {
    ## (1) Draw latent U|Y,B,S=1  (aka U|Y,Z,A) -- selected sample
    linpred <- intZ_s %*% B[1:(1 + nZvars)] + X_s
    U_s[Y_s == 1] <- rnorm.lt(sum(Y_s), mv = linpred[Y_s == 1]) # Draws for Y=1 --> truncated at left by 0
    U_s[Y_s == 0] <- rnorm.rt(sum(1 - Y_s), mv = linpred[Y_s == 0]) # Draws for Y=0 --> truncated at right by 0

    ## (2) Draw B|Z,A,U
    B <- rmnorm(1, ZA_sTZA_sInv %*% t(intZA_s) %*% U_s, ZA_sTZA_sInv)

    ## (3) Create proxy given current B -- selected sample
    X_s <- as.matrix(A_s) %*% B[-c(1:(1 + nZvars))]

    #### STARTING FROM HERE, CODE BELOW ADAPTED FROM mub_reg_bayes.R with U in place of Y

    #### Regression of U|Z,S=1 (target regression of interest, using latent U in place of Y)
    mean_UZ_s <- colMeans(cbind(U_s, Z_s))
    var_UZ_s <- var(cbind(U_s, Z_s))
    # MLEs of regression coefficients for U|Z,S=1
    beta_UZ.Z_s <- drop(var_UZ_s[1, -1] %*% solve(var_UZ_s[-1, -1])) # slopes
    beta_UZ.Z_s <- c(mean_UZ_s[1] - beta_UZ.Z_s %*% mean_UZ_s[-1], beta_UZ.Z_s) # intercept
    # Residual variance
    var_U.Z_s <- (n_s / (n_s - nZvars)) * drop(var_UZ_s[1, 1] - var_UZ_s[1, -1] %*% solve(var_UZ_s[-1, -1]) %*% var_UZ_s[-1, 1])

    pd <- 0
    while (pd == 0) # Checking draws to ensure positive-definite covariance matrix
    {
      smat_invertible <- 0
      while (smat_invertible == 0) # Checking to ensure residual covariance matrix (X,U|Z,S=1) is invertible
      {
        #### Means and variances for (X,Z|S)
        # Selected
        if (nZvars == 1) {
          mean_XZ_s_d <- c(mean(X_s), mean(Z_s))
        } else {
          mean_XZ_s_d <- c(mean(X_s), colMeans(Z_s))
        }
        var_XZ_s_d <- var(cbind(X_s, Z_s))
        # Non-selected
        mean_XZ_ns_d <- c(mean_A_ns %*% B[-c(1:(nZvars + 1))], mean_Z_ns)
        var_XZ_ns_d <- matrix(nrow = 1 + nZvars, ncol = 1 + nZvars)
        var_XZ_ns_d[1, 1] <- drop(crossprod(B[-c(1:(nZvars + 1))], var_A_ns) %*% B[-c(1:(nZvars + 1))])
        var_XZ_ns_d[1, -1] <- var_XZ_ns_d[-1, 1] <- drop(crossprod(B[-c(1:(nZvars + 1))], var_ZA_ns[-c(1:nZvars), 1:nZvars]))
        var_XZ_ns_d[-1, -1] <- var_ZA_ns[1:nZvars, 1:nZvars]

        #### Regression of X|Z,S=1
        beta_XZ.Z_s_d <- drop(var_XZ_s_d[1, -1] %*% solve(var_XZ_s_d[-1, -1])) # slopes
        beta_XZ.Z_s_d <- c(mean_XZ_s_d[1] - beta_XZ.Z_s_d %*% mean_XZ_s_d[-1], beta_XZ.Z_s_d) # intercept
        # Residual variance
        var_X.Z_s_d <- (n_s / (n_s - nZvars)) * drop(var_XZ_s_d[1, 1] - var_XZ_s_d[1, -1] %*% solve(var_XZ_s_d[-1, -1]) %*% var_XZ_s_d[-1, 1])

        #### Residual covariance of X,U|Z,S=1
        var_UA_s <- var(cbind(U_s, A_s))
        var_UZ_s <- var(cbind(U_s, Z_s))
        var_XU.Z_s_d <- drop(B[-c(1:(nZvars + 1))] %*% var_UA_s[1, -1] - var_XZ_s_d[1, -1] %*% solve(var_XZ_s_d[-1, -1]) %*% var_UZ_s[1, -1])

        #### Regression of X|Z,S=0
        beta_XZ.Z_ns_d <- drop(var_XZ_ns_d[1, -1] %*% solve(var_XZ_ns_d[-1, -1])) # slopes
        beta_XZ.Z_ns_d <- c(mean_XZ_ns_d[1] - beta_XZ.Z_ns_d %*% mean_XZ_ns_d[-1], beta_XZ.Z_ns_d) # intercept
        # Residual variance
        var_X.Z_ns_d <- (n_ns / (n_ns - nZvars)) * drop(var_XZ_ns_d[1, 1] - var_XZ_ns_d[1, -1] %*% solve(var_XZ_ns_d[-1, -1]) %*% var_XZ_ns_d[-1, 1])

        # Selected cases residual covariance matrix: Cov(X_d,U|Z,S=1)
        Smat_s_d <- matrix(c(var_X.Z_s_d, var_XU.Z_s_d, var_XU.Z_s_d, var_U.Z_s), nrow = 2, ncol = 2, byrow = TRUE) # 2x2 matrix

        # Check to see if residual covariance matrix for S=1 will be invertible
        # If not, re-draw the regression coefs from Y|Z,A,S=1 and go back to top of loop (recreate means/variances based on new proxy)
        check <- tryCatch(solve(Smat_s_d), silent = T)
        if (class(check)[1] == "try-error") {
          print(paste("Residual covariance matrix non-invertible (draw = ", d, ")", sep = ""))
          print("Redrawing regression coefficients that create the proxy")
          # Redraw betas in Y|Z,A,S==1 (coefficients that create the proxy)
          B <- t(rmnorm(1, ZA_sTZA_sInv %*% t(intZA_s) %*% U_s, ZA_sTZA_sInv))
        } else if (class(check)[1] == "matrix") {
          smat_invertible <- 1
        }
      }

      ###############
      # DRAWS STEP 2a: Draw residual covariance matrices for (X_d,U|Z,S)
      ###############
      DRAW_var_XU.Z_s <- riwish(n_s - nZvars - 1, Smat_s_d) * (n_s - nZvars - 1)
      # non-selected: Var(X_d|Z,S=0)
      DRAW_var_X.Z_ns <- (n_ns - nZvars - 1) * var_X.Z_ns_d / rchisq(1, n_ns - nZvars - 1) # scalar

      ###############
      # DRAWS STEP 2b: Draw regression coefficients for X_d,U|Z,S=1
      ###############
      # Covariance matrix for betas
      coef.vc.s_d <- matrix(0, nrow = (nZvars + 1) * 2, ncol = (nZvars + 1) * 2)
      # upper left (X_d|Z,S=1)
      coef.vc.s_d[1:(nZvars + 1), 1:(nZvars + 1)] <- DRAW_var_XU.Z_s[1, 1] * mult.mat.s
      # lower right (Y|Z,S=1)
      coef.vc.s_d[-c(1:(nZvars + 1)), -c(1:(nZvars + 1))] <- DRAW_var_XU.Z_s[2, 2] * mult.mat.s
      # upper right and lower left
      coef.vc.s_d[1:(nZvars + 1), -c(1:(nZvars + 1))] <- coef.vc.s_d[-c(1:(nZvars + 1)), 1:(nZvars + 1)] <- DRAW_var_XU.Z_s[1, 2] * mult.mat.s
      # Draw of the betas for selected: X_d|Z,S=1, Y|Z,S=1
      DRAW_coef.s <- rmnorm(mean = c(beta_XZ.Z_s_d, beta_UZ.Z_s), varcov = coef.vc.s_d)
      DRAW_beta_X0.Z_s <- DRAW_coef.s[1]
      DRAW_beta_XZ.Z_s <- DRAW_coef.s[2:(nZvars + 1)]
      DRAW_beta_U0.Z_s <- DRAW_coef.s[nZvars + 2]
      DRAW_beta_UZ.Z_s <- DRAW_coef.s[-c(1:(nZvars + 2))]

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
        # "discrete" for draws from {0,0.5,1} with 1/3 prob each
        if (phi_character == "discrete") {
          phi <- sample(c(0, 0.5, 1), size = 1, prob = c(1 / 3, 1 / 3, 1 / 3))
        } else {
          # Draw phi using provided string
          phi <- eval(parse(text = phi_character))
          phi <- pmax(.Machine$double.eps, phi)
        }
      }

      # Cor(X_d,U|Z,S=1)
      DRAW_rho_XU.Z_s <- DRAW_var_XU.Z_s[1, 2] / sqrt(DRAW_var_XU.Z_s[1, 1] * DRAW_var_XU.Z_s[2, 2])
      # g(phi) multiplier
      g_d <- ((phi + (1 - phi) * DRAW_rho_XU.Z_s) / ((1 - phi) + phi * DRAW_rho_XU.Z_s))
      # ratio of variance components for U|Z, X|Z
      vratio_d <- DRAW_var_XU.Z_s[2, 2] / DRAW_var_XU.Z_s[1, 1]
      # Regression coefficients for non-selected cases, U|Z,S=0
      DRAW_beta_U0.Z_ns <- DRAW_beta_U0.Z_s + g_d * sqrt(vratio_d) * (DRAW_beta_X0.Z_ns - DRAW_beta_X0.Z_s)
      DRAW_beta_UZ.Z_ns <- DRAW_beta_UZ.Z_s + g_d * sqrt(vratio_d) * (DRAW_beta_XZ.Z_ns - DRAW_beta_XZ.Z_s)
      DRAW_var_U.Z_ns <- DRAW_var_XU.Z_s[2, 2] + g_d^2 * vratio_d * (DRAW_var_X.Z_ns - DRAW_var_XU.Z_s[1, 1])
      DRAW_covar_XU.Z_ns <- DRAW_var_XU.Z_s[1, 2] + g_d * sqrt(vratio_d) * (DRAW_var_X.Z_ns - DRAW_var_XU.Z_s[1, 1])
      DRAW_beta_UX.XZ_ns <- DRAW_covar_XU.Z_ns / DRAW_var_X.Z_ns
      DRAW_beta_U0.XZ_ns <- DRAW_beta_U0.Z_ns - DRAW_beta_UX.XZ_ns * DRAW_beta_X0.Z_ns
      if (DRAW_var_U.Z_ns - DRAW_beta_UX.XZ_ns^2 * DRAW_var_X.Z_ns > 0) pd <- 1
    }
    ###############
    # DRAWS STEP 4: Compute draws of MUB(phi) - selected minus non-selected coefficients
    ###############
    mubdraws[d, 1] <- DRAW_beta_U0.Z_s / sqrt(DRAW_var_XU.Z_s[2, 2]) - DRAW_beta_U0.Z_ns / sqrt(DRAW_var_U.Z_ns)
    mubdraws[d, 2:(nZvars + 1)] <- DRAW_beta_UZ.Z_s / sqrt(DRAW_var_XU.Z_s[2, 2]) - DRAW_beta_UZ.Z_ns / sqrt(DRAW_var_U.Z_ns)
  }
  mubdraws <- as.data.frame(mubdraws)
  # add column names
  for (j in 0:nZvars)
  {
    colnames(mubdraws)[1 + j] <- paste("index.B", j, sep = "")
  }
  return(mubdraws)
}
