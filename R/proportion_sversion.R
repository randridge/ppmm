#' Calculate Bayes estimates for binary outcomes
#'
#' @param data_YZ_s list for selected sample (Microdata). See details below for more details.
#' @param data_Z_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param phi value of phi, fixed to a value between 0 and 1 (defaults to 0)
#' @param phi_character draw phi from specific distribution (i.e."`runif(1)`", defaults to `NULL`)
#' * If `phi_character=NULL`, fix at value set to `phi`;
#' * If `phi_character=` specific distribution (i.e. `runif(1)` or `rbeta(1,1,1)`), then the value of `phi` input to the function is ignored
#' @param ndraws number of draws
#'
#' @returns A `ndraws*17` matrix of Bayes estimates with each row containing:
#' \itemize{
#' \item `phi` : drawn phi in this draw
#' \item `pi`: drawn proportion of selected sample
#' \item `muX_s` : drawn mean of the proxy for selected sample
#' \item `muU_s` : drawn mean of the underlying latent variable for selected sample
#' \item `muY_s` : drawn mean of outcome variable for selected sample
#' \item `muX_ns` : drawn mean of the proxy for non-selected sample
#' \item `muU_ns` : drawn mean of the underlying latent variable for non-selected sample
#' \item `muY_ns` : drawn mean of outcome variable for non-selected sample
#' \item `sigmaXX_s` : drawn variance of the proxy for selected sample
#' \item `sigmaXU_s` : drawn covariance between the proxy and the underlying latent variable for selected sample
#' \item `sigmaUU_s` : drawn variance of the underlying latent variable for selected sample
#' \item `sigmaXX_ns` : drawn variance of the proxy for non-selected sample
#' \item `sigmaXU_ns` : drawn covariance between the proxy and the underlying latent variable for non-selected sample
#' \item `sigmaUU_ns` : drawn variance of the underlying latent variable for non-selected sample
#' \item `rho_s` : drawn biserial correlation for selected sample
#' \item `muY` : drawn population mean of outcome variable
#' \item `msb` : drawn MSB
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
#' @seealso [means_bayes()]
#' @seealso [means_mle()]
#' @seealso [prop_mle()]
#'
#' @references
#' Andridge, R. R., West, B. T., Little, R. J. A., Boonstra, P. S., & Alvarado-Leiton, F. (2019). Indices of non-ignorable selection bias for proportions estimated from non-probability samples. \emph{Journal of the Royal Statistical Society. Series C, Applied statistics, 68(5)}, 1465–1483.
#' @examples
#' set.seed(123)
#' # num of total sample
#' N <- 10000
#' # generate auxiliary variables Z
#' Z <- rnorm(N, 0, 1)
#' # define the biserial correlation rho between outcome variable Y and the proxy X
#' rho <- 0.5
#' # define the mean of the outcome variable Y
#' mu_y <- 0.3
#' # define a0, a1 such that E[Y] = mu_y
#' a1 <- rho / sqrt(1 - rho^2)
#' a0 <- qnorm(mu_y) * sqrt(1 + a1^2)
#' # generate the latent variable U
#' U <- rnorm(Z, a0 + a1 * Z, 1)
#' # generate the observed variable Y
#' Y <- ifelse(U > 0, 1, 0)
#' # generate the sample selection indicator S with 0.05 of the population selected
#' prob <- plogis(-3 + 0.3 * Z + 0.1 * U)
#' S <- rbinom(prob, 1, prob)
#' index_s <- S == 1
#' index_ns <- S == 0
#' # generate the list of microdata for selected (data_YZ_s) and non-selected sample (data_Z_ns)
#' data_YZ_s <- list(Y = Y[index_s], Z = Z[index_s])
#' data_Z_ns <- list(Z = Z[index_ns])
#'
#' # draw bayesian result from ppmm model
#' result <- prop_bayes(data_YZ_s, data_Z_ns, phi_character = "runif(1)")
#'
#' # Use summary statistics for non-selected sample ------------------------------------------
#' # generate the list of summary statistics for non-selected sample (data_Z_sumry)
#' data_Z_sumry <- list(
#'   mean_Z = mean(Z[index_ns]),
#'   var_Z = var(Z[index_ns]),
#'   n_Z = length(Z[index_ns])
#' )
#' result_sumry <- prop_bayes(data_YZ_s, data_Z_sumry, phi_character = "runif(1)")
#'
prop_bayes <- function(data_YZ_s,
                       data_Z_ns,
                       phi = 0,
                       phi_character = NULL,
                       ndraws = 1000) {
  ### Some functions to do truncated normals in R
  rnorm.lt <- function(n, lv = rep(0, n), mv = rep(0.5, n), sv = rep(1, n)) {
    lstd <- (lv - mv) / sv
    sv * qnorm(runif(n, pnorm(lstd), rep(1, n))) + mv
  }
  rnorm.rt <- function(n, rv = rep(0, n), mv = rep(0.5, n), sv = rep(1, n)) {
    rstd <- (rv - mv) / sv
    sv * qnorm(runif(n, rep(0, n), pnorm(rstd))) + mv
  }

  sumry_list <- ErrorCheck(data_YZ_s, data_Z_ns, prop_check = TRUE)
  data_YZ_s <- sumry_list$data_YZ_s
  sumry_Z_ns <- sumry_list$sumry_Z_ns

  # Make sure sumry_Z_ns$mean_Z (mean of Z for non-selected sample) is a vector (in case it is a scalar)
  if (!is.vector(sumry_Z_ns$mean_Z)) sumry_Z_ns$mean_Z <- as.matrix(sumry_Z_ns$mean_Z)

  # Total n
  n_s <- length(data_YZ_s$Y) # sample size for selected/non-missing sample
  n_ns <- sumry_Z_ns$n_Z # sample size for non-selected/missing sample
  n <- n_s + n_ns

  # Number of selected cases where Y=1 and where Y=0
  n_s_y1 <- sum(data_YZ_s$Y)
  n_s_y0 <- n_s - n_s_y1

  # Probit regression of Y|Z for selected cases to find starting point for Gibbs sampler
  fit <- glm(data_YZ_s$Y ~ data_YZ_s$Z, family = binomial(link = "probit"))
  betaHat <- as.vector(fit$coef)
  z_s <- model.matrix(fit) # attaches column of 1s (Selected sample)
  zmean_ns <- c(1, sumry_Z_ns$mean_Z) # attaches 1 to mean vector (Non-selected sample)
  z_sTz_sInv <- solve(t(z_s) %*% z_s) # so only have to calculate once

  # Starting point from probit fit
  B <- as.matrix(betaHat)

  # Starting proxy value
  X_s <- z_s %*% B # Selected sample (microdata)
  Xmean_ns <- zmean_ns %*% B # Non-selected sample (sum stats only)
  Xvar_ns <- t(B[-1]) %*% sumry_Z_ns$var_Z %*% B[-1] # Non-selected sample (sum stats only)

  # Initialize matrix to hold results
  draws <- matrix(nrow = ndraws, ncol = 17)

  # Start looping
  for (j in 1:ndraws)
  {
    # Intialize latent U vector (Selected sample only)
    u_s <- rep(NA, n_s)

    ## (0) Draw phi from prior specified by phi_character
    if (!is.null(phi_character)) {
      # Draw phi using provided string
      phi <- eval(parse(text = phi_character))
      phi <- pmax(.Machine$double.eps, phi)
    }

    ## (1) Draw latent U|Y,B for selected sample (aka U|Y,X)
    u_s[data_YZ_s$Y == 1] <- rnorm.lt(n_s_y1, mv = X_s[data_YZ_s$Y == 1]) # Draws for Y=1 --> truncated at left by 0
    u_s[data_YZ_s$Y == 0] <- rnorm.rt(n_s_y0, mv = X_s[data_YZ_s$Y == 0]) # Draws for Y=0 --> truncated at right by 0

    ## (2) Draw B|Z,U
    B <- t(rmvnorm(1, z_sTz_sInv %*% t(z_s) %*% u_s, z_sTz_sInv))

    ## (3) Create proxy given current B
    X_s <- z_s %*% B # Selected sample (microdata)
    Xmean_ns <- zmean_ns %*% B # Non-selected sample (sum stats only)
    Xvar_ns <- t(B[-1]) %*% sumry_Z_ns$var_Z %*% B[-1] # Non-selected sample (sum stats only)

    ## (4) Scale proxy X to have same variance as latent U among respondents
    # Draw the population variances of X, U from posterior
    varXdraw <- sum((X_s - mean(X_s))^2) / rchisq(1, n_s - 1)
    varUdraw <- sum((u_s - mean(u_s))^2) / rchisq(1, n_s - 1)
    # Use draws to scale the proxy
    x_s <- X_s * sqrt(varUdraw / varXdraw) # Selected sample (microdata)
    xmean_ns <- Xmean_ns * sqrt(varUdraw / varXdraw) # Non-selected sample (sum stats only)
    xvar_ns <- Xvar_ns * varUdraw / varXdraw # Non-selected sample (sum stats only)


    ## (5) Draw from PPM dependent on value of phi, using (X,U)
    if (phi == 0) {
      y1_s <- x_s # Selected sample (microdata)
      y2_s <- u_s # Selected sample (microdata)
      y1_ns_mean <- xmean_ns # Non-selected sample (sum stats only)
      y1_ns_var <- xvar_ns # Non-selected sample (sum stats only)

      # Calculate means, sums of squares
      # Selected sample
      y1Bar_s <- mean(y1_s)
      y2Bar_s <- mean(y2_s)
      s11_s <- sum((y1_s - y1Bar_s)^2) / n_s
      s22_s <- sum((y2_s - y2Bar_s)^2) / n_s
      s12_s <- sum((y1_s - y1Bar_s) * (y2_s - y2Bar_s)) / n_s
      b21.1_s <- s12_s / s11_s
      s22.1_s <- s22_s - (s12_s)^2 / s11_s
      # Nonrespondent data
      y1Bar_ns <- y1_ns_mean
      s11_ns <- y1_ns_var
      # Draws
      # (1) PI
      PI <- rbeta(1, n_s + 0.5, n_ns + 0.5)
      # (2) SIGMA11_s
      SIGMA11_s <- n_s * s11_s / rchisq(1, n_s - 1)
      # (3) MU1_s | SIGMA11_s
      MU1_s <- rnorm(1, y1Bar_s, sqrt(SIGMA11_s / n_s))
      # (4) SIGMA22.1_s
      SIGMA22.1_s <- n_s * s22.1_s / rchisq(1, n_s - 2)
      # (5) BETA21.1_s | SIGMA22.1_s
      BETA21.1_s <- rnorm(1, b21.1_s, sqrt(SIGMA22.1_s / (n_s * s11_s)))
      # (6) BETA20.1_s | BETA21.1_s, SIGMA22.1_s
      BETA20.1_s <- rnorm(1, y2Bar_s - BETA21.1_s * y1Bar_s, sqrt(SIGMA22.1_s / n_s))
      # (7) SIGMA11_ns
      SIGMA11_ns <- n_ns * s11_ns / rchisq(1, n_ns - 1)
      # (8) MU1_ns | SIGMA11_ns
      MU1_ns <- rnorm(1, y1Bar_ns, sqrt(SIGMA11_ns / n_ns))
      # Transform draws to get other parameters
      # (a) MU2_s
      MU2_s <- BETA20.1_s + BETA21.1_s * MU1_s
      # (b) MU2_ns
      MU2_ns <- BETA20.1_s + BETA21.1_s * MU1_ns
      # (c) SIGMA12_s
      SIGMA12_s <- BETA21.1_s * SIGMA11_s
      # (d) SIGMA22_s
      SIGMA22_s <- SIGMA22.1_s + BETA21.1_s^2 * SIGMA11_s
      # (e) SIGMA22_ns
      SIGMA22_ns <- SIGMA22.1_s + BETA21.1_s^2 * SIGMA11_ns
      # (f) SIGMA12_ns
      SIGMA12_ns <- SIGMA11_ns * BETA21.1_s
      # All Draws
      drawsPPM <- list(
        pi = PI, muX_s = MU1_s, muU_s = MU2_s, muX_ns = MU1_ns, muU_ns = MU2_ns,
        sigmaXX_s = SIGMA11_s, sigmaXU_s = SIGMA12_s, sigmaUU_s = SIGMA22_s,
        sigmaXX_ns = SIGMA11_ns, sigmaXU_ns = SIGMA12_ns, sigmaUU_ns = SIGMA22_ns
      )
    } else {
      if (phi == 1) {
        y1_s <- x_s # Selected sample (microdata)
        y2_s <- u_s # Selected sample (microdata)
        y1_ns_mean <- xmean_ns # Non-selected sample (sum stats only)
        y1_ns_var <- xvar_ns # Non-selected sample (sum stats only)
      } else {
        y1_s <- x_s # Selected sample (microdata)
        y2_s <- (1 - phi) * x_s + phi * u_s # Selected sample (microdata)
        y1_ns_mean <- xmean_ns # Non-selected sample (sum stats only)
        y1_ns_var <- xvar_ns # Non-selected sample (sum stats only)
      }
      # Calculate means, sums of squares
      # Selected sample
      y1Bar_s <- mean(y1_s)
      y2Bar_s <- mean(y2_s)
      s11_s <- sum((y1_s - y1Bar_s)^2) / n_s
      s22_s <- sum((y2_s - y2Bar_s)^2) / n_s
      s12_s <- sum((y1_s - y1Bar_s) * (y2_s - y2Bar_s)) / n_s
      b12.2_s <- s12_s / s22_s
      s11.2_s <- s11_s - (s12_s)^2 / s22_s
      # Nonrespondent data
      y1Bar_ns <- y1_ns_mean
      s11_ns <- y1_ns_var
      # Draws
      # (1) PI
      PI <- rbeta(1, n_s + 0.5, n_ns + 0.5)
      # (2) SIGMA22_s
      SIGMA22_s <- n_s * s22_s / rchisq(1, n_s - 1)
      # (3) MU2_s | SIGMA22_s
      MU2_s <- rnorm(1, y2Bar_s, sqrt(SIGMA22_s / n_s))
      # (4, 5) SIGMA11.2_s, SIGMA11_ns with constraint
      goodDraw <- FALSE
      ct <- 1
      while (!goodDraw) {
        # Repeat these 2 draws until SIGMA11_ns > SIGMA11.2_s
        # (4) SIGMA11.2_s
        SIGMA11.2_s <- n_s * s11.2_s / rchisq(1, n_s - 2)
        # (5) SIGMA11_ns
        SIGMA11_ns <- n_ns * s11_ns / rchisq(1, n_ns - 1)
        # Check to see if draws meet the condition
        goodDraw <- (SIGMA11_ns >= SIGMA11.2_s)
        if (ct > 20) {
          goodDraw <- TRUE
          SIGMA11.2_s <- SIGMA11_ns
        }
        ct <- ct + 1
      }
      # (6) BETA12.2_s | SIGMA11.2_s
      BETA12.2_s <- rnorm(1, b12.2_s, sqrt(SIGMA11.2_s / (n_s * s22_s)))
      # (7) BETA10.2_s | BETA12.2_s, SIGMA11.2_s
      BETA10.2_s <- rnorm(1, y1Bar_s - BETA12.2_s * y2Bar_s, sqrt(SIGMA11.2_s / n_s))
      # (8) MU1_ns | SIGMA11_ns
      MU1_ns <- rnorm(1, y1Bar_ns, sqrt(SIGMA11_ns / n_ns))
      # Transform draws to get other parameters
      # (a) MU2_ns
      MU2_ns <- (MU1_ns - BETA10.2_s) / BETA12.2_s
      # (b) MU1_s
      MU1_s <- BETA10.2_s + BETA12.2_s * MU2_s
      # (c) SIGMA12_s
      SIGMA12_s <- BETA12.2_s * SIGMA22_s
      # (d) SIGMA11_s
      SIGMA11_s <- SIGMA11.2_s + BETA12.2_s^2 * SIGMA22_s
      # (e) SIGMA22_ns
      SIGMA22_ns <- (SIGMA11_ns - SIGMA11.2_s) / BETA12.2_s^2
      # (f) SIGMA12_ns
      SIGMA12_ns <- SIGMA22_ns * BETA12.2_s
      # All Draws
      drawsPPM <- list(
        pi = PI, muX_s = MU1_s, muU_s = MU2_s, muX_ns = MU1_ns, muU_ns = MU2_ns,
        sigmaXX_s = SIGMA11_s, sigmaXU_s = SIGMA12_s, sigmaUU_s = SIGMA22_s,
        sigmaXX_ns = SIGMA11_ns, sigmaXU_ns = SIGMA12_ns, sigmaUU_ns = SIGMA22_ns
      )
    }
    if (phi != 0 & phi != 1) {
      # Transform draws of [X,W] to get draws from [X,U]
      # W = (1-phi)*X + phi*U --> U = (W - (1-phi)*X)/phi
      # Start with draws of [X,W] and then overwrite parms relating to U
      drawsXW <- drawsPPM
      drawsPPM$muU_s <- (drawsXW$muU_s - (1 - phi) * drawsXW$muX_s) / phi
      drawsPPM$muU_ns <- (drawsXW$muU_ns - (1 - phi) * drawsXW$muX_ns) / phi
      drawsPPM$sigmaUU_s <- (drawsXW$sigmaUU_s + (1 - phi)^2 * drawsXW$sigmaXX_s - 2 * (1 - phi) * drawsXW$sigmaXU_s) / phi^2
      drawsPPM$sigmaUU_ns <- (drawsXW$sigmaUU_ns + (1 - phi)^2 * drawsXW$sigmaXX_ns - 2 * (1 - phi) * drawsXW$sigmaXU_ns) / phi^2
      drawsPPM$sigmaXU_s <- (drawsXW$sigmaXU_s - (1 - phi) * drawsXW$sigmaXX_s) / phi
      drawsPPM$sigmaXU_ns <- (drawsXW$sigmaXU_ns - (1 - phi) * drawsXW$sigmaXX_ns) / phi
    }

    ## (8) Transform parameters to get draws of the mean of Y (unmodified Bayes method)
    #### Estimate of mean of Y for selected sample
    drawsPPM$muY_s <- pnorm(drawsPPM$muU_s / sqrt(drawsPPM$sigmaUU_s))
    #### Estimate of mean of Y for non-selected sample
    drawsPPM$muY_ns <- pnorm(drawsPPM$muU_ns / sqrt(drawsPPM$sigmaUU_ns))
    #### Estimate of mean of Y (overall)
    drawsPPM$muY <- drawsPPM$pi * drawsPPM$muY_s + (1 - drawsPPM$pi) * drawsPPM$muY_ns

    #### Estimate of MUBP index
    drawsPPM$msb <- drawsPPM$muY_s - drawsPPM$muY

    #### Value of phi (important to keep if drawing phi)
    drawsPPM$phi <- phi

    #### Calculate the biserial correlation for selected sample
    drawsPPM$rho_s <- drawsPPM$sigmaXU_s / sqrt(drawsPPM$sigmaXX_s * drawsPPM$sigmaUU_s)

    # Save draws
    draws[j, ] <- unlist(drawsPPM)
  }
  # End looping
  draws <- as.data.frame(draws)
  names(draws) <- names(unlist(drawsPPM))
  output_order <- c(
    "phi", "pi",
    "muX_s", "muU_s", "muY_s",
    "muX_ns", "muU_ns", "muY_ns",
    "sigmaXX_s", "sigmaXU_s", "sigmaUU_s",
    "sigmaXX_ns", "sigmaXU_ns", "sigmaUU_ns",
    "rho_s",
    "muY", "msb"
  )
  return(draws[, output_order])
}


#' Calculate MUBP estimates for binary outcomes with maximum likelihood estimation
#'
#' @param data_YZ_s list for selected sample (Microdata). See details below for more details.
#' @param data_Z_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param phi value of phi with a value between 0 and 1
#' @param sfrac if sfrac is not NULL, the proportion of selected sample is fixed to `sfrac`, else it is estimated from the data (the proportion of selected sample).
#'
#' @returns A `length(phi)*16` matrix of Bayes estimates with each row containing:
#' \itemize{
#' \item `phi` : corresponding phi
#' \item `pi`: estimated proportion of selected sample (or defined by user)
#' \item `muX_s` : ML estimate of the mean of the proxy for selected sample
#' \item `muU_s` : ML estimate of the mean of the underlying latent variable for selected sample
#' \item `muY_s` : ML estimate of the mean of outcome variable for selected sample
#' \item `muX_ns` : ML estimate of the mean of the proxy for non-selected sample
#' \item `muU_ns` : ML estimate of the mean of the underlying latent variable for non-selected sample
#' \item `muY_ns` : ML estimate of the mean of outcome variable for non-selected sample
#' \item `sigmaXX_s` : ML estimate of the variance of the proxy for selected sample
#' \item `sigmaXU_s` : ML estimate of the covariance between the proxy and the underlying latent variable for selected sample
#' \item `sigmaUU_s` : ML estimate of the variance of the underlying latent variable for selected sample
#' \item `sigmaXX_ns` : ML estimate of the variance of the proxy for non-selected sample
#' \item `sigmaUU_ns` : ML estimate of the variance of the underlying latent variable for non-selected sample
#' \item `rho_s` : ML estimate of the biserial correlation for selected sample
#' \item `muY` : ML estimate of the population mean of outcome variable
#' \item `msb` : ML estimate of the MSB
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
#' @seealso [means_mle()]
#' @seealso [prop_bayes()]
#'
#' @references
#' Andridge, R. R., West, B. T., Little, R. J. A., Boonstra, P. S., & Alvarado-Leiton, F. (2019). Indices of non-ignorable selection bias for proportions estimated from non-probability samples. \emph{Journal of the Royal Statistical Society. Series C, Applied statistics, 68(5)}, 1465–1483.
#' @examples
#' set.seed(123)
#' # num of total sample
#' N <- 10000
#' # generate auxiliary variables Z
#' Z <- rnorm(N, 0, 1)
#' # define the biserial correlation rho between outcome variable Y and the proxy X
#' rho <- 0.5
#' # define the mean of the outcome variable Y
#' mu_y <- 0.3
#' # define a0, a1 such that E[Y] = mu_y
#' a1 <- rho / sqrt(1 - rho^2)
#' a0 <- qnorm(mu_y) * sqrt(1 + a1^2)
#' # generate the latent variable U
#' U <- rnorm(Z, a0 + a1 * Z, 1)
#' # generate the observed variable Y
#' Y <- ifelse(U > 0, 1, 0)
#' # generate the sample selection indicator S with 0.05 of the population selected
#' prob <- plogis(-3 + 0.3 * Z + 0.1 * U)
#' S <- rbinom(prob, 1, prob)
#' index_s <- S == 1
#' index_ns <- S == 0
#' # generate the list of microdata for selected (data_YZ_s) and non-selected sample (data_Z_ns)
#' data_YZ_s <- list(Y = Y[index_s], Z = Z[index_s])
#' data_Z_ns <- list(Z = Z[index_ns])
#'
#' ## calculate MUB result from ppmm model with maximum likelihood estimation
#' result <- prop_mle(data_YZ_s, data_Z_ns)
#'
#' # Use summary statistics for non-selected sample ------------------------------------------
#' # generate the list of summary statistics for non-selected sample (data_Z_sumry)
#' data_Z_sumry <- list(
#'   mean_Z = mean(Z[index_ns]),
#'   var_Z = var(Z[index_ns]),
#'   n_Z = length(Z[index_ns])
#' )
#' result_sumry <- prop_mle(data_YZ_s, data_Z_sumry)
#'
prop_mle <- function(data_YZ_s,
                     data_Z_ns,
                     phi = c(0, 0.5, 1),
                     sfrac = NULL) {
  sumry_list <- ErrorCheck(data_YZ_s, data_Z_ns, prop_check = TRUE)
  data_YZ_s <- sumry_list$data_YZ_s
  sumry_Z_ns <- sumry_list$sumry_Z_ns

  # Make sure sumry_Z_ns$mean_Z (mean of Z for non-selected sample) is a vector (in case it is a scalar)
  if (!is.vector(sumry_Z_ns$mean_Z)) sumry_Z_ns$mean_Z <- as.matrix(sumry_Z_ns$mean_Z)

  # Probit regression of Y|Z for selected cases to find starting point for Gibbs sampler
  fit <- glm(data_YZ_s$Y ~ data_YZ_s$Z, family = binomial(link = "probit"))
  # Summary stats for the proxy X for the non-selected sample (mMLEs)
  # Coefficients from probit fit
  B <- as.vector(coef(fit))
  # Mean of proxy in non-selected sample (using population mean vector of Z)
  Xmean_ns <- drop(c(1, sumry_Z_ns$mean_Z) %*% B) # Non-selected sample (sum stats only)
  # Calculate variance in non-selected sample (using population covariance matrix of Z)
  Xvar_ns <- drop(t(B[-1]) %*% sumry_Z_ns$var_Z %*% B[-1]) # Non-selected sample (sum stats only)
  # Create proxy for selected sample (linear predictor from Probit model)
  X_s <- predict(fit, type = "link")
  # outcome vector Y for selected sample
  y_s <- data_YZ_s$Y

  # Proxy vector, mean, and variance for selected sample
  Xmean_s <- mean(X_s)
  Xvar_s <- sum((X_s - Xmean_s)^2) / length(X_s)

  # Polyserial correlation and threshold
  # TWO-STEP METHOD
  # Cutpoint fixed
  w <- qnorm(1 - mean(y_s))
  # Maximize likelihood wrt p, holding w constant
  # Likelihood containing (p)
  f <- function(pars) {
    p <- pars[1]
    a <- -(w + Xmean_s * p / sqrt(Xvar_s)) / sqrt(1 - p^2)
    b <- (p / sqrt(Xvar_s)) / sqrt(1 - p^2)
    logPhi <- pnorm(a + b * X_s, log.p = T)
    log1minusPhi <- pnorm(a + b * X_s, log.p = T, lower.tail = F)
    -sum(y_s * logPhi + (1 - y_s) * log1minusPhi)
  }
  result <- optimize(f, interval = c(-0.99, 0.99))
  rho_s <- result$minimum

  # MLEs for distribution of U
  g <- (phi + (1 - phi) * rho_s) / (phi * rho_s + (1 - phi))

  umean_s <- -w
  uvar_s <- 1
  xucov_s <- rho_s * sqrt(Xvar_s)
  umean_ns <- umean_s + g * (Xmean_ns - Xmean_s) / sqrt(Xvar_s)
  uvar_ns <- 1 + g^2 * (Xvar_ns - Xvar_s) / Xvar_s
  # If uvar_ns < 0 replace with boundary value .Machine$double.eps
  # This will cause pnorm(umean_ns/sqrt(uvar_ns)) = +1 or -1 depending on sign of umean_ns
  uvar_ns <- ifelse(uvar_ns < 0, .Machine$double.eps, uvar_ns)

  # Calculate sfrac if not provided
  if (is.null(sfrac)) {
    sfrac <- length(y_s) / (length(y_s) + sumry_Z_ns$n_Z)
  }

  # MLEs for distribution of Y and the MSB
  ## E[Y|M=0]
  ymean_s <- mean(y_s) # same as pnorm(umean_s) b/c 2-step estimator
  ## E[Y|M=1]
  ymean_ns <- pnorm(umean_ns / sqrt(uvar_ns))
  ## E[Y]
  ymean <- sfrac * ymean_s + (1 - sfrac) * ymean_ns
  ## MSB
  msb <- ymean_s - ymean

  sumry_table <- data.frame(
    phi = phi,
    pi = sfrac,
    muX_s = Xmean_s,
    muU_s = umean_s,
    muY_s = ymean_s,
    muX_ns = Xmean_ns,
    muU_ns = umean_ns,
    muY_ns = ymean_ns,
    sigmaXX_s = Xvar_s,
    sigmaXU_s = xucov_s,
    sigmaUU_s = uvar_s,
    sigmaXX_ns = Xvar_ns,
    sigmaUU_ns = uvar_ns,
    rho_s = rho_s,
    muY = ymean,
    msb = msb
  )

  return(sumry_table)
}
