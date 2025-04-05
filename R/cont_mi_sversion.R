#' Function for Proxy Pattern-Mixture Analysis (Continuous Outcomes)
#'
#' `conti_mi` calculate multiply imputed outcome for continuous outcomes.
#'
#' @param data_YZ_s list of Y, Z for selected sample (Microdata only)
#' @param data_Z_ns list of Z for non-selected sample (Microdata only)
#' @param phi value of phi, fixed to a value between 0 and 1 (defaults to 0)
#' @param phi_character draw phi from specific distribution (i.e."`runif(1)`", defaults to `NULL`)
#' * If `phi_character=NULL`, fix at value set to `phi`;
#' * If `phi_character=` specific distribution (i.e. `runif(1)` or `rbeta(1,1,1)`), then the value of `phi` input to the function is ignored
#' @param D number of multiply imputed data sets (defaults to 10)
#' @param burnin number of burn-in replicates (defaults to 0)
#' @param thin number of replicates between imputations (defaults to 1)
#'
#' @returns A `N*D` matrix of multiply imputed Y vectors, with each column an imputed outcome data set.
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
#' }
#'
#' @seealso [prop_mi()]
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
#' # impute the missing outcome data
#' imp <- cont_mi(data_YZ_s, data_Z_ns)

cont_mi <- function(data_YZ_s,
                    data_Z_ns,
                    phi = 0,
                    phi_character = NULL,
                    D = 10,
                    burnin = 0,
                    thin = 1)
{
  sumry_list <- ErrorCheck(data_YZ_s, data_Z_ns, prop_check = FALSE, mi_check = TRUE)
  data_YZ_s <- sumry_list$data_YZ_s
  data_Z_ns <- sumry_list$data_Z_ns

  # Total n
  n_s <- length(data_YZ_s$Y) # sample size for Selected/Respondents/Non-Missing sample
  n_ns <- nrow(data_Z_ns$Z) # sample size for Non-selected/Missing sample
  n <- n_s + n_ns

  # Regression of Y|Z
  fit <- lm(data_YZ_s$Y ~ data_YZ_s$Z)
  betaHat <- as.vector(fit$coef)
  Vb <- summary(fit)$cov.unscaled
  s2 <- summary(fit)$sigma^2
  dfResid <- fit$df.residual

  ##############################
  z_s <- model.matrix(fit) # attaches column of 1s (Non-Missing sample)
  z_ns <- cbind(rep(1, nrow(data_Z_ns$Z)), data_Z_ns$Z) # attaches column of 1s (Missing sample)
  zmean_ns <- colMeans(z_ns) # column mean of missing sample
  zvar_ns <- var(data_Z_ns$Z) * (n_ns - 1) / n_ns # variance of missing sample

  ##############################

  # Initialize vector(matrix) to hold multiply-imputed Y
  imp <- matrix(nrow = n, ncol = D)
  imp[1:n_s, ] <- data_YZ_s$Y # Fill in observed values of Y

  # Total number of replicates needed is burnin+1+(D-1)*thin
  # Draw U, B on all replicates
  # Only draw PPM parameters on the replicates to perform MI
  # Keep track of which iterate you're on
  d <- 0
  # Burn-in iterations for drawing (U,B)
  for (i in 1:(burnin + 1 + (D - 1) * thin))
  {
    #    if ((j %% 100)==0) print(paste("Iteration",j))

    ## (0) Draw phi from prior specified by phi_character
    if (!is.null(phi_character)) {
      # Draw phi using provided string
      phi <- eval(parse(text = phi_character))
      phi <- pmax(.Machine$double.eps, phi)
    }

    # Draw parameter values from posterior dependent on value of lambda

    ## (1) Draw SIGMA^2 from inverse-ChiSq
    SIGMA2 <- dfResid * s2 / rchisq(1, dfResid)

    ## (2) Draw BETA | SIGMA, Y
    B <- t(rmvnorm(1, betaHat, Vb * SIGMA2))

    ## (3) Create proxy given current B
    X_s <- z_s %*% B # Selected/Non-Missing sample (microdata)
    X_ns <- z_ns %*% B # Non-selected/Missing sample (microdata)
    Xmean_ns <- zmean_ns %*% B # Non-selected/Missing sample (sum stats only)
    Xvar_ns <- t(B[-1]) %*% zvar_ns %*% B[-1] # Non-selected/Missing sample (sum stats only)

    ##### IF REPLICATE IS TO BE USED FOR IMPUTATION ####
    # Proceed with drawing PPM parameters conditional on drawn U,B
    if ((i > burnin) & (((i - 1 - burnin) %% thin) == 0)) {
      d <- d + 1
      # print(c("INNER LOOP :",d))
      ############ Step 1 ############
      ## (4) Scale proxy X to have same variance as latent U among respondents
      # Draw the population variances of X, Y from posterior
      varYdraw <- sum((data_YZ_s$Y - mean(data_YZ_s$Y))^2) / rchisq(1, n_s - 1)
      varXdraw <- sum((X_s - mean(X_s))^2) / rchisq(1, n_s - 1)
      # Use draws to scale the proxy
      X_s <- X_s * sqrt(varYdraw / varXdraw) # Non-Missing sample (microdata)
      X_ns <- X_ns * sqrt(varYdraw / varXdraw) # Missing sample (microdata)
      Xmean_ns <- Xmean_ns * sqrt(varYdraw / varXdraw) # Non-selected sample (sum stats only)
      Xvar_ns <- Xvar_ns * varYdraw / varXdraw # Non-selected sample (sum stats only)


      ## (5) Draw from PPM dependent on value of phi
      if (phi == 0) {
        y1_0 <- X_s # Selected sample (microdata)
        y2_0 <- data_YZ_s$Y # Selected sample (microdata)
        y1_1_mean <- Xmean_ns # Non-selected sample (sum stats only)
        y1_1_var <- Xvar_ns # Non-selected sample (sum stats only)

        # Calculate means, sums of squares
        # Selected sample
        y1Bar_0 <- mean(y1_0)
        y2Bar_0 <- mean(y2_0)
        s11_0 <- sum((y1_0 - y1Bar_0)^2) / n_s
        s22_0 <- sum((y2_0 - y2Bar_0)^2) / n_s
        s12_0 <- sum((y1_0 - y1Bar_0) * (y2_0 - y2Bar_0)) / n_s
        b21.1_0 <- s12_0 / s11_0
        s22.1_0 <- s22_0 - (s12_0)^2 / s11_0
        # Nonrespondent data
        y1Bar_1 <- y1_1_mean
        s11_1 <- y1_1_var
        # Draws
        # (1) PI
        PI <- rbeta(1, n_s + 0.5, n_ns + 0.5)
        # (2) SIGMA11_0
        SIGMA11_0 <- n_s * s11_0 / rchisq(1, n_s - 1)
        # (3) MU1_0 | SIGMA11_0
        MU1_0 <- rnorm(1, y1Bar_0, sqrt(SIGMA11_0 / n_s))
        # (4) SIGMA22.1_0
        SIGMA22.1_0 <- n_s * s22.1_0 / rchisq(1, n_s - 2)
        # (5) BETA21.1_0 | SIGMA22.1_0
        BETA21.1_0 <- rnorm(1, b21.1_0, sqrt(SIGMA22.1_0 / (n_s * s11_0)))
        # (6) BETA20.1_0 | BETA21.1_0, SIGMA22.1_0
        BETA20.1_0 <- rnorm(1, y2Bar_0 - BETA21.1_0 * y1Bar_0, sqrt(SIGMA22.1_0 / n_s))
        # (7) SIGMA11_1
        SIGMA11_1 <- n_ns * s11_1 / rchisq(1, n_ns - 1)
        # (8) MU1_1 | SIGMA11_1
        MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1 / n_ns))
        # Transform draws to get other parameters
        # (a) MU2_0
        MU2_0 <- BETA20.1_0 + BETA21.1_0 * MU1_0
        # (b) MU2_1
        MU2_1 <- BETA20.1_0 + BETA21.1_0 * MU1_1
        # (c) SIGMA12_0
        SIGMA12_0 <- BETA21.1_0 * SIGMA11_0
        # (d) SIGMA22_0
        SIGMA22_0 <- SIGMA22.1_0 + BETA21.1_0^2 * SIGMA11_0
        # (e) SIGMA22_1
        SIGMA22_1 <- SIGMA22.1_0 + BETA21.1_0^2 * SIGMA11_1
        # (f) SIGMA12_1
        SIGMA12_1 <- SIGMA11_1 * BETA21.1_0
        # All Draws
        drawsPPM <- list(
          pi = PI, mu1_0 = MU1_0, mu2_0 = MU2_0, mu1_1 = MU1_1, mu2_1 = MU2_1,
          sigma11_0 = SIGMA11_0, sigma12_0 = SIGMA12_0, sigma22_0 = SIGMA22_0,
          sigma11_1 = SIGMA11_1, sigma12_1 = SIGMA12_1, sigma22_1 = SIGMA22_1
        )
      } else {
        if (phi == 1) {
          y1_0 <- X_s # Selected sample (microdata)
          y2_0 <- data_YZ_s$Y # Selected sample (microdata)
          y1_1_mean <- Xmean_ns # Non-selected sample (sum stats only)
          y1_1_var <- Xvar_ns # Non-selected sample (sum stats only)
        } else {
          y1_0 <- X_s # Selected sample (microdata)
          y2_0 <- (1 - phi) * X_s + phi * (data_YZ_s$Y) # Selected sample (microdata)
          y1_1_mean <- Xmean_ns # Non-selected sample (sum stats only)
          y1_1_var <- Xvar_ns # Non-selected sample (sum stats only)
        }
        # Calculate means, sums of squares
        # Selected sample
        y1Bar_0 <- mean(y1_0)
        y2Bar_0 <- mean(y2_0)
        s11_0 <- sum((y1_0 - y1Bar_0)^2) / n_s
        s22_0 <- sum((y2_0 - y2Bar_0)^2) / n_s
        s12_0 <- sum((y1_0 - y1Bar_0) * (y2_0 - y2Bar_0)) / n_s
        b12.2_0 <- s12_0 / s22_0
        s11.2_0 <- s11_0 - (s12_0)^2 / s22_0
        # Nonrespondent data
        y1Bar_1 <- y1_1_mean
        s11_1 <- y1_1_var
        # Draws
        # (1) PI
        PI <- rbeta(1, n_s + 0.5, n_ns + 0.5)
        # (2) SIGMA22_0
        SIGMA22_0 <- n_s * s22_0 / rchisq(1, n_s - 1)
        # (3) MU2_0 | SIGMA22_0
        MU2_0 <- rnorm(1, y2Bar_0, sqrt(SIGMA22_0 / n_s))
        # (4, 5) SIGMA11.2_0, SIGMA11_1 with constraint
        goodDraw <- FALSE
        ct <- 1
        while (!goodDraw) {
          # Repeat these 2 draws until SIGMA11_1 > SIGMA11.2_0
          # (4) SIGMA11.2_0
          SIGMA11.2_0 <- n_s * s11.2_0 / rchisq(1, n_s - 2)
          # (5) SIGMA11_1
          SIGMA11_1 <- n_ns * s11_1 / rchisq(1, n_ns - 1)
          # Check to see if draws meet the condition
          goodDraw <- (SIGMA11_1 >= SIGMA11.2_0)
          if (ct > 20) {
            goodDraw <- TRUE
            SIGMA11.2_0 <- SIGMA11_1
          }
          ct <- ct + 1
        }
        # (6) BETA12.2_0 | SIGMA11.2_0
        BETA12.2_0 <- rnorm(1, b12.2_0, sqrt(SIGMA11.2_0 / (n_s * s22_0)))
        # (7) BETA10.2_0 | BETA12.2_0, SIGMA11.2_0
        BETA10.2_0 <- rnorm(1, y1Bar_0 - BETA12.2_0 * y2Bar_0, sqrt(SIGMA11.2_0 / n_s))
        # (8) MU1_1 | SIGMA11_1
        MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1 / n_ns))
        # Transform draws to get other parameters
        # (a) MU2_1
        MU2_1 <- (MU1_1 - BETA10.2_0) / BETA12.2_0
        # (b) MU1_0
        MU1_0 <- BETA10.2_0 + BETA12.2_0 * MU2_0
        # (c) SIGMA12_0
        SIGMA12_0 <- BETA12.2_0 * SIGMA22_0
        # (d) SIGMA11_0
        SIGMA11_0 <- SIGMA11.2_0 + BETA12.2_0^2 * SIGMA22_0
        # (e) SIGMA22_1
        SIGMA22_1 <- (SIGMA11_1 - SIGMA11.2_0) / BETA12.2_0^2
        # (f) SIGMA12_1
        SIGMA12_1 <- SIGMA22_1 * BETA12.2_0
        # All Draws
        drawsPPM <- list(
          pi = PI, mu1_0 = MU1_0, mu2_0 = MU2_0, mu1_1 = MU1_1, mu2_1 = MU2_1,
          sigma11_0 = SIGMA11_0, sigma12_0 = SIGMA12_0, sigma22_0 = SIGMA22_0,
          sigma11_1 = SIGMA11_1, sigma12_1 = SIGMA12_1, sigma22_1 = SIGMA22_1
        )
      }

      ############ Step 2 ############
      # Get parameters of conditional dist'n of Y|X,M=1
      if (phi == 0 || phi == 1) {
        cmean <- drawsPPM$mu2_1 + drop(drawsPPM$sigma12_1 / drawsPPM$sigma11_1) * (X_ns - drawsPPM$mu1_1)
        cvar <- drawsPPM$sigma22_1 - (drawsPPM$sigma12_1^2) / drawsPPM$sigma11_1
      } else {
        # Transform draws of [X,W] to get draws from [X,Y]
        # W = (1-phi)*X + phi*Y --> Y = (W - (1-phi)*X)/phi
        mu2_1 <- (drawsPPM$mu2_1 - (1 - phi) * drawsPPM$mu1_1) / phi
        sigma12_1 <- (drawsPPM$sigma12_1 - (1 - phi) * drawsPPM$sigma11_1) / phi
        sigma22_1 <- (drawsPPM$sigma22_1 + (1 - phi)^2 * drawsPPM$sigma11_1 - 2 * (1 - phi) * drawsPPM$sigma12_1) / phi^2
        cmean <- mu2_1 + drop(sigma12_1 / drawsPPM$sigma11_1) * (X_ns - drawsPPM$mu1_1)
        cvar <- sigma22_1 - (sigma12_1^2) / drawsPPM$sigma11_1
      }

      # Impute Y from conditional dist'n of Y|X
      yimp <- rnorm(n_ns, mean = cmean, sd = sqrt(cvar))
      imp[(n_s + 1):n, d] <- yimp
    }
  }
  colnames(imp) <- NULL
  return(imp)
}
