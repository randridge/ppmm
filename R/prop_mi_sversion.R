#' Function for Proxy Pattern-Mixture Analysis (Binary Outcomes)
#'
#' `prop_mi` calculate multiply imputed outcome for binary outcomes.
#'
#' @param data_YZ_s list of Y, Z for selected sample (Microdata only)
#' @param data_Z_ns list of Z for non-selected sample (Microdata only)
#' @param phi value of phi, fixed to a value between 0 and 1 (defaults to 0)
#' @param phi_character draw phi from specific distribution (i.e."`runif(1)`", defaults to `NULL`)
#' * If `phi_character=NULL`, fix at value set to `phi`;
#' * If `phi_character=` specific distribution (i.e. `runif(1)` or `rbeta(1,1,1)`), then the value of `phi` input to the function is ignored
#' @param D number of multiply imputed data sets (defaults to 10)
#' @param burnin number of burn-in replicates (defaults to 500)
#' @param thin number of replicates between imputations (defaults to 100)
#'
#' @returns A `N*D` matrix of multiply imputed binary Y vectors, with each column an imputed outcome data set.
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
#' @seealso [cont_mi()]
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
#' # impute the missing outcome data
#' imp <- prop_mi(data_YZ_s, data_Z_ns, phi = 0, D = 10, burnin = 500, thin = 100)

prop_mi <- function(data_YZ_s,
                    data_Z_ns,
                    phi = 0,
                    phi_character = NULL,
                    D = 10,
                    burnin = 500,
                    thin = 100) {
  ### Some functions to do truncated normals in R
  rnorm.lt <- function(n, lv = rep(0, n), mv = rep(0.5, n), sv = rep(1, n)) {
    lstd <- (lv - mv) / sv
    sv * qnorm(runif(n, pnorm(lstd), rep(1, n))) + mv
  }
  rnorm.rt <- function(n, rv = rep(0, n), mv = rep(0.5, n), sv = rep(1, n)) {
    rstd <- (rv - mv) / sv
    sv * qnorm(runif(n, rep(0, n), pnorm(rstd))) + mv
  }

  sumry_list <- ErrorCheck(data_YZ_s, data_Z_ns, prop_check = TRUE, mi_check = TRUE)
  data_YZ_s <- sumry_list$data_YZ_s
  data_Z_ns <- sumry_list$data_Z_ns

  # Total n
  n_s <- length(data_YZ_s$Y) # sample size for selected sample
  n_ns <- nrow(data_Z_ns$Z) # sample size for non-selected sample
  n <- n_s + n_ns

  # Number of selected cases where Y=1 and where Y=0
  n_s_y1 <- sum(data_YZ_s$Y)
  n_s_y0 <- n_s - n_s_y1

  # Probit regression of Y|Z for selected cases to find starting point for Gibbs sampler
  fit <- glm(data_YZ_s$Y ~ data_YZ_s$Z, family = binomial(link = "probit"))
  betaHat <- as.vector(fit$coef)
  z_s <- model.matrix(fit) # attaches column of 1s (selected sample)
  z_ns <- cbind(rep(1, nrow(data_Z_ns$Z)), data_Z_ns$Z) # attaches column of 1s (non-selected sample)
  zmean_ns <- colMeans(z_ns) # column mean of non-selected sample
  zvar_ns <- var(data_Z_ns$Z) * (n_ns - 1) / n_ns # variance of non-selected sample
  z_sTz_sInv <- solve(t(z_s) %*% z_s) # so only have to calculate once

  # Starting point from probit fit
  B <- as.matrix(betaHat)

  # Starting proxy value
  X_s <- z_s %*% B # selected sample (microdata)
  X_ns <- z_ns %*% B # non-selected sample (microdata)
  Xmean_ns <- zmean_ns %*% B # non-selected sample (sum stats only)
  Xvar_ns <- t(B[-1]) %*% zvar_ns %*% B[-1] # non-selected (sum stats only)

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
    # Intialize latent U vector (Selected sampleonly)
    u <- rep(NA, n_s)

    ## (0) Draw phi from prior specified by phi_character
    if (!is.null(phi_character)) {
      # Draw phi using provided string
      phi <- eval(parse(text = phi_character))
      phi <- pmax(.Machine$double.eps, phi)
    }

    ## (1) Draw latent U|Y,B for selected sample (aka U|Y,X)
    u[data_YZ_s$Y == 1] <- rnorm.lt(n_s_y1, mv = X_s[data_YZ_s$Y == 1]) # Draws for Y=1 --> truncated at left by 0
    u[data_YZ_s$Y == 0] <- rnorm.rt(n_s_y0, mv = X_s[data_YZ_s$Y == 0]) # Draws for Y=0 --> truncated at right by 0

    ## (2) Draw B|Z,U
    B <- t(rmvnorm(1, z_sTz_sInv %*% t(z_s) %*% u, z_sTz_sInv))

    ## (3) Create proxy given current B
    X_s <- z_s %*% B # selected sample (microdata)
    X_ns <- z_ns %*% B # non-selected sample (microdata)
    Xmean_ns <- zmean_ns %*% B # non-selected sample (sum stats only)
    Xvar_ns <- t(B[-1]) %*% zvar_ns %*% B[-1] # Missing sample (sum stats only)

    ##### IF REPLICATE IS TO BE USED FOR IMPUTATION ####
    # Proceed with drawing PPM parameters conditional on drawn U,B
    if ((i > burnin) & (((i - 1 - burnin) %% thin) == 0)) {
      d <- d + 1
      # print(c("INNER LOOP :",d))
      ############ Step 1 ############
      ## (4) Scale proxy X to have same variance as latent U among respondents
      # Draw the population variances of X, U from posterior
      varXdraw <- sum((X_s - mean(X_s))^2) / rchisq(1, n_s - 1)
      varUdraw <- sum((u - mean(u))^2) / rchisq(1, n_s - 1)
      # Use draws to scale the proxy
      X_s <- X_s * sqrt(varUdraw / varXdraw) # Non-Missing sample (microdata)
      X_ns <- X_ns * sqrt(varUdraw / varXdraw) # Missing sample (microdata)
      Xmean_ns <- Xmean_ns * sqrt(varUdraw / varXdraw) # Non-selected sample (sum stats only)
      Xvar_ns <- Xvar_ns * varUdraw / varXdraw # Non-selected sample (sum stats only)


      ## (5) Draw from PPM dependent on value of phi, using (X,U)
      if (phi == 0) {
        y1_0 <- X_s # Selected sample (microdata)
        y2_0 <- u # Selected sample (microdata)
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
          y2_0 <- u # Selected sample (microdata)
          y1_1_mean <- Xmean_ns # Non-selected sample (sum stats only)
          y1_1_var <- Xvar_ns # Non-selected sample (sum stats only)
        } else {
          y1_0 <- X_s # Selected sample (microdata)
          y2_0 <- (1 - phi) * X_s + phi * u # Selected sample (microdata)
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
      if (phi != 0 & phi != 1) {
        # Transform draws of [X,W] to get draws from [X,U]
        # W = (1-phi)*X + phi*U --> U = (W - (1-phi)*X)/phi
        # Start with draws of [X,W] and then overwrite parms relating to U
        drawsXW <- drawsPPM
        drawsPPM$mu2_0 <- (drawsXW$mu2_0 - (1 - phi) * drawsXW$mu1_0) / phi
        drawsPPM$mu2_1 <- (drawsXW$mu2_1 - (1 - phi) * drawsXW$mu1_1) / phi
        drawsPPM$sigma22_0 <- (drawsXW$sigma22_0 + (1 - phi)^2 * drawsXW$sigma11_0 - 2 * (1 - phi) * drawsXW$sigma12_0) / phi^2
        drawsPPM$sigma22_1 <- (drawsXW$sigma22_1 + (1 - phi)^2 * drawsXW$sigma11_1 - 2 * (1 - phi) * drawsXW$sigma12_1) / phi^2
        drawsPPM$sigma12_0 <- (drawsXW$sigma12_0 - (1 - phi) * drawsXW$sigma11_0) / phi
        drawsPPM$sigma12_1 <- (drawsXW$sigma12_1 - (1 - phi) * drawsXW$sigma11_1) / phi
      }

      #########################
      ## (8) Transform parameters to get draws of the mean of Y (unmodified Bayes method)
      #### Estimate of mean of Y for selected sample
      # drawsPPM$muY_0 <- pnorm(drawsPPM$mu2_0/sqrt(drawsPPM$sigma22_0))
      #### Estimate of mean of Y for non-selected sample
      # drawsPPM$muY_1 <- pnorm(drawsPPM$mu2_1/sqrt(drawsPPM$sigma22_1))
      #### Estimate of mean of Y (overall)
      # drawsPPM$muY <- drawsPPM$pi*drawsPPM$muY_0 + (1-drawsPPM$pi)*drawsPPM$muY_1
      #########################

      ############ Step 2 ############

      ## (6) Draw U for nonrespondents given PPM parameters and current value of X
      ### Calculate conditional mean and variance for [U|X,M=1]
      cmean <- drawsPPM$mu2_1 + (drop(drawsPPM$sigma12_1 / drawsPPM$sigma11_1)) * (X_ns - drawsPPM$mu1_1)
      cvar <- drawsPPM$sigma22_1 - drawsPPM$sigma12_1^2 / drawsPPM$sigma11_1
      ### Draw U for nonrespondents/nonselected
      u_ns <- rnorm(n_ns, mean = cmean, sd = sqrt(cvar))
      ### imputed Y=1 if U>0 ##
      imp[(n_s + 1):n, d] <- (u_ns > 0)
    }
  }
  # End looping

  return(imp)
}
