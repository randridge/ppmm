################# Error Check for Coefficient in Regression##########################
#' Error Check for Coefficient in Regression
#'
#' `Coeff_ErrorCheck` examinate if the input data for calculating coefficient in regression under proxy pattern-mixture model is correct.
#'
#' @param data_YZA_s list for selected sample (Microdata or Summary statistics). See details below for more details.
#' @param data_ZA_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param prop_check if `TRUE`, check if Y is binary for binary outcome. Default is `FALSE`.
#'
#' @returns A list of data_YZA_s, sumry_YZA_s and sumry_ZA_ns with
#' \itemize{
#' \item `data_YZA_s=data_YZA_s` : the same as the input data_YZA_s
#' \item `sumry_YZA_s=list(mean_YZA, var_YZA, n_YZA)` : summary statistics for selected sample
#' \item `sumry_ZA_ns=list(mean_ZA, var_ZA, n_ZA)` : summary statistics for non-selected sample
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
#'

Coeff_ErrorCheck <- function(data_YZA_s,
                             data_ZA_ns,
                             prop_check = FALSE) {
  is_YZA_micro <- FALSE
  is_ZA_micro <- FALSE
  ## test if list of Y or Z or A missing in data_YZA_s
  if (all(c("Y", "Z", "A") %in% names(data_YZA_s))) { # microdata
    is_YZA_micro <- TRUE
    # check if data_YZA_s has missing (need Y no missing & Z no missing & A no missing)
    if (any(is.na(data_YZA_s$Y)) | any(is.na(data_YZA_s$Z)) | any(is.na(data_YZA_s$A))) { # if Y missing or Z missing or A missing
      stopfn("Requires complete Microdata data_YZA_s (list of Y, Z and A)")
    }
    data_YZA_s$Y <- as.matrix(data_YZA_s$Y)
    data_YZA_s$Z <- as.matrix(data_YZA_s$Z)
    data_YZA_s$A <- as.matrix(data_YZA_s$A)

    if (ncol(data_YZA_s$Y) != 1) { # check if the column of the matrix is not 1
      stopfn("Microdata data_YZA_s need to have only one column of Y (both vector and matrix work)")
    }
    if (ncol(data_YZA_s$Z) == 1) {
      colnames(data_YZA_s$Z) <- c("Z")
    }
    if (ncol(data_YZA_s$A) == 1) {
      colnames(data_YZA_s$A) <- c("A")
    }
    # check if Y is not a number -> change to 0/1 indicator
    if (!is.numeric(data_YZA_s$Y)) {
      if (length(unique(data_YZA_s$Y)) != 2) {
        stopfn("Requires binary Y")
      } else {
        # change factor Y -> binary Y
        Binary_Y <- fastDummies::dummy_cols(data.frame(Y = c(data_YZA_s$Y)),
          "Y",
          remove_first_dummy = T,
          remove_selected_columns = T
        )
        # change Y to matrix Y
        data_YZA_s$Y <- as.matrix(Binary_Y)
        cat(
          "Input Y in data_YZA_s is string -> make Y to be indicator of",
          str_sub(names(Binary_Y), 3), "\n"
        )
      }
    }
    # if prop -> binary
    if (prop_check) {
      if (!identical(c(0, 1), as.double(sort(unique(data_YZA_s$Y))))) {
        stopfn("Requires binary Y")
      }
    }

    # check if Y, Z and A have same dimension
    if (nrow(data_YZA_s$Y) != nrow(data_YZA_s$Z) |
      nrow(data_YZA_s$Y) != nrow(data_YZA_s$A)) {
      stopfn("Microdata data_YZA_s need to have same observation numbers for Y, Z and A")
    }
    # turn Microdata data_YZA_s to Summary Statistics sumry_YZA_s
    sumry_YZA_s <- Coeff_MicroToSummary(data_YZA_s$Z, data_YZA_s$A, data_YZA_s$Y)
  } else if (all(c("mean_YZA", "var_YZA", "n_YZA") %in% names(data_YZA_s))) { # summary
    if (prop_check) { # check for proportion
      stopfn("Requires Microdata data_YZA_s (list of Y, Z and A)")
    }
    # check if data_YZA_s has missing
    if (any(is.na(data_YZA_s$mean_YZ)) | any(is.na(data_YZA_s$var_YZ)) |
      is.na(data_YZA_s$n_YZ)) { # if mean/var/n missing
      stopfn("Requires complete Summary Statistics data_YZA_s (list of mean_YZA, var_YZA and n_YZA)")
    }
    # check if length(mean) = ncol(covariance matrix) = nrows(covariance matrix)
    if ((length(data_YZA_s$mean_YZA) != nrow(data_YZA_s$var_YZA)) | (length(data_YZA_s$mean_YZA) != ncol(data_YZA_s$var_YZA))) {
      stopfn("Summary Statistics data_YZA_s need to have same number of covariates for Y, Z and A")
    }
    # return Summary Statistics
    sumry_YZA_s <- data_YZA_s
  } else {
    stopfn("Requires data_YZA_s to be Microdata (list of Y, Z and A) or Summary Statistics(list of mean_YZA, var_YZA and n_YZA)")
  }



  ## test if list of Z or A missing in data_ZA_ns
  if (all(c("Z", "A") %in% names(data_ZA_ns))) { # microdata data_ZA_ns
    is_ZA_micro <- TRUE
    # check if data_ZA_ns has missing (need Z no missing, A no missing)
    if (any(is.na(data_ZA_ns$Z)) | any(is.na(data_ZA_ns$A))) { # if Z missing or A missing
      stopfn("Requires complete Microdata data_ZA_ns (list of Z and A)")
    }
    data_ZA_ns$Z <- as.matrix(data_ZA_ns$Z)
    data_ZA_ns$A <- as.matrix(data_ZA_ns$A)

    if (ncol(data_ZA_ns$Z) == 1) {
      colnames(data_ZA_ns$Z) <- c("Z")
    }
    if (ncol(data_ZA_ns$A) == 1) {
      colnames(data_ZA_ns$A) <- c("A")
    }

    # turn Microdata data_ZA_ns to Summary Statistics sumry_ZA_ns
    sumry_ZA_ns <- Coeff_MicroToSummary(data_ZA_ns$Z, data_ZA_ns$A)
  } else if (all(c("mean_ZA", "var_ZA", "n_ZA") %in% names(data_ZA_ns))) { # summary
    # check if data_ZA_ns has missing
    if (any(is.na(data_ZA_ns$mean_ZA)) | any(is.na(data_ZA_ns$var_ZA)) | is.na(data_ZA_ns$n_ZA)) { # if mean/var/n missing
      stopfn("Requires complete Summary Statistics data_ZA_ns (list of mean_ZA, var_ZA and n_ZA)")
    }
    data_ZA_ns$var_ZA <- as.matrix(data_ZA_ns$var_ZA)
    if (ncol(data_ZA_ns$var_ZA) == 1) {
      colnames(data_ZA_ns$var_ZA) <- colnames(sumry_YZA_s$var_YZA)
      rownames(data_ZA_ns$var_ZA) <- rownames(sumry_YZA_s$var_YZA)
    }
    # check if length(mean) = ncol(covariance matrix) = nrows(covariance matrix)
    if ((length(data_ZA_ns$mean_ZA) != nrow(data_ZA_ns$var_ZA)) |
      (length(data_ZA_ns$mean_ZA) != ncol(data_ZA_ns$var_ZA))) {
      stopfn("Summary Statistics data_ZA_ns need to have same number of covariates for Z and A")
    }
    # return Summary Statistics
    sumry_ZA_ns <- data_ZA_ns
  } else {
    stopfn("Requires data_ZA_ns to be Microdata (list of Z and A) or Summary Statistics(list of mean_ZA, var_ZA and n_ZA)")
  }

  # check if covariates ZA (in order) matched in data_YZA_s and data_ZA_ns
  if (!identical((colnames(sumry_YZA_s$var_YZA)[-1]), colnames(sumry_ZA_ns$var_ZA))) {
    stopfn("data_YZA_s and data_ZA_ns need to have exactly same covariates for Z and A in order")
  }

  if (!is_YZA_micro) {
    data_YZA_s <- NULL
  }
  if (!is_ZA_micro) {
    data_ZA_ns <- NULL
  }
  # return data_YZA_s, sumry_YZA_s, sumry_ZA_ns
  return(list(
    data_YZA_s = data_YZA_s,
    data_ZA_ns = data_ZA_ns,
    sumry_YZA_s = sumry_YZA_s,
    sumry_ZA_ns = sumry_ZA_ns
  ))
}

Coeff_MicroToSummary <- function(Z, #
                                 A,
                                 Y = NULL) {
  if (is.null(Y)) {
    # Y is NULL, combine Z, A
    dt <- cbind(Z, A)
  } else {
    # combine Z, A and Y if Y has value
    dt <- cbind(Y, Z, A)
  }
  # Mean for (Y,Z,A)
  dt_mean <- colMeans(dt)
  # Covariance matrix for (Y,Z,A)
  dt_var <- stats::var(dt)

  if (ncol(dt_var) == 1) {
    colnames(dt_var) <- c("Z")
    rownames(dt_var) <- c("Z")
  }
  # Sample size
  n <- nrow(dt)

  # Summary statistics for non-selected/missing
  if (is.null(Y)) {
    sumry <- list(
      mean_ZA = dt_mean,
      var_ZA = dt_var,
      n_ZA = n
    )
  } else { # Summary statistics for selected/non-missing
    sumry <- list(
      mean_YZA = dt_mean,
      var_YZA = dt_var,
      n_YZA = n
    )
  }
  return(sumry)
}
