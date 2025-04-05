################## Error Check #############################
#' Check if error in the input data for calculating outcome mean estimates under proxy pattern-mixture model
#'
#' `ErrorCheck` examinate if the input data for calculating outcome mean estimates under proxy pattern-mixture model is correct.
#' @param data_YZ_s list for selected sample (Microdata or Summary statistics). See details below for more details.
#' @param data_Z_ns list for non-selected sample (Microdata or Summary statistics). See details below for more details.
#' @param prop_check if `TRUE`, check if Y is binary for binary outcome. Default is `FALSE`.
#' @param mi_check if `TRUE`, check if data_YZ_s and data_Z_ns are both microdata for multiple imputation. Default is `FALSE`.
#'
#' @returns A list of data_YZ_s, data_Z_ns, sumry_YZ_s and sumry_Z_ns with
#' \itemize{
#' \item `data_YZ_s=data_YZ_s` : the same as the input data_YZ_s
#' \item `data_Z_ns=data_Z_ns` : the same as the input data_Z_ns
#' \item `sumry_YZ_s=list(mean_YZ, var_YZ, n_YZ)` : summary statistics for selected sample
#' \item `sumry_Z_ns=list(mean_Z, var_Z, n_Z)` : summary statistics for non-selected sample
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
#'  \item{`data_YZA=list(Y,Z)`(Microdata)}{
#'    \itemize{
#'      \item `Y` : `n*1` matrix of outcome variable `Y`
#'      \item `Z` : `n*nZparams` matrix of auxiliary variables `Z`
#'    }
#'  }
#'  \item{`data_YZA=list(mean_YZ, var_YZ, n_YZ)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_YZ` : `(1+nZparams)*1` vector of mean (in order) of outcome variable `Y` and auxiliary variables `Z`
#'      \item `var_YZ` : `(1+nZparams+nAparams)*(1+nZparams+nAparams)` matrix of variance (in order) of outcome variable `Y` and auxiliary variables `Z`
#'      \item `n_YZ` : number of selected sample
#'    }
#'  }
#'
#'  \item{`data_Z_ns=list(Z)`(Microdata)}{
#'    \itemize{
#'      \item `Z` : `(N-n)*nZparams` matrix of auxiliary variables `Z`
#'    }
#'  }
#'  \item{`data_Z_ns=list(mean_Z, var_Z, n_Z)`(Summary Statistics)}{
#'    \itemize{
#'      \item `mean_Z` : `(nZparams)*1` vector of mean (in order) of auxiliary variables `Z`
#'      \item `var_Z` : `(nZparams)*(nZparams)` matrix of variance (in order) of auxiliary variables `Z`
#'      \item `n_Z` : number of non-selected sample
#'    }
#'  }
#'
#' }
#'
ErrorCheck <- function(data_YZ_s, # selected/non-missing sample
                       data_Z_ns, # non-selected/missing sample
                       prop_check = FALSE, # Y default not binary
                       mi_check = FALSE) { # default not multiple imputation

  is_YZ_micro <- FALSE
  is_Z_micro <- FALSE
  ## test if list of Y or Z missing in data_YZ_s
  if (all(c("Y", "Z") %in% names(data_YZ_s))) { # microdata
    is_YZ_micro <- TRUE
    # check if data_YZ_s has missing (need Y no missing & Z no missing)
    if (any(is.na(data_YZ_s$Y)) | any(is.na(data_YZ_s$Z))) { # if Y missing or Z missing
      stopfn("Requires complete Microdata data_YZ_s (list of Y and Z)")
    }
    data_YZ_s$Y <- as.matrix(data_YZ_s$Y)
    data_YZ_s$Z <- as.matrix(data_YZ_s$Z)

    if (ncol(data_YZ_s$Y) != 1) { # check if the column of the matrix is not 1
      stopfn("Microdata data_YZ_s need to have only one column of Y (both vector and matrix work)")
    }
    if (ncol(data_YZ_s$Z) == 1) {
      colnames(data_YZ_s$Z) <- c("Z")
    }

    # check if Y is not a number -> change to 0/1 indicator
    if (!is.numeric(data_YZ_s$Y)) {
      if (length(unique(data_YZ_s$Y)) != 2) {
        stopfn("Requires binary Y")
      } else {
        # change factor Y -> binary Y
        Binary_Y <- fastDummies::dummy_cols(data.frame(Y = c(data_YZ_s$Y)),
          "Y",
          remove_first_dummy = T,
          remove_selected_columns = T
        )
        # change Y to matrix Y
        data_YZ_s$Y <- as.matrix(Binary_Y)
        cat(
          "Note: Modeling probability Y =",
          str_sub(names(Binary_Y), 3), "\n"
        )
      }
    }
    # if prop -> binary
    if (prop_check) {
      if (!identical(c(0, 1), as.double(sort(unique(data_YZ_s$Y))))) {
        stopfn("Requires binary Y")
      }
    }
    # check if Y and Z have same dimension
    if (nrow(data_YZ_s$Y) != nrow(data_YZ_s$Z)) {
      stopfn("Microdata data_YZ_s need to have same observation numbers for Y and Z")
    }
    # turn Microdata data_YZ_s to Summary Statistics sumry_YZ_s
    sumry_YZ_s <- MicroToSummary(data_YZ_s$Z, data_YZ_s$Y)
  } else if (all(c("mean_YZ", "var_YZ", "n_YZ") %in% names(data_YZ_s))) { # summary
    if (mi_check) { # check for multiple imputation
      stopfn("Requires Microdata data_YZ_s (list of Y and Z)")
    }
    if (prop_check) { # check for proportion
      stopfn("Requires Microdata data_YZ_s (list of Y and Z)")
    }
    # check if data_YZ_s has missing
    if (any(is.na(data_YZ_s$mean_YZ)) | any(is.na(data_YZ_s$var_YZ)) | is.na(data_YZ_s$n_YZ)) { # if mean/var/n missing
      stopfn("Requires complete Summary Statistics data_YZ_s (list of mean_YZ, var_YZ and n_YZ)")
    }
    # check if length(mean) = ncol(covariance matrix) = nrows(covariance matrix)
    if ((length(data_YZ_s$mean_YZ) != nrow(data_YZ_s$var_YZ)) | (length(data_YZ_s$mean_YZ) != ncol(data_YZ_s$var_YZ))) {
      stopfn("Summary Statistics data_YZ_s need to have same number of covariates for Y and Z")
    }
    # return Summary Statistics
    sumry_YZ_s <- data_YZ_s
  } else {
    stopfn("Requires data_YZ_s to be Microdata (list of Y and Z) or Summary Statistics(list of mean_YZ, var_YZ and n_YZ)")
  }


  ## test if list of Z missing in data_Z_ns
  if ("Z" %in% names(data_Z_ns)) { # microdata data_Z_ns
    is_Z_micro <- TRUE
    # check if data_Z_ns has missing (need Z no missing)
    if (any(is.na(data_Z_ns$Z))) { # if Z missing
      stopfn("Requires complete Microdata data_Z_ns (list of Z)")
    }
    data_Z_ns$Z <- as.matrix(data_Z_ns$Z)

    if (ncol(data_Z_ns$Z) == 1) {
      colnames(data_Z_ns$Z) <- c("Z")
    }
    # turn Microdata data_Z_ns to Summary Statistics sumry_Z_ns
    sumry_Z_ns <- MicroToSummary(data_Z_ns$Z)
  } else if (all(c("mean_Z", "var_Z", "n_Z") %in% names(data_Z_ns))) { # summary
    if (mi_check) { # check for multiple imputation
      stopfn("Requires Microdata data_Z_ns (list of Z)")
    }
    # check if data_Z_ns has missing
    if (any(is.na(data_Z_ns$mean_Z)) | any(is.na(data_Z_ns$var_Z)) | is.na(data_Z_ns$n_Z)) { # if mean/var/n missing
      stopfn("Requires complete Summary Statistics data_Z_ns (list of mean_Z, var_Z and n_Z)")
    }
    data_Z_ns$var_Z <- as.matrix(data_Z_ns$var_Z)
    if (ncol(data_Z_ns$var_Z) == 1) { # check if the column of the matrix is 1
      colnames(data_Z_ns$var_Z) <- c("Z")
      rownames(data_Z_ns$var_Z) <- c("Z")
    }

    # check if length(mean) = ncol(covariance matrix) = nrows(covariance matrix)
    if ((length(data_Z_ns$mean_Z) != nrow(data_Z_ns$var_Z)) | (length(data_Z_ns$mean_Z) != ncol(data_Z_ns$var_Z))) {
      stopfn("Summary Statistics data_Z_ns need to have same number of covariates for Z")
    }
    # return Summary Statistics
    sumry_Z_ns <- data_Z_ns
  } else {
    stopfn("Requires data_Z_ns to be Microdata (list of Z) or Summary Statistics(list of mean_Z, var_Z and n_Z)")
  }

  ## check if Z match in data_YZ_s and data_Z_ns
  # check if covariates Z matched in data_YZ_s and data_Z_ns
  if (!identical((colnames(sumry_YZ_s$var_YZ)[-1]), colnames(sumry_Z_ns$var_Z))) {
    stopfn("data_YZ_s and data_Z_ns need to have exactly same covariates for Z")
  }
  if (!is_YZ_micro) {
    data_YZ_s <- NULL
  }
  if (!is_Z_micro) {
    data_Z_ns <- NULL
  }
  # return data_YZ_s, data_Z_ns, sumry_YZ_s, sumry_Z_ns
  return(list(
    data_YZ_s = data_YZ_s,
    data_Z_ns = data_Z_ns,
    sumry_YZ_s = sumry_YZ_s,
    sumry_Z_ns = sumry_Z_ns
  ))
}

############### Turn Microdata to Summary Statistics ###############
MicroToSummary <- function(Z, # microdata of Z
                           Y = NULL) {
  if (is.null(Y)) {
    # Y is NULL
    dt <- Z
  } else {
    # combine Z and Y if Y has value
    dt <- cbind(Y, Z)
  }
  # Mean for (Y,Z)
  dt_mean <- colMeans(dt)
  # Covariance matrix for (Y,Z)
  dt_var <- as.matrix(stats::var(dt))

  if (ncol(dt_var) == 1) {
    colnames(dt_var) <- c("Z")
    rownames(dt_var) <- c("Z")
  }

  # Sample size
  n <- nrow(dt)
  # Summary statistics for non-selected/missing
  if (is.null(Y)) {
    sumry <- list(
      mean_Z = dt_mean,
      var_Z = dt_var,
      n_Z = n
    )
  } else { # Summary statistics for selected/non-missing
    sumry <- list(
      mean_YZ = dt_mean,
      var_YZ = dt_var,
      n_YZ = n
    )
  }
  return(sumry)
}
