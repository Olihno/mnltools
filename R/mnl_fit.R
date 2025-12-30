#' Log-likelihood for a multinomial logit (MNL) model
#'
#' @param Beta Numeric vector of parameters.
#' @param Data Data frame or data.table in long format.
#' @param X Character vector of column names for explanatory variables.
#' @param choiceVar Name of the binary choice indicator variable.
#' @param GroupID Vector of group identifiers (e.g. individual or choice set).
#'
#' @return Vector of individual log-likelihood contributions.
#' @keywords internal
logLik_mnl <- function(Beta,
                       Data,
                       X,
                       choiceVar,
                       GroupID) {

  # ensure matrix for X
  Xmat <- as.matrix(Data[, X, drop = FALSE])

  V    <- Xmat %*% as.vector(Beta[seq_along(X)])   # utilities
  num  <- exp(V)                                   # numerator
  den  <- tapply(num, GroupID, sum)                # denominator by choice set

  # probabilities for chosen alternatives
  chosen <- Data[[choiceVar]] == 1
  prob   <- num[chosen] / den

  llik <- log(prob)
  return(llik)
}

#' Fit a multinomial logit (MNL) model by maximum likelihood
#'
#' This function estimates a multinomial logit model for one choice occasion per individual
#' using the \pkg{maxLik} package. It can optionally compute the Hessian for
#' standard errors.
#'
#' @param Data Data frame (or data.table) in long format, containing one row per alternative.
#' @param X Character vector of column names for explanatory variables.
#' @param choiceVar Name of the binary choice indicator variable. Defaults to `"choice"`.
#' @param idVar Name of the grouping variable that identifies the choice set
#'   (e.g. individual ID). Defaults to `"i_NUMBER"`.
#' @param starting.values Numeric vector of starting values for the parameters.
#'   Length must match \code{X}.
#' @param fixed Optional integer vector of indices of parameters that should be
#'   held fixed at their starting values (passed to \code{maxLik}).
#' @param hessian Logical. If \code{TRUE}, compute the final Hessian (slower, but
#'   needed for standard errors). If \code{FALSE}, skip Hessian computation
#'   (useful for initial optimization or many repetitions).
#'
#' @return An object of class \code{"maxLik"} as returned by \code{maxLik::maxLik}
#'   or \code{maxLik::maxBFGS}, containing coefficient estimates, log-likelihood,
#'   and (optionally) the Hessian.
#'
#' @examples
#' \dontrun{
#'   # Suppose SRS is a data frame with columns:
#'   #   i_NUMBER  (individual ID)
#'   #   choice    (1 for chosen alternative, 0 otherwise)
#'   #   x_1, ..., x_5 (explanatory variables)
#'   fit <- mnl_fit(
#'     Data = SRS,
#'     X    = c("x_1", "x_2", "x_3", "x_4", "x_5"),
#'     starting.values = rep(0, 5),
#'     hessian = FALSE
#'   )
#' }
#' @export
mnl_fit <- function(Data,
                    X,
                    choiceVar       = "choice",
                    idVar           = "i_NUMBER",
                    starting.values = NULL,
                    fixed           = NULL,
                    hessian         = TRUE) {

  if (is.null(starting.values)) {
    starting.values <- rep(0, length(X))
  }

  # Ensure Data is a plain data.frame for maxLik
  Data_df  <- as.data.frame(Data)
  GroupID  <- as.factor(Data_df[[idVar]])

  if (hessian) {
    # BHHH with final Hessian (your MNL_st.err version)
    res <- maxLik::maxLik(
      logLik = logLik_mnl,
      Data   = Data_df,
      X      = X,
      choiceVar = choiceVar,
      GroupID   = GroupID,
      start = starting.values,
      print.level = 0,
      method = "BHHH",
      grad  = NULL,
      hess  = NULL,
      fixed = fixed,
      iterlim  = 200,
      constraints = NULL,
      tol    = 1e-6,
      reltol = 1e-6,
      finalHessian = TRUE
    )
  } else {
    # BFGS without Hessian (your MNL_no_st.err version)
    res <- maxLik::maxBFGS(
      fn = logLik_mnl,
      Data   = Data_df,
      X      = X,
      choiceVar = choiceVar,
      GroupID   = GroupID,
      start = starting.values,
      grad  = NULL,
      hess  = NULL,
      fixed = fixed,
      print.level = 0,
      iterlim     = 100,
      constraints = NULL,
      tol    = 1e-8,
      reltol = 1e-8,
      finalHessian = FALSE,
      parscale = rep(1, length(starting.values))
    )
  }

  return(res)
}

#' Fit MNL in two steps (fast optimization + final Hessian)
#'
#' Convenience wrapper that first estimates the MNL model without computing
#' the Hessian and then re-estimates it using the previous estimates as
#' starting values while computing the Hessian for standard errors.
#'
#' @inheritParams mnl_fit
#'
#' @return An object of class \code{"maxLik"} as returned by \code{mnl_fit}
#'   with the Hessian evaluated at the final estimates.
#' @export
mnl_full <- function(Data,
                     X,
                     choiceVar       = "choice",
                     idVar           = "i_NUMBER",
                     starting.values = NULL,
                     fixed           = NULL) {

  # First step: no Hessian
  fit1 <- mnl_fit(
    Data            = Data,
    X               = X,
    choiceVar       = choiceVar,
    idVar           = idVar,
    starting.values = starting.values,
    fixed           = fixed,
    hessian         = FALSE
  )

  # Second step: Hessian at optimum
  fit2 <- mnl_fit(
    Data            = Data,
    X               = X,
    choiceVar       = choiceVar,
    idVar           = idVar,
    starting.values = stats::coef(fit1),
    fixed           = fixed,
    hessian         = TRUE
  )

  fit2
}
