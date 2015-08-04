#' Calculate Spline Terms Based on rcs()
#'
#' Calculate X' values for restricted cubic spline term, eg, linear predictor = B0 + B1*X + B2*X'.
#' Based on Frank Harrell's Regression Modeling Strategies.
#'
#' @param x Numeric vector
#' @param k1 Value of knot 1
#' @param k2 Value of knot 2
#' @param k3 Value of knot 3
#' @return Numeric vector of X' terms.
#' @seealso \code{\link[rms]{rcs}} which this function assumes you use,
#' \code{\link[Hmisc]{rcspline.eval}} for getting knot locations
#' @export
#' @examples
#' calc.spline(rnorm(n = 100), k1 = -1.27, k2 = 0.02, k3 = 1.04)

calc.spline <- function(x, k1, k2, k3){
  y1 <- (x - k1)**3
  y2 <- (x - k2)**3
  y3 <- (x - k3)**3
  y <- ifelse(y1 < 0, 0, y1) -
    ifelse(y2 < 0, 0, y2 * ((k3 - k1) / (k3 - k2))) +
    ifelse(y3 < 0, 0, y3 * ((k2 - k1) / (k3 - k2)))
  return(y / (k3 - k1)**2)
}
