#' Odds Ratios for Cluster Bootstrapped Multinomial Logistic Regression
#'
#' From a matrix of coefficient estimates from B vglm() model fits on cluster bootstrapped data,
#' calculate and plot odds ratios and 95% CIs for one-unit change in each variable. Written to be
#' used within multi.plot.ors(), but can be called explicitly. Note: does not currently support
#' group effects of variables with nonlinear terms; will calculate one odds ratio per coefficient.
#'
#' @param coef.matrix Matrix of coefficient estimates, with rows = number of bootstrapped data sets,
#' columns = number of coefficients.
#' @param remove.vars Character vector of variable names to **not** include in calculations/plots.
#' Defaults to NULL (show all variables).
#' @param round.vars Character vector of variable names whose results should be rounded to something
#' other than two decimal places. Useful for variables with very small changes in odds for one-unit
#' change in variable. Defaults to NULL.
#' @param round.digits Integer; number of digits to round [round.vars] to.
#' @param out.strings List of character strings to label outcome comparisons. Defaults to B vs. A,
#' C vs. A, etc. Note that onus is on the user to supply correct labels.
#' @return Data frame with one record per coefficient, including odds ratio estimate, lower and
#' upper confidence limits, character string of results (format: "OR (LCL, UCL)") and text
#' describing comparison.
#' @seealso \code{\link[VGAM]{vglm}}, which this function assumes you are using;
#' \code{multi.plot.ors}, which calculates p-values and plots all results using ggplot2.
#' @export
#' @examples
#' ## Create data frame
#' df <- data.frame(id = sample(1:20, size = 100, replace = TRUE),
#'                  x1 = rnorm(n = 100),
#'                  x2 = rbinom(p = 0.75, n = 100, size = 1),
#'                  y = sample(LETTERS[1:3], size = 100, replace = TRUE))
#' df <- df[order(df$id),]
#' df$time <- unlist(lapply(1:length(unique(df$id)),
#'                          FUN = function(idnum){ 1:nrow(df[df$id == unique(df$id)[idnum],]) }))
#'
#' ## Using create.sampdata(), generate list of cluster bootstrapped data sets
#' bootdata.list <- create.sampdata(org.data = df,
#'                                  id.var = 'id',
#'                                  n.sets = 25)
#'
#' ## Fit model to original and bootstrapped data frame, saving errors and warnings to .txt file
#' boot.fits.a <- multi.bootstrap(org.data = df,
#'                                data.sets = bootdata.list,
#'                                ref.outcome = grep('A', levels(df$y)),
#'                                multi.form = as.formula('y ~ x1 + x2'))
#'
#' ## Create matrices of coefficients for all bootstrap fits
#' boot.matrix.a <- do.call(rbind,
#'                          lapply(boot.fits.a$boot.models,
#'                                 FUN = function(x){ x@@coefficients }))
#'
#' ## Get odds ratios and CIs for x2
#' ors <- get.or.results(boot.matrix.a, remove.vars = 'x1')

get.or.results <- function(coef.matrix,          ## matrix of bootstrapped model coefficients
                           remove.vars = NULL,   ## variables NOT to calculate ORs for
                           round.vars = NULL,    ## variable names to round to more places
                           round.digits = NULL,  ## number of decimal places to round to if not 2
                           out.strings = NULL){  ## labels for each set of coefficients

  ## Remove columns for intercepts and any other variables requested
  remove.cols <- grep(paste0(c('Intercept', remove.vars), collapse = '|'), colnames(coef.matrix))
  coef.matrix <- coef.matrix[,setdiff(1:ncol(coef.matrix), remove.cols)]

  ## Get final bootstrapped coefficients and CI for each covariate
  ##  (mean, 2.5th & 97.5th percentiles of all bootstrapped coefficients)
  coefs <- colMeans(coef.matrix, na.rm = TRUE)
  cis <- apply(coef.matrix, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  coefs.cis <- rbind(coefs, cis)
  rownames(coefs.cis) <- c('estimate', 'lcl', 'ucl')

  ## Transpose to rows = variables, columns = quantities and exponentiate to get ORs/CIs
  ors <- as.data.frame(t(exp(coefs.cis)))

  ## Create text string of results for each variable, rounded to round.digits places if requested
  ors$round.places <- 2
  if(!is.null(round.vars)){
    ors$round.places[grep(paste0(round.vars, collapse = '|'), rownames(ors))] <- round.digits
  }
  ors$results <- unlist(lapply(1:nrow(ors), FUN = function(x){
    with(ors[x,], {
      paste0(format(round(estimate, round.places), nsmall = round.places), ' (',
             format(round(lcl, round.places), nsmall = round.places), ', ',
             format(round(ucl, round.places), nsmall = round.places), ')')
    })
  }))

  ## Add label for which outcome is which
  ## Get coefficient set for each row
  ors$coefset <-
    as.numeric(unlist(lapply(rownames(ors), FUN = function(x){ strsplit(x, ':')[[1]][2] })))

  ## Get number of unique sets of coefficients
  n.coefsets <- length(unique(ors$coefset))

  ## Does number of out.strings match number of unique sets of coefficients?
  if(!is.null(out.strings) & length(out.strings) != n.coefsets){
    stop('Error: number of outcome labels must match number of unique sets of coefficients')
  } else{
    ## If no out.strings argument provided, create defaults of B, C, ... vs. A
    if(is.null(out.strings)){
      out.strings <- paste(LETTERS[2:(n.coefsets + 1)], 'vs. A')
    }
  }

  ors$outlabel <- out.strings[ors$coefset]

  return(subset(ors, select = -c(round.places, coefset)))
}
