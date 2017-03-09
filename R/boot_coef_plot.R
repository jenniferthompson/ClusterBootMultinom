#' Plot Distribution of Bootstrapped Coefficients
#'
#' Create histograms of bootstrapped estimates for all coefficients in a model, adding
#' reference lines at values of coefficients from model fit to original data.
#'
#' @param coef.matrix Matrix of bootstrapped coefficients (rows = data sets, columns =
#'   coefficients).
#' @param org.coefs Numeric vector to plot as reference, usually coefficients from model fit on
#'   original data.
#' @param plot.ints Whether to plot intercept terms. Defaults to FALSE.
#' @return List of data frame of all estimates and final plot, faceted by coefficient.
#' @import ggplot2 dplyr
#' @export
#' @seealso Uses ggplot2 and dplyr.
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
#' ## Check distribution of bootstrapped estimates
#' coefplot.a <- boot.coef.plot(coef.matrix = boot.matrix.a,
#'                              org.coefs = boot.fits.a$org.model@@coefficients,
#'                              plot.ints = FALSE)
#' coefplot.a$coef.plot + ggplot2::ggtitle('Reference = A')

boot.coef.plot <- function(coef.matrix,         ## matrix of bootstrapped coefficients (columns = coef)
                           org.coefs,           ## vector of coefficients from original model
                           plot.ints = FALSE){  ## whether to plot intercept terms

  if(!plot.ints){
    keep.cols <- setdiff(1:ncol(coef.matrix), grep('Intercept', colnames(coef.matrix)))
    coef.matrix <- coef.matrix[,keep.cols]
    org.coefs <- org.coefs[keep.cols]
  }

  ## Create data set in format for ggplot faceting
  coef.data <-
    bind_rows(lapply(1:ncol(coef.matrix),
                     FUN = function(x){data.frame(variable = rep(colnames(coef.matrix)[x],
                                                                 nrow(coef.matrix)),
                                                  bootstrapped = coef.matrix[,x],
                                                  original = org.coefs[x],
                                                  mean.boot = mean(coef.matrix[,x]))} ))
  coef.data$variable <- as.factor(coef.data$variable)

  ## Plot distributions of all coefficients, including intercept, using facets
  coef.plot <- ggplot(aes(x = bootstrapped), data = coef.data) +
    facet_wrap( ~ variable) +
    geom_histogram(colour = 'grey50', fill = 'grey50') +
    ## Reference lines for original model coefficients (red),
    ##  mean(bootstrapped coefficients) (blue)
    geom_vline(aes(xintercept = original), colour = 'red', linetype = 'solid', alpha = 0.5) +
    geom_vline(aes(xintercept = mean.boot), colour = 'blue', linetype = 'dashed', alpha = 0.5) +
    theme(strip.text = element_text(size = 7))

  return(list('coef.data' = coef.data, 'coef.plot' = coef.plot))
}
