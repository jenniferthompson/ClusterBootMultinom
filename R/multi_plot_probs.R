#' Calculate and Plot Predicted Probabilities and 95\% Confidence Intervals for Each
#' Outcome Level vs. Continuous Covariate
#'
#' Using matrices of bootstrapped multinomial logistic regression coefficients, calculate
#' and plot predicted probabilities and 95\% confidence limits for each outcome at a
#' given level of a continuous covariate. Also includes Wald p-values for each comparison.
#' This function can handle covariates which are either linear or are modeled using
#' restricted cubic splines with three knots. This function is commonly used with two
#' model fits and corresponding coefficient matrices, representing the outcome with the
#' lowest and highest level as reference, respectively. This allows for calculation of
#' confidence limits for all levels of the outcome. Note that in this case, the coefficient
#' matrices, variance-covariance matrices, and model objects should all be specified in a
#' consistent order (eg, reference = lowest, followed by reference = highest).
#'
#' @param xval String; name of variable to plot on X axis.
#' @param xval.limits Numeric vector; percentiles of \code{xval} to include. Defaults to
#'   excluding lower and upper 10\% of X values.
#' @param xval.knots Numeric vector of length 3; percentile values of knot placement if
#'   \code{xval} is modeled with restricted cubic splines. Defaults to c(0.05, 0.50, 0.95) -
#'   \code{rcspline.eval()} defaults.
#' @param data.set Data frame to use for quantiles and label of \code{xval}.
#' @param mod.objs List of \code{vglm()} model objects - original model fit(s).
#' @param design.mats List of numeric vectors representing values all covariates are
#'   adjusted to, corresponding to placement of coefficients. Assumes that \code{xval}
#'   is the last variable in the model and should be placed at the end of these vectors.
#' @param coef.list List of matrices containing bootstrapped coefficient estimates.
#' @param vcov.list List of variance-covariance matrices of coefficient estimates.
#' @param plot.raw Add raw data to plots? Defaults to TRUE.
#' @return List of 1) Data frame of all predicted probabilities and confidence
#'   limits (\code{pp.data}); 2) Data frame of Wald test results for each outcome
#'   comparison (\code{results.data}); 3) ggplot2 object which contains plots of
#'   \code{xval} vs. predicted probability of each outcome level.
#' @seealso \code{\link[VGAM]{vglm}}, which this function assumes you are using;
#'   \code{rcs()}; Hmisc (\code{label()}; aod (\code{wald.test()}); ggplot2; dplyr; tidyr.
#' @import tidyr
#' @importFrom Hmisc label
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
#' ## Calculate predicted probs and CIs for x1 at outcomes B, C
#' ## Design matrix: first two columns = intercepts, second two set X2 = 1
#' design.tmp <- matrix(c(1, 0, 0, 1, 1, 0, 0, 1), nrow = 2)
#'
#' ## Calculate and plot predicted probabilities for outcomes at levels of x1
#' ## To get probabilities and CIs for all outcome levels, run models twice,
#' ##  with lowest and highest outcome levels as reference
#' ## Fit model to original and bootstrapped data frame, saving errors and warnings to .txt file
#' boot.fits.c <- multi.bootstrap(org.data = df,
#'                                data.sets = bootdata.list,
#'                                ref.outcome = grep('C', levels(df$y)),
#'                                multi.form = as.formula('y ~ x1 + x2'))
#'
#' boot.matrix.c <- do.call(rbind,
#'                          lapply(boot.fits.c$boot.models,
#'                                 FUN = function(x){ x@@coefficients }))
#'
#' x1.prob.results <- multi.plot.probs(xval = 'x1',
#'                                     data.set = df,
#'                                     design.mats = list(design.tmp[1,],
#'                                                        design.tmp[2,]),
#'                                     mod.objs = list(boot.fits.a$org.mod,
#'                                                     boot.fits.c$org.mod),
#'                                     coef.list = list(boot.matrix.a,
#'                                                      boot.matrix.c),
#'                                     vcov.list = list(var(boot.matrix.a),
#'                                                      var(boot.matrix.c)))

multi.plot.probs <-
  function(
    xval, ## variable name of variable on X axis
    xval.limits = c(0.1, 0.9), ## quantiles to include in figures; defaults to 10th-90th pctiles
    xval.knots = c(0.05, 0.5, 0.95), ## values of knot placement for xval; defaults are rcspline.eval defaults
    data.set, ## data set to use
    mod.objs, ## list of model objects
    design.mats, ## list of design matrices to start with
    coef.list, ## list of bootstrapped coefficient matrices
    vcov.list, ## list of variance-covariance matrices
    plot.raw = TRUE ## Add raw data to final plot?
  ){

  ## -- Create multinomial design matrix -- ##
  xval.vals <- unique(data.set[!is.na(data.set[,xval]), xval])
  ## Get actual numeric limits for plot X axis
  xval.quant <- quantile(data.set[,xval], probs = xval.limits, na.rm = TRUE)

  ## Subset xval values to those between quantiles of interest
  xval.vals <- xval.vals[xval.vals >= xval.quant[1] & xval.vals <= xval.quant[2]]

  ## If a nonlinear term is involved... ##
  if(length(grep(paste('rcs(', xval, sep = ''),
                 dimnames(mod.objs[[1]]@x)[[2]],
                 fixed = TRUE)) > 0){
    ## Calculate spline terms for unique values of bioxval
    xval.k <- quantile(xval.vals, probs = xval.knots, na.rm = TRUE)
    xval.splines <- calc.spline(xval.vals,
                                  k1 = xval.k[1],
                                  k2 = xval.k[2],
                                  k3 = xval.k[3])

    ## Create design matrices for each outcome: all values the same except xval
    ## Design matrices must be in order by outcome for this to be correct
    ## Also, currently only works for three knots
    xval.design <- do.call(rbind,
                             lapply(1:length(design.mats), FUN = function(designmat){
                               tmp <- tmpspline <- rep(0, length(design.mats))
                               do.call(rbind, lapply(1:length(xval.vals), FUN = function(x){
                                 tmp[designmat] <- xval.vals[x]
                                 tmpspline[designmat] <- xval.splines[x]
                                 return(c(design.mats[[designmat]], tmp, tmpspline))
                               }))
                             }))
  } else{ ## If no nonlinear term involved...
    xval.design <- do.call(rbind,
                             lapply(1:length(design.mats), FUN = function(designmat){
                               tmp <- rep(0, length(design.mats))
                               do.call(rbind, lapply(1:length(xval.vals), FUN = function(x){
                                 tmp[designmat] <- xval.vals[x]
                                 return(c(design.mats[[designmat]], tmp))
                               }))
                             }))
  }

  ## Calculate linear predictors, predicted probabilities, SEs, CIs for all values of bioxval
  ##  for each model object
  pp.list <- lapply(1:length(mod.objs), FUN = function(modnum){
    calc.ppci(design.matrix = xval.design,
              model.obj = mod.objs[[modnum]],
              use.coefs = colMeans(coef.list[[modnum]]),
              use.vcov = vcov.list[[modnum]])
  })

  ## -- Combine results for plotting -- ##
  ## Need results for all combinations of (1 vs. last), (all others vs. 1)
  n.outcomes <- ncol(pp.list[[1]]$PredictedProbs) + 1
  out.combs <- c(paste0('log(mu[,1]/mu[,', n.outcomes, '])'),
                 paste0('log(mu[,', 2:n.outcomes, ']/mu[,1])'))

  ## Get point estimates, CIs for outcome comparisons
  ppest <- do.call(cbind, lapply(pp.list, FUN = function(x){ x$PredictedProbs }))[,out.combs]
  pplcl <- do.call(cbind, lapply(pp.list, FUN = function(x){ x$ProbsLCLs }))[,out.combs]
  ppucl <- do.call(cbind, lapply(pp.list, FUN = function(x){ x$ProbsUCLs }))[,out.combs]

  ## Substitute text of outcome levels for number values
  out.levels <- unlist(lapply(colnames(ppest), FUN = function(x){
      mod.objs[[1]]@misc$ynames[as.numeric(gsub('log(mu[,', '', strsplit(x, ']/')[[1]][1],
                                                   fixed = TRUE))]
      }))
  colnames(ppest) <- paste(out.levels, 'pp', sep = '.')
  colnames(pplcl) <- paste(out.levels, 'lcl', sep = '.')
  colnames(ppucl) <- paste(out.levels, 'ucl', sep = '.')

  ## Combine into single data set with one row per xval value + outcome level
  ppest <- data.frame(cbind(xval.vals, ppest))
  pplcl <- data.frame(cbind(xval.vals, pplcl))
  ppucl <- data.frame(cbind(xval.vals, ppucl))

  pp.data <- left_join(ppest, pplcl, by = 'xval.vals') %>%
    left_join(ppucl, by = 'xval.vals') %>%
    gather(key = out.quant, value = quantity, 2:ncol(.)) %>%
    separate(out.quant, into = c('outcome', 'prob.quant'), sep = '\\.') %>%
    spread(key = prob.quant, value = quantity)

  ## -- Calculate p-values -- ##
  ## First, figure out which comparison each set of coefs in each model is making
  comparisons <- do.call(rbind, lapply(1:length(mod.objs), FUN = function(x){
    prednames <- gsub('[', '(',
                      gsub(']', ')', mod.objs[[x]]@misc$predictors.names, fixed = TRUE),
                      fixed = TRUE)
    tmp <- strsplit(gsub('[a-z(),]', '', prednames), '/')
    do.call(rbind, lapply(1:length(tmp), FUN = function(y){
      data.frame(modobj = x,
                 compnum = y,
                 lev1 = out.levels[as.numeric(tmp[[y]][1])],
                 lev2 = out.levels[as.numeric(tmp[[y]][2])],
                 comp = paste(out.levels[as.numeric(tmp[[y]][1])],
                              out.levels[as.numeric(tmp[[y]][2])],
                              sep = ' vs. '),
                 sorted = paste0(min(out.levels[as.numeric(tmp[[y]][1])],
                                     out.levels[as.numeric(tmp[[y]][2])]),
                                 max(out.levels[as.numeric(tmp[[y]][1])],
                                     out.levels[as.numeric(tmp[[y]][2])])))
    }))
  }))

  ## Keep only unique comparisons - don't need B vs. A from model 2 if A vs. B is in model 1
  comparisons <-
    subset(comparisons, modobj == 1 | !(sorted %in% subset(comparisons, modobj == 1)$sorted))

  ## For each row in comparisons, calculate Wald p-value for bioxval-related term(s)
  wald.vals <- lapply(1:nrow(comparisons), FUN = function(x){
    use.mod <- comparisons[x, 'modobj']
    use.comp <- comparisons[x, 'compnum']

    ## Get which columns of coefficient matrix to test
    use.cols <- grep(paste0(xval, '.*:', use.comp, '$'), colnames(coef.list[[use.mod]]))

    ## Get means: if one term involved, mean of that column; otherwise, colMeans(specified columns)
    col.means <- NULL
    if(length(use.cols) == 1){
      col.means <- mean(coef.list[[use.mod]][,use.cols])
    } else{
      col.means <- colMeans(coef.list[[use.mod]][,use.cols])
    }

    aod::wald.test(var(coef.list[[use.mod]][,use.cols]), col.means, Terms = 1:length(use.cols))
  })

  ## Add chi square, df, p-value to comparisons data set
  comparisons <- cbind(comparisons,
                       do.call(rbind, lapply(wald.vals, FUN = function(x){ x$result$chi2 })))

  ## Create string to place at bottom of plot
  p.text <- paste0('Bootstrapped Wald P-Values:\n',
                   paste(unlist(lapply(1:nrow(comparisons), FUN = function(x){
                     paste0(comparisons[x, 'comp'], ': p',
                            ifelse(comparisons[x, 'P'] < 0.0001, ' < 0.0001',
                            ifelse(comparisons[x, 'P'] < 0.001, ' < 0.001',
                                   paste(' =', format(round(comparisons[x, 'P'], 3), nsmall = 3)))))
                     })),
                     collapse = '\n'))

  ## -- Plot predicted probabilities + CIs -- ##
  ## Is xval labeled in original data set?
  xval.lab <- Hmisc::label(data.set[,xval])
  xval.lab <- ifelse(xval.lab == '', xval, xval.lab)

  xval.probs <- ggplot(aes(x = xval.vals, y = pp), data = pp.data) +
    facet_wrap( ~ outcome) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = 'grey80') +
    geom_line(aes(colour = outcome)) +
    scale_colour_discrete(guide = FALSE) +
    scale_x_continuous(name = paste(xval.lab, p.text, sep = '\n\n')) +
    scale_y_continuous(limits = 0:1,
                       name='Adjusted Probability of Outcome\nAdjusted to Median or Mode of Covariates')

  ## -- Plot raw data -- ##
  if(plot.raw){
    data.raw <-
      data.set[!is.na(data.set[,xval]) &
                 data.set[,xval] >= xval.quant[1] &
                 data.set[,xval] <= xval.quant[2],
               c(xval, as.character(mod.objs[[1]]@terms$terms[1])[2])]
    names(data.raw) <- c('xval', 'outcome')

    xval.probs <- xval.probs +
      geom_point(aes(x = xval, colour = outcome, y = 1), data = data.raw,
                 alpha = 0.2, position = position_jitter(height = 0.1))
  }

  return(list('pp.data' = pp.data,
              'results.data' = comparisons,
              'prob.line.plot' = xval.probs))

}
