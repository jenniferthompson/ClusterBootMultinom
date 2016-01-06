#' Calculate and plot odds ratios, 95% CIs and Wald p-values for bootstrapped multinomial
#' models
#'
#' From a matrix of bootstrapped coefficients, calculate and plot odds ratios, 95% CIs,
#' and Wald p-values for all coefficients, then plot.
#'
#' @param coef.list List of matrices including coefficients from bootstrapped
#'   models (columns = coefficients).
#' @param label.data If desired, data frame with two columns, variable and var.label,
#'   containing variable names and strings to use in plot labels, respectively. Default is NULL.
#' @param remove.vars Character vector of variable names to **not** include in calculations/plots.
#' Defaults to NULL (show all variables). Passed to \code{get.or.results}.
#' @param round.vars Character vector of variable names whose results should be rounded to something
#' other than two decimal places. Useful for variables with very small changes in odds for one-unit
#' change in variable. Defaults to NULL. Passed to \code{get.or.results}.
#' @param round.digits Integer; number of digits to round [round.vars] to.
#'   Passed to \code{get.or.results}.
#' @param out.strings.list List of character vectors to label outcome comparisons.
#' @param delete.row Row to delete from plots and calculations. Used in situations where
#'   models are run twice with different reference levels; in this case, one comparison is
#'   redundant (eg, 'B vs. A' is reciprocal of 'A vs. B').
#' @param yval.offset Numeric; amount to offset lines for each outcome level in final plot.
#' @return List of 1) \code{or.data}, a data frame containing odds ratios, confidence
#'   limits, p-values and accompanying information; 2) \code{or.plot}, a ggplot2 object
#'   which plots ORs and CIs for all variables and outcome comparisons included, adding
#'   p-values to axis labels.
#' @seealso ggplot2.
#' @import ggplot2
#' @export
#' @examples
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
#' ## Fit model to original and bootstrapped data frame,
#' ##   saving errors and warnings to .txt file
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
#' ## Calculate and plot odds ratios and CIs, Wald p-values for x2
#' covariate.ors <- multi.plot.ors(coef.list = list(boot.matrix.a),
#'                                 out.strings = list(c('B vs A', 'C vs A')),
#'                                 remove.vars = 'x1')

multi.plot.ors <-
  function(coef.list,               ## list of bootstrapped coef matrices
           label.data = NULL,       ## data set to use for variable labels
           remove.vars = NULL,      ## variables to *not* include in plots
           round.vars = NULL,       ## variables to round to digits other than 2
           round.digits = NULL,     ## number of places to round "special" variables to
           out.strings.list,        ## list of character vectors describing outcome comparisons for each element of coef.list
           delete.row = 'none',     ## outcome comparison that is redundant
           yval.offset = 0.25){     ## amount to offset lines for each outcome level in final plot

  ## Calculate ORs, get results strings for plots
  if(is.null(out.strings.list) | (length(coef.list) != length(out.strings.list))){
    stop('Error: length of outcome strings list (out.strings.list) must match number of coefficient matrices (coef.list)')
  } else{
    ## Calculate odds ratios, CIs for each matrix of coefficients using get.or.results()
    ors.list <- lapply(1:length(coef.list), FUN = function(x){
      get.or.results(coef.list[[x]], remove.vars, round.vars, round.digits, out.strings.list[[x]])
    })

    ## Calculate p-values for each matrix of coefficients using Wald test
    ## Currently does not handle chunk tests for overall effect of nonlinear terms/interactions
    p.list <- lapply(coef.list, FUN = function(x){
      apply(x,
            MARGIN = 2,
            FUN = function(x){ aod::wald.test(var(x), mean(x), Terms = 1)$result$chi2['P'] })
    })

    ## Combine ORs, p-values for each set of bootstraps
    plot.data.list <- lapply(1:length(coef.list), FUN = function(x){
      or.data <- ors.list[[x]]
      or.data$variable <- rownames(or.data)
      p.data <- data.frame(variable = names(p.list[[x]]), pvalue = p.list[[x]])
      merge(or.data, p.data, by = 'variable', all.x = TRUE, all.y = FALSE)
    })

    ## Combine data sets, removing redundant row (reciprocal of another row)
    plot.data <- filter(bind_rows(plot.data.list), outlabel != delete.row) %>%
      mutate(## outlabel must be a factor for ggplot
             outlabel = factor(as.character(outlabel)),
             ## Create string for each row for printing numeric results (OR, CI, P) on Y axis
             results = ifelse(pvalue < 0.0001, paste0(outlabel, ': ', results, '; p <0.0001'),
                       ifelse(pvalue < 0.001, paste0(outlabel, ': ', results, '; p <0.001'),
                              paste0(outlabel, ': ', results, '; p = ',
                                     format(round(pvalue, 3), nsmall = 3)))),
             ## Get variable label for more clarity
             varlabel = as.character(unlist(lapply(variable, FUN = function(x){
               varname <- strsplit(x, ':')[[1]][1]
               if(is.null(label.data)){
                 return(varname)
               } else{
                 return(label.data$var.label[match(varname, label.data$variable)])
               }
             }))),
             ## Order variables alphabetically by label
             out.order = unlist(lapply(outlabel, FUN = function(x){
               grep(x, setdiff(unlist(out.strings.list), delete.row)) })))

    ## Set Y value - start at integer, add 0.25 for each outcome comparison
    plot.data$yval <- unlist(lapply(1:nrow(plot.data), FUN = function(x){
      match(plot.data$varlabel[x], sort(unique(plot.data$varlabel))) +
        (plot.data$out.order[x] - 1) * -yval.offset }))

    ## Create vector of strings for Y axis
    plot.ytext <- unlist(lapply(sort(unique(plot.data$varlabel)), FUN = function(x){
      tmp <- plot.data[plot.data$varlabel == x,]
      paste(c(x, tmp$results), collapse = '\n')
    }))

    ## Create plot of odds ratios for non-biomarker covariates for one-unit increase ##
    plot.ors <- ggplot(aes(x = estimate, y = yval), data = plot.data) +
      geom_vline(xintercept = 1, linetype = 'dashed', colour = 'grey70') +
      geom_segment(aes(x = lcl, xend = ucl, y = yval, yend = yval,
                       colour = outlabel, width = 1.25)) +
      geom_point() +
      scale_colour_discrete(guide = FALSE) +
      scale_x_continuous(name = 'Odds Ratio (95% CI) for 1-Unit Increase') +
      scale_y_continuous(name = '',
                         breaks = sort(unique(ceiling(plot.data$yval))),
                         labels = plot.ytext) +
      theme(axis.text.y = element_text(vjust = 0.8))

    return(list('or.data' = plot.data,
                'or.plot' = plot.ors))
  }
}
