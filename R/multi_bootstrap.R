#' Run Multinomial Regression Models on a List of Data Sets
#'
#' Run \code{vglm()} on a master data set and each of a list of data sets created from it (often
#' using \code{create.sampdata()}), saving original model object, all successful model objects,
#' and number of models with errors.
#'
#' @param org.data Original data set, of class \code{data.frame}.
#' @param data.sets List of data sets with same variables as \code{org.data}.
#' @param ref.outcome Integer; level of outcome variable to use as reference.
#' @param multi.form Formula used for all models.
#' @param n.boot Integer representing number of successful model fits required. Defaults to
#'   80\% of length of \code{data.sets}.
#' @param xvar String to include in printed status updates. Defaults to "Exposure."
#' @return List of 1) \code{org.model} (model fit to \code{org.data}); 2) \code{boot.models}
#'   (fits for all successful models); and 3) \code{num.failed} (number of models which failed).
#' @seealso \code{\link[VGAM]vglm()}.
#' @import VGAM
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

multi.bootstrap <- function(org.data,                          ## original data set
                            data.sets,                         ## list of data sets for bootstrapping
                            ref.outcome,                       ## outcome level to use as reference
                            multi.form,                        ## model formula
                            n.boot = length(data.sets) / 1.25, ## # successful bootstraps desired
                            xvar = 'Exposure'){                ## name of exposure for status updates

  ## Run model using original data, without accounting for repeated measures;
  ##  if this doesn't converge, stop function
  bio.multi <- do.call(try.vglm,
                       list(formula = multi.form,
                            data = org.data,
                            family = multinomial(refLevel = ref.outcome)))

  if(inherits(bio.multi, 'try-error')){
    stop('Error: Original model does not converge; rethink approach')
  } else{
    ## If original model does converge:
    ## Initialize list of successful model runs, counts of convergence successes and failures
    mod.list <- vector('list', n.boot)
    n.succ <- n.fail <- 0

    ## For each replication...
    iter <- 1
    while(n.succ < n.boot){
      ## Run model on sampled patients
      cur.model <- do.call(try.vglm,
                           list(formula = multi.form,
                                data = data.sets[[iter]],
                                family = multinomial(refLevel = ref.outcome)))

      ## Did model have an error/warning?
      curmod.failed <- inherits(cur.model, 'try-error')

      ## If model failed, increment number of failures and save error messages
      if(curmod.failed){
        n.fail <- n.fail + 1
        ## With first failure, begin a text file to store error messages; with following failues,
        ##  add subsequent errors to it
        if(n.fail == 1){
          sink(file = paste('vglm_errors_', xvar, '.txt', sep = ''), append = FALSE)
          cat('Errors and warnings for ', xvar, '\n\nFailure ', n.fail, ': ', cur.model[1], '\n',
              sep = '')
          sink()
        } else{
          sink(file = paste('vglm_errors_', xvar, '.txt', sep = ''), append = TRUE)
          cat('Failure ', n.fail, ': ', cur.model[1], '\n', sep = '')
          sink()
        }
        ## Otherwise, increment # successes, save model fit, coefficients, SEs, LR test if needed
      } else{
        n.succ <- n.succ + 1
        mod.list[[n.succ]] <- cur.model
      }

      print(paste("Finished bootstrap iteration", iter, "for", xvar))
      iter <- iter + 1
    }

    return(list('org.model' = bio.multi,  ## Original model fit - does not account for repeats
                'boot.models' = mod.list, ## list of all bootstrapped model objects
                'num.failed' = n.fail))   ## number of failures to get to n.boot successful runs
  }
}
