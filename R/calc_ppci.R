#' Calculate Predicted Probabilities and 95\% Confidence Intervals for Each Outcome Level vs.
#' Continuous Covariate
#'
#' Given a matrix of bootstrapped coefficient estimates from multinomial regression using
#' \code{vglm()}, calculates linear predictors, predicted probabilities, and SEs and
#' confidence limits for a continuous exposure variable at all non-reference levels of the
#' outcome. Has the capabilitiy to include restricted cubic splines using \code{rcs()}.
#'
#' @param design.matrix Design matrix with covariate values, stacked by outcome level.
#' @param model.obj \code{vglm()} model object, used to get number of outcome levels and
#'   label linear predictors.
#' @param use.coefs Numeric vector of coefficients, if not taken from original model (eg,
#'   \code{colMeans}(matrix of bootstrapped coefficients).
#' @param use.vcov Numeric matrix to use as variance-covariance matrix (eg,
#'   var(matrix of bootstrapped coefficients)).
#' @return List of linear predictors (LinearPredictors); variance-covariance matrix of
#'   linear predictors (VarianceLP), predicted probabilities (PredictedProbs) and
#'   their variance-covariance matrix (VariancePPs); SEs (ProbsSes), lower and upper
#'   confidence limits of predicted probabilities (ProbsSEs, ProbsLCLs and ProbsUCLs,
#'   respectively). Each has (number of outcome levels - 1) columns, representing the
#'   quantities for all outcome levels except the reference.
#' @seealso \code{\link[VGAM]{vglm}}, which this function assumes you are using;
#' \code{multi.plot.probs}, which calls this function; \code{rcs()}.
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
#' x1.vals <- sort(unique(df$x1))
#'
#' ## Add all unique x1 values to get complete design matrix
#' x1.design <- do.call(rbind,
#'                      lapply(1:nrow(design.tmp), FUN = function(r){
#'                        tmp <- matrix(rep(c(design.tmp[r,],
#'                                            rep(0, (length(unique(df$y)) - 1))),
#'                                          length(x1.vals)),
#'                                      nrow = length(x1.vals), byrow = TRUE)
#'                        tmp[,(ncol(design.tmp) + r)] <- x1.vals
#'                        tmp
#'                      }))
#'
#' ## Calculate linear predictors, predicted probabilities, etc
#' x1probs <- calc.ppci(design.matrix = x1.design,
#'                      model.obj = boot.fits.a$org.model,
#'                      use.coefs = colMeans(boot.matrix.a),
#'                      use.vcov = var(boot.matrix.a))

calc.ppci <- function(design.matrix, ## design matrix with covariate values, stacked by outcome level
                      model.obj, ## vglm() model object (original model fit for bootstrapped)
                      use.coefs, ## vector of coefficients, if not from original model
                      use.vcov){ ## variance-covariance matrix to use, if not from original model

  ## Get number of unique outcomes and number of unique combinations of covariate values
  n.levels <- ncol(model.obj@y) ## Number of possible outcomes
  n.comb <- nrow(design.matrix) / (n.levels - 1) ## Number of combinations

  ## Calculate linear predictors for all combinations
  lin.pred <- apply(design.matrix, MARGIN = 1, FUN = function(x){ sum(use.coefs * as.numeric(x)) })

  ## Split design matrix, one matrix per non-reference outcome level
  list.design <- lapply(1:(n.levels - 1), FUN = function(x){design.matrix[1:n.comb + (x - 1)*n.comb,]})

  ## Create matrix of linear predictors (rows = # design matrices, columns = # unique marker values)
  lp.mat <- matrix(lin.pred, nrow = (n.levels - 1), ncol = n.comb, byrow = TRUE)

  ## Create **variance-covariance matrix for linear predictor** for each combination
  ## - #rows = #columns = number of unique outcomes - 1
  ## x = row in design matrix (signifying a specific combination of covariate values)
  var.lp <- lapply(1:n.comb,
                   FUN = function(c){
                     vlp <- matrix(NA, nrow = (n.levels - 1), ncol = (n.levels - 1))
                     for(x in 1:nrow(vlp)){
                       for(y in 1:ncol(vlp)){
                         vlp[x,y] <-
                           t(as.numeric(list.design[[x]][c,])) %*%
                            use.vcov %*%
                            as.numeric(list.design[[y]][c,])
                       }
                     }
                     return(vlp)
                   })

  ## Create **beta matrix** (per Appendix B of Liu) for each combination
  ## - #rows = #columns = number of unique outcomes - 1
  ## x = row in design matrix (signifying a specific combination of covariate values)
  beta.mat <-
    lapply(1:n.comb,
           FUN = function(c){
             denom <- (1 + sum(exp(lp.mat[,c])))**2

             betamat <- matrix(NA, nrow = (n.levels - 1), ncol = (n.levels - 1))
             for(x in 1:nrow(betamat)){
               for(y in 1:ncol(betamat)){
                 if(x == y){
                   betamat[x,y] <-
                     (exp(lp.mat[x,c]) * (1 + sum(exp(lp.mat[setdiff(1:nrow(lp.mat), x), c])))) / denom
                 } else{
                   betamat[x,y] <- (-exp(lp.mat[x,c])*exp(lp.mat[y,c])) / denom
                 }
               }
             }

             return(betamat)

           })

  ## Calculate variance-covariance matrix for each set of **predicted probabilities**
  var.pp <- lapply(1:n.comb,
                   FUN = function(x){ beta.mat[[x]] %*% var.lp[[x]] %*% t(beta.mat[[x]]) })

  ## Calculate predicted probabilities of each outcome and CIs for each combination
  pred.probs <- do.call(rbind,
                        lapply(1:n.comb,
                               FUN = function(c){
                                 probs <- rep(NA, (n.levels - 1))
                                 for(i in 1:(n.levels - 1)){
                                   probs[i] <- exp(lp.mat[i, c]) / (1 + sum(exp(lp.mat[,c])))
                                 }
                                 return(probs)
                               }))

  pp.ses <- do.call(rbind,
                    lapply(1:n.comb,
                           FUN = function(c){
                             sqrt(diag(var.pp[[c]]))
                           }))

  pp.lcls <- pred.probs - (qnorm(.975)*pp.ses)
  pp.ucls <- pred.probs + (qnorm(.975)*pp.ses)

  ## Truncate at 0, 1 if needed (sometimes necessary with delta method)
  pp.lcls <- ifelse(pp.lcls < 0, 0, pp.lcls)
  pp.ucls <- ifelse(pp.ucls > 1, 1, pp.ucls)

  rownames(pred.probs) <- rownames(pp.ses) <- rownames(pp.lcls) <- rownames(pp.ucls) <-
    colnames(lp.mat) <- names(var.lp) <- names(var.pp) <- paste('Combination', 1:n.comb, sep = '')
  colnames(pred.probs) <- colnames(pp.ses) <- colnames(pp.lcls) <- colnames(pp.ucls) <-
    rownames(lp.mat) <- colnames(summary(model.obj)@pearson.resid)

  ## Return all possibly relevant info
  return(list('LinearPredictors' = lp.mat,
              'VarianceLP' = var.lp,
              'VariancePP' = var.pp,
              'PredictedProbs' = pred.probs,
              'ProbsSEs' = pp.ses,
              'ProbsLCLs' = pp.lcls,
              'ProbsUCLs' = pp.ucls))
}
