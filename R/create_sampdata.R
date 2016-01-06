#' Create Cluster Bootstrapped Data Sets
#'
#' Using an original data frame with N subjects, creates data sets with all records from each of
#' N subject IDs sampled with replacement. Final data sets are saved in a list.
#'
#' @param org.data Original data set, of class \code{data.frame}.
#' @param id.var Character string; name of subject identifier variable to sample.
#' @param n.sets Integer; number of final data sets desired.
#' @param first.seed Numeric (defaults to 56); set for reproducibility.
#' @return List of \code{n.sets} data frames.
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

create.sampdata <- function(org.data,         ## original data set
                            id.var,           ## name of ID variable in original data set
                            n.sets,           ## number of data sets needed
                            first.seed = 56){ ## seed to start with to ensure reproducibility

  ## Get vector of all unique IDs in data set
  all.ids <- unique(org.data[,id.var])

  ## Create vector of seeds for every bootstrap
  set.seed(first.seed)
  seed.vec <- sample(0:100000, size = n.sets)

  ## For each seed, create a data set with all records from each patient in sample with replacement
  ##  using that seed
  data.list <- lapply(seed.vec, FUN = function(x){
    ## Set initial seed
    set.seed(x)

    ## Create vector of sampled IDs
    samp.ids <- sample(all.ids, size = length(all.ids), replace = TRUE)

    final.data <- bind_rows(lapply(samp.ids, FUN = function(y){
      org.data[org.data[,id.var] == y,]
    }))

    return(final.data)
  })

  return(data.list)
}
