#'
#' Returns standard Laplace variates.
#'
#' \code{rlaplace} function generates a random sample of size n from a standard Laplace distribution.
#'
#' @param n size of the generated sample.
#' @return vector of length n of standard Laplace variates.
#'
#'
rlaplace <-
function(n){
   x = rexp(n)
   ind = 2*(rbinom(n,1,.5)-.5)
   rlaplace = x*ind
   rlaplace
}
