#' A \code{gee_screenr object}
#'
#' The result of
#' \preformatted{
#'  set.seed(20220401)
#'  uniclus <- unicorns %>%
#'  mutate(cluster = sample(1:25, size = dim(unicorns)[1], replace = TRUE ))
#'  uniobj3 <- gee_screenr(testresult ~ Q1 + Q2 + Q3 + Q5 + Q6 + Q7, id = cluster,
#'                      data = uniclus, link = "logit", Nfolds = 10, seed = 123)
#' uniobj3 <- gee_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7, data = uniclus, link = "logit", Nfolds = 10, seed = 123)
#' }
#' @docType data
#'
#' @format An object of class \code{gee_screenr}
#'
#' @examples
#' \dontrun{
#' summary(uniobj3)
#'}
"uniobj3"
