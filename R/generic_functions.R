#'
#' Get data from mgcviz smooths.
#'
#' @description Generic function for producing QQ-plots.
#' @export get_data
get_data <- function(...) UseMethod("get_data")


###### Internal generics
.prepare <- function(...) UseMethod(".prepare")
