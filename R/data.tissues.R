#' tissues dataset
#'
#' The tissues data of colon cancer dataset,
#' which includes gene expression values from 22 colon tumor tissues and 40 non-tumor tissues.
#' @name tissues
#' @format `numeric`
#' @docType data
#' @source http://genomics-pubs.princeton.edu/oncology/affydata/.
#' @keywords data
#' @examples
#'
#' \dontrun{
#' tissues <- scan("tissues.txt")
#' devtools::use_data(tissues, overwrite = TRUE)
#' data("tissues")
#' }
NULL