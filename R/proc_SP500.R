#' Preprocess SP500 data
#'
#' The \code{proc_SP500} function preprocesses the SP500 stock data by calculating monthly
#' returns for selected sectors and generating idiosyncratic errors.
#'
#' @param SP500data A data frame containing SP500 stock data with columns including:
#'   \describe{
#'     \item{symbol}{Stock symbol.}
#'     \item{date}{Date of the stock data.}
#'     \item{adjusted}{Adjusted closing price of the stock.}
#'     \item{sector}{Sector of the stock.}
#'   }
#' @param sectors A character vector specifying the sectors to include in the analysis.
#' @return A list with components:
#'   \item{Uhat}{A matrix of idiosyncratic errors.}
#'   \item{Khat}{Estimated number of factors.}
#'   \item{factorparthat}{Estimated factor returns.}
#'   \item{sectornames}{Sector for each column of \code{Uhat}.}
#'
#' @details
#' \enumerate{
#'    \item Calculates monthly returns for each stock in the specified sectors  
#'    \item Estimates the number of factors via \code{hdbinseg::get.factor.model(ic="ah")}  
#'    \item Uses \code{POET::POET()} to extract factor loadings/factors and form idiosyncratic errors  
#' }
#'
#' @importFrom dplyr arrange group_by filter select distinct left_join
#' @importFrom purrr map
#' @export
#' @examples
#'
#' data("SP500")
#' set.seed(1234)
#' sectors <- c("Energy", "Financials", "Materials", "Real Estate", "Utilities", "Information Technology")
#' Uhat <- proc_SP500(SP500, sectors)$Uhat
#' \donttest{
#' PPPres <- thresPPP(Uhat, eps = 0, thres = list(value = 0.0020, fun = "hard"), nsample = 100)
#' postmean <- estimate(PPPres)
#' diag(postmean) <- 0 # hide color for diagonal
#' plot(postmean)}
#'
proc_SP500 <- function(SP500data, sectors) {
  preproc_stock <- function(data, period = "yearly") {
    tmo_returns <- data %>%
      dplyr::arrange(sector) %>%
      dplyr::group_by(symbol) %>%
      tidyquant::tq_transmute(
        select = adjusted, mutate_fun = periodReturn,
        period = period, col_rename = "returns"
      ) %>%
      tidyr::spread(key = symbol, value = returns) %>%
      timetk::tk_xts(date_var = date)


    tmo_returns[, (as.data.frame(tmo_returns) %>%
      purrr::map(~ sum(is.na(.x))) %>% unlist()) == 0]
  }
  generateUhat <- function(Ymat, Khat = NULL) {
    if (is.null(Khat)) Khat <- hdbinseg::get.factor.model(Ymat, ic = "ah")$r.hat


    POETest <- POET::POET(t(Ymat), K = Khat, C = 0.1, thres = "hard")
    factorparthat <- t(POETest$loadings %*% POETest$factors)

    list(Uhat = Ymat - factorparthat, Khat = Khat, factorparthat = factorparthat)
  }

  data <- SP500data %>%
    na.omit() %>%
    dplyr::filter(sector %in% sectors)

  Ystock <- data %>%
    split(.$sector) %>%
    purrr::map(~ preproc_stock(.x, "monthly")) %>%
    do.call("cbind", .)

  sector_df <- data %>%
    dplyr::select(symbol, sector) %>%
    dplyr::distinct() %>%
    dplyr::arrange(sector) %>%
    dplyr::mutate(symbol = as.character(symbol)) %>%
    dplyr::mutate(sector = as.character(sector))
  data.frame(symbol = colnames(Ystock)) %>%
    dplyr::left_join(sector_df, by = "symbol") %>%
    .$sector -> sectornames

  # image(cov(Ystock))
  UKhat <- generateUhat(scale(Ystock, scale = FALSE, center = T))
  UKhat$sectornames <- sectornames
  return(UKhat)
}
