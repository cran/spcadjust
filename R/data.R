#' Cardiac surgery data
#'
#' A dataset describing the results of cardiac surgery. The data give
#' information about the date of surgery, surgeon, Parsonnet score and
#' outcome up to 90 days after surgery.
# 
#
#'
#' @format A data frame with  5595 rows and 5 variables:
#' \tabular{ll}{
#'   date:\tab date of the operation in days since the beginning of study \cr
#'   time:\tab number of days between surgery and the earlier of death and 90 days\cr
#'   status:\tab status at endpoint, 0/1 for censored/dead\cr
#'   Parsonnet:\tab Parsonnet score of the patient\cr
#'   surgeon:\tab surgeon performing the operation, numbered from 1 to 7\cr
#' }
#' @keywords datasets
#'
#' @usage data(cardiacsurgery)
#' 
#' @source Based on the data described in Steiner et al (2000). A
#' subset of the data has been selected, some noise has been
#' introduced and the follow-up was censored at 90 days.
#' 
#' @references Steiner SH, Cook RJ, Farewell VT, Treasure T (2000).
#' Monitoring surgical performance using risk-adjusted cumulative sum
#' charts. Biostat 1(4) 441-452.
"cardiacsurgery"
