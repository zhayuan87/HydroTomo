#' The pumping tests information conducted in Bulanghe Basin, China by Xi'an Geological Survey.
#'
#' @format A data frame with 12 rows and 7 variables:
#' \describe{
#'   \item{test_name}{character, the pumping test name of the 12 tests}
#'   \item{pumping_well}{A1 to A6, the pumping well name.}
#'   \item{Q_(L_per_s)}{numeric, the pumping rate in cubic Liter per second.}
#'   \item{X}{numeric, the x coordinate of the pumping well in meter.}
#'   \item{Y}{numeric, the y coordinate of the pumping well in meter.}
#'   \item{start} {POSIXct, the start time of the pumping test.}
#'   \item{end} {POSIXct, the end time of the pumping test.}
#' }
"pumpingtest_information"

#' The observation wells information conducted in Bulanghe Basin, China by Xi'an Geological Survey.
#' @format A data frame with 30 rows and 5 variables:
#' \describe{
#'   \item{obs_well}{character, the observation well name of the 30 wells}
#'   \item{X}{numeric, the x coordinate of the observation well in meter.}
#'   \item{Y}{numeric, the y coordinate of the observation well in meter.}
#'   \item{value}{numeric, the steady-state drawdown value in centi-meter.}
#'   \item{pumping_test}{character, the pumping test name corresponding to the observation well.}
#'   \item{element_number}{integer, the default element number in a 80*80 element.}
#'   \item{element_pump}{integer, the default element number for pumping well in a 80*80 element.}
#'   }
"dfl_ss_unique"

#' The pumping test information in package form.
#' @format A list of 12 data frames, each data frame contains pumping rate and location.
#' \describe{
#'   \item{A1-1}{data frame of pumping test 1, with columns Qp(m3/d), x, y.}
#'   \item{A1-2}{data frame of pumping test 2, with columns Qp, x, y.}
#'   ...
#'   }
"qHT"

#' The pumping test information in package form.
#' @format A list of 12 data frames, each data frame contains observation data and location.
#' \describe{
#'   \item{A1-1}{data frame of pumping test 1, with columns data(drawdown in m), x, y.}
#'   \item{A1-2}{data frame of pumping test 2, with columns data, x, y.}
#'   ...
#'   }
"oHT"

#' The synthetic data generated using the same blh pumping and observation wells.
#' @format A list of oHT1, qHT1, trueK, and grid.
#' \describe{
#'   \item{trueK}{data frame of true transmissivity field used in the synthetic data generation.}
#'   \item{qHT1}{list of 6 pumping tests data frames used in the synthetic data generation.}
#'   \item{oHT1}{list of 6 observation data frames used in the synthetic data generation.}
#'   \item{grid}{the grid information used in the synthetic data generation.}
#'   }
"blh_synthetic"
