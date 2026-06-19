#' Pumping test information from Bulanghe Basin, China
#'
#' Details of 12 pumping tests conducted at the Bulanghe Basin by the Xi'an
#' Geological Survey, China.  Used as reference data for hydraulic tomography
#' inversion.
#'
#' @format A data frame with 12 rows and 7 variables:
#' \describe{
#'   \item{test_name}{character, the pumping test name (e.g. \code{"A1-1"}).}
#'   \item{pumping_well}{character, the pumping well name (A1 to A6).}
#'   \item{Q_(L_per_s)}{numeric, the pumping rate in L/s \code{[L\eqn{^3}/T]}.}
#'   \item{X}{numeric, the x coordinate of the pumping well \code{[L]}.}
#'   \item{Y}{numeric, the y coordinate of the pumping well \code{[L]}.}
#'   \item{start}{POSIXct, the start time of the pumping test.}
#'   \item{end}{POSIXct, the end time of the pumping test.}
#' }
"pumpingtest_information"

#' Observation well data from Bulanghe Basin, China
#'
#' Steady-state drawdown measurements at 30 observation wells during pumping
#' tests at the Bulanghe Basin, conducted by the Xi'an Geological Survey, China.
#'
#' @format A data frame with 30 rows and 7 variables:
#' \describe{
#'   \item{obs_well}{character, the observation well name.}
#'   \item{X}{numeric, the x coordinate of the observation well \code{[L]}.}
#'   \item{Y}{numeric, the y coordinate of the observation well \code{[L]}.}
#'   \item{value}{numeric, the steady-state drawdown \code{[L]} (originally in cm;
#'     convert to m by dividing by 100).}
#'   \item{pumping_test}{character, the corresponding pumping test name.}
#'   \item{element_number}{integer, the grid element number in an 80×80 grid.}
#'   \item{element_pump}{integer, the grid element number for the pumping well
#'     in an 80×80 grid.}
#' }
"dfl_ss_unique"

#' Pumping test information in list form
#'
#' A list of 12 data frames, one per pumping test.  Each data frame contains
#' the pumping rate and well location.  This is the format expected by all
#' forward and inverse functions.
#'
#' @format A list of 12 data frames, each with columns:
#' \describe{
#'   \item{Qp}{numeric, the pumping rate in m\eqn{^3}/d \code{[L\eqn{^3}/T]}.}
#'   \item{x}{numeric, the x coordinate of the pumping well \code{[L]}.}
#'   \item{y}{numeric, the y coordinate of the pumping well \code{[L]}.}
#' }
"qHT"

#' Observation data in list form
#'
#' A list of 12 data frames, one per pumping test.  Each data frame contains
#' the observed drawdown and observation well location.  This is the format
#' expected by all forward and inverse functions.
#'
#' @format A list of 12 data frames, each with columns:
#' \describe{
#'   \item{data}{numeric, the observed drawdown \code{[L]}.}
#'   \item{x}{numeric, the x coordinate of the observation well \code{[L]}.}
#'   \item{y}{numeric, the y coordinate of the observation well \code{[L]}.}
#' }
"oHT"

#' Synthetic data for the Bulanghe Basin test case
#'
#' A synthetic test dataset generated using the same pumping and observation
#' well configuration as the real Bulanghe Basin data, with a known true
#' transmissivity field for validating inversion results.
#'
#' @format A list with components:
#' \describe{
#'   \item{trueK}{data frame, the true transmissivity field (columns \code{x},
#'     \code{y}, \code{Tp}) \code{[L\eqn{^2}/T]}.}
#'   \item{qHT1}{list of 6 pumping test data frames (columns \code{Qp, x, y}).}
#'   \item{oHT1}{list of 6 observation data frames (columns \code{data, x, y}).}
#'   \item{grid}{the grid list from \code{GenGrid()}.}
#' }
"blh_synthetic"

#' Best-truth synthetic data for the Bulanghe Basin test case
#'
#' A refined synthetic dataset with the best-estimate true transmissivity field
#' for validating inversion results against the Bulanghe Basin configuration.
#'
#' @format A list with components:
#' \describe{
#'   \item{trueK}{data frame, the true transmissivity field (columns \code{x},
#'     \code{y}, \code{Tp}) \code{[L\eqn{^2}/T]}.}
#'   \item{qHT1}{list of 6 pumping test data frames (columns \code{Qp, x, y}).}
#'   \item{oHT1}{list of 6 observation data frames (columns \code{data, x, y}).}
#'   \item{grid}{the grid list from \code{GenGrid()}.}
#' }
"blh_syn_best_truth"
