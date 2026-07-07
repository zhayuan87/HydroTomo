# =============================================================================
# 2D Transient Groundwater Flow Simulation
#
# Solves  S ∂h/∂t = ∇·(T ∇h) + Q  using deSolve::ode.2D().
#
# Key relationship to steady-state code:
#   - diffusion2D() (steady2D.R) already returns ∂h/∂t = (1/S)·[flux - Q],
#     so it is reused unchanged as the ODE right-hand side.
#   - Ftransient2dsim() replaces steady.2D() with ode.2D() and exposes the
#     storage coefficient S and the time vector as new parameters.
#   - samDataTr() samples the ODE output at arbitrary (x, y, time) triplets.
# =============================================================================


#' Solve 2D transient groundwater flow equation
#'
#' Solves the 2D groundwater flow equation
#' \eqn{S \partial h / \partial t = \nabla \cdot (T \nabla h) + Q}
#' using \code{deSolve::ode.2D()}.  Boundary conditions are zero drawdown
#' on all four edges.  The PDE is discretised by the same finite-difference
#' scheme as the steady-state solver; see \code{\link{diffusion2D}}.
#'
#' @param domain  6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.
#' @param grid    grid from \code{GenGrid()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param TT      transmissivity \eqn{T} L\eqn{^2}/T — scalar or length-n
#'   vector (one value per cell).
#' @param SS      storage coefficient \eqn{S} [-] — scalar or length-n
#'   vector.  Default \code{1e-4}.
#' @param h0      initial head vector (length \code{n}).  Default: all zeros.
#' @param Qinf    data frame with columns \code{Qp} L\eqn{^3}/T, \code{x},
#'   \code{y} — pumping rate and location.  Positive = extraction.
#' @param times   numeric vector of output times.  **Include all observation
#'   times** so that \code{\link{samDataTr}} can extract exact values without
#'   interpolation.
#' @param lrw     real work array length passed to \code{ode.2D}.  Increase
#'   for larger grids (rule of thumb: \eqn{\geq 1000 \times n}).
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{out}}{ODE output matrix (\code{nt × (n+1)}): column 1 is
#'       time, columns 2–(n+1) are head values ordered as
#'       \code{expand.grid(x, y)}.}
#'     \item{\code{grid}}{grid list from \code{GenGrid()}.}
#'     \item{\code{times}}{the time vector actually used (= \code{out[,1]}).}
#'   }
#' @export
#' @examples
#' # Homogeneous T and S, single pumping well
#' times <- c(0, 1, 5, 10, 50, 100)
#' res   <- Ftransient2dsim(times = times)
#' # Head at the final time step, reshaped to matrix
#' grid2d <- res$grid
#' h_end  <- matrix(res$out[nrow(res$out), -1], grid2d$nx, grid2d$ny)
#'
#' # Heterogeneous T field
#' set.seed(1)
#' TT  <- random2d(nsim = 1)$Tp
#' res <- Ftransient2dsim(TT = TT, SS = 1e-3, times = seq(0, 100, by = 10))
#'
#' # ---- Compare with Theis analytical solution ----
#' \dontrun{
#' # Large domain, pumping at centre, observation well 200 m away
#' Q  <- 1.3e-3; T_val <- 1.5e-3; S_val <- 2e-5; rw <- 200
#' domain <- c(80, 80, 0, 800, 0, 800)
#' grid   <- GenGrid(domain)
#' times  <- seq(0, 1000, by = 1)
#' res2d  <- Ftransient2dsim(domain = domain, grid = grid,
#'                           TT = T_val, SS = S_val,
#'                           Qinf = data.frame(Qp = Q, x = 400, y = 400),
#'                           times = times, lrw = 2000000)
#' # Extract drawdown at observation well
#' obs_t  <- c(10, 30, 100, 300, 1000)
#' Oinf   <- data.frame(data = NA, x = 600, y = 400, time = obs_t)
#' Oinf   <- samDataTr(Oinf = Oinf, grid = grid, result_tr = res2d)
#' s_2d   <- -Oinf$data
#' # Theis analytical solution
#' s_theis <- theis_drawdown(Q, rw, T_val, S_val, obs_t)
#' # Plot comparison with ggplot2
#' library(ggplot2)
#' ggplot() +
#'   geom_line(aes(obs_t, s_theis, color = "Theis"), linewidth = 1) +
#'   geom_point(aes(obs_t, s_2d, color = "2D"), size = 2.5) +
#'   scale_x_log10() +
#'   labs(x = "t (s)", y = "s (m)", color = NULL) +
#'   theme_bw()
#' }
Ftransient2dsim <- function(domain = c(40, 40, 0, 40, 0, 40),
                             grid   = NULL,
                             TT     = 0.1,
                             SS     = 1e-4,
                             h0     = NULL,
                             Qinf   = data.frame(Qp = 10, x = 20.5, y = 20.5),
                             times  = seq(0, 100, by = 10),
                             lrw    = 1600000) {

  require('deSolve')

  if (is.null(grid)) grid <- GenGrid(domain)
  n  <- grid$n
  nx <- grid$nx; ny <- grid$ny
  dx <- grid$dx; dy <- grid$dy

  Tp <- if (length(TT) == 1L) rep(TT, n) else TT
  Sp <- if (length(SS) == 1L) rep(SS, n) else SS
  Tpm <- matrix(Tp, nx, ny)
  Spm <- matrix(Sp, nx, ny)
  if (is.null(h0)) h0 <- rep(0, n)

  Qseq <- getQseq(grid = grid, Qinf = Qinf)
  Nxp  <- Qseq$Nxp; Nyp <- Qseq$Nyp; Qp <- Qinf$Qp

  para <- list(dx = dx, dy = dy, nx = nx, ny = ny,
               Tpm = Tpm, Spm = Spm, Qp = Qp, Nxp = Nxp, Nyp = Nyp)

  out <- ode.2D(y      = h0,
                times  = times,
                func   = diffusion2D,
                parms  = para,
                dimens = c(nx, ny),
                lrw    = lrw)

  list(out   = unclass(out),   # plain matrix; drop deSolve class for easier indexing
       grid  = grid,
       times = out[, 1])
}


#' Sample transient 2D head values at observation locations and times
#'
#' Extracts simulated head values from the output of
#' \code{\link{Ftransient2dsim}} at the spatial locations and observation
#' times listed in \code{Oinf}.
#'
#' If an observation time does not exactly match a row in \code{result_tr},
#' the nearest available time step is used.  To get exact values, include all
#' observation times in the \code{times} argument of
#' \code{\link{Ftransient2dsim}}.
#'
#' @param Oinf      data frame with columns \code{data}, \code{x}, \code{y},
#'   \code{time}.
#' @param grid      grid from \code{GenGrid()}.
#' @param result_tr list returned by \code{\link{Ftransient2dsim}}.
#' @return \code{Oinf} with the \code{data} column filled.
#' @export
#' @examples
#' times  <- seq(0, 100, by = 10)
#' res    <- Ftransient2dsim(times = times)
#' Oinf   <- data.frame(data = NA,
#'                      x    = c(10.5, 25.5, 30.5),
#'                      y    = c(20.5, 15.5, 25.5),
#'                      time = c(10,   50,   100))
#' Oinf   <- samDataTr(Oinf = Oinf, grid = res$grid, result_tr = res)
samDataTr <- function(Oinf,
                      grid,
                      result_tr) {

  out_mat   <- result_tr$out        # nt × (n+1); col 1 = time
  sim_times <- out_mat[, 1]
  Oinf_elem <- getOelem(grid = grid, Oinf = Oinf)

  Oinf$data <- mapply(
    function(nelem, t_obs) {
      it <- which.min(abs(sim_times - t_obs))
      out_mat[it, nelem + 1L]       # +1 offset: col 1 is time
    },
    Oinf_elem$nelem,
    Oinf$time
  )
  Oinf
}
