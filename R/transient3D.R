# =============================================================================
# 3D Transient Groundwater Flow Simulation
#
# Solves  Ss ∂h/∂t = ∇·(K ∇h) + Q  using deSolve::ode.3D().
#
# Key relationship to steady-state code:
#   - diffusion3D_GW() (steady3D.R) does NOT divide by Ss, so a new function
#     diffusion3D_tr() is introduced that adds the Ss storage term.
#   - Ftransient3dsim() calls ode.3D() with diffusion3D_tr() as the ODE
#     right-hand side, exposing Ss, h0, and times as new parameters.
#   - samData3DTr() samples the ode.3D output at (x, y, z, time) triplets.
# =============================================================================


#' 3D transient groundwater flow ODE residual
#'
#' Computes \eqn{\partial h / \partial t = S_s^{-1} \left[\nabla \cdot
#' (K \nabla h) - Q_\text{src}\right]} for use with \code{deSolve::ode.3D()}.
#' Extends \code{\link{diffusion3D_GW}} by adding the specific storage term.
#'
#' @param t   time (passed by the ODE solver; may be used for time-varying Q
#'   in future extensions).
#' @param h   current head vector (length \code{nx * ny * nz}).
#' @param par named list of parameters built inside
#'   \code{\link{Ftransient3dsim}}: \code{dx, dy, dz, nx, ny, nz, Kp, Ssp,
#'   Qp, Nxp, Nyp, Nzp}.
#' @return A list containing the \eqn{\partial h / \partial t} vector.
#' @export
diffusion3D_tr <- function(t, h, par) {
  with(par, {

    Km  <- array(Kp,  dim = c(nx, ny, nz))
    Ssm <- array(Ssp, dim = c(nx, ny, nz))
    y   <- array(h,   dim = c(nx, ny, nz))

    # ---- face-centred K values (arithmetic mean, same as diffusion3D_GW) ----
    Kx <- array(0, dim = c(nx + 1L, ny, nz))
    Kx[1L, , ]      <- Km[1L, , ]
    Kx[nx + 1L, , ] <- Km[nx, , ]
    if (nx > 1L) Kx[2L:nx, , ] <- (Km[1L:(nx-1L), , ] + Km[2L:nx, , ]) * 0.5

    Ky <- array(0, dim = c(nx, ny + 1L, nz))
    Ky[, 1L, ]      <- Km[, 1L, ]
    Ky[, ny + 1L, ] <- Km[, ny, ]
    if (ny > 1L) Ky[, 2L:ny, ] <- (Km[, 1L:(ny-1L), ] + Km[, 2L:ny, ]) * 0.5

    Kz <- array(0, dim = c(nx, ny, nz + 1L))
    Kz[, , 1L]      <- Km[, , 1L]
    Kz[, , nz + 1L] <- Km[, , nz]
    if (nz > 1L) Kz[, , 2L:nz] <- (Km[, , 1L:(nz-1L)] + Km[, , 2L:nz]) * 0.5

    # ---- zero-drawdown boundary conditions on all six faces -----------------
    BNDx <- matrix(0, nrow = ny, ncol = nz)
    BNDy <- matrix(0, nrow = nx, ncol = nz)
    BNDz <- matrix(0, nrow = nx, ncol = ny)

    # ---- flux divergence  ∇·(K ∇h) -----------------------------------------
    dY <- ReacTran::tran.3D(
      C        = y,
      C.x.up   = BNDx, C.x.down = BNDx,
      C.y.up   = BNDy, C.y.down = BNDy,
      C.z.up   = BNDz, C.z.down = BNDz,
      D.x      = Kx, D.y = Ky, D.z = Kz,
      dx       = dx, dy = dy, dz = dz
    )$dC

    # ---- pumping source terms  Q [L³/T] / cell volume -----------------------
    for (p in seq_along(Qp)) {
      if (is.list(Nzp)) {
        # well-screen pumping: distribute Qp evenly across z layers
        nzp_vec <- Nzp[[p]]
        nlay    <- length(nzp_vec)
        for (iz in nzp_vec) {
          dY[Nxp[p], Nyp[p], iz] <-
            dY[Nxp[p], Nyp[p], iz] - Qp[p] / dx / dy / dz / nlay
        }
      } else {
        dY[Nxp[p], Nyp[p], Nzp[p]] <-
          dY[Nxp[p], Nyp[p], Nzp[p]] - Qp[p] / dx / dy / dz
      }
    }

    # ---- divide by specific storage to obtain ∂h/∂t -------------------------
    dY <- dY / Ssm

    list(as.vector(dY))
  })
}


#' Solve 3D transient groundwater flow equation
#'
#' Solves \eqn{S_s \partial h / \partial t = \nabla \cdot (K \nabla h) + Q}
#' in a 3D confined aquifer using \code{deSolve::ode.3D()}.  Zero-drawdown
#' boundary conditions on all six faces.
#'
#' @param domain  9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2,
#'   z1, z2)}.
#' @param grid    grid from \code{GenGrid3D()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param KK      hydraulic conductivity \eqn{K} L/T — scalar or length-n
#'   vector (one value per cell, ordered as \code{expand.grid(x, y, z)}).
#' @param Ss      specific storage \eqn{S_s} 1/L — scalar or length-n
#'   vector.  Default \code{1e-4}.
#' @param h0      initial head vector (length \code{n}).  Default: all zeros.
#' @param Qinf    data frame with columns \code{Qp} L\eqn{^3}/T, \code{x},
#'   \code{y}, \code{z} (point pumping) or \code{Qp}, \code{x}, \code{y},
#'   \code{z_top}, \code{z_bottom} (well-screen pumping).
#' @param times   numeric vector of output times.  **Include all observation
#'   times** to allow exact sampling via \code{\link{samData3DTr}}.
#' @param lrw     real work array length.  Default 20000000; increase for
#'   larger grids.
#' @return A list:
#'   \describe{
#'     \item{\code{out}}{ODE output matrix (\code{nt × (n+1)}): column 1 is
#'       time, columns 2–(n+1) are head values ordered as
#'       \code{expand.grid(x, y, z)}.}
#'     \item{\code{grid}}{grid list from \code{GenGrid3D()}.}
#'     \item{\code{times}}{time vector (\code{= out[,1]}).}
#'   }
#' @export
#' @examples
#' times  <- c(0, 1, 10, 50, 100)
#' res3d  <- Ftransient3dsim(times = times)
#' # Head at the last time step as a 3D array
#' g      <- res3d$grid
#' h_end  <- array(res3d$out[nrow(res3d$out), -1], dim = c(g$nx, g$ny, g$nz))
#'
#' # Heterogeneous K field (small domain for speed)
#' grid3d <- GenGrid3D(c(15, 15, 5, 0, 15, 0, 15, 0, 5))
#' set.seed(7)
#' KK     <- random3d(nsim = 1, grid = grid3d)$Kp
#' res3d  <- Ftransient3dsim(grid = grid3d, KK = KK, Ss = 1e-4,
#'                           times = c(0, 10, 50, 100))
Ftransient3dsim <- function(domain = c(40, 40, 10, 0, 40, 0, 40, 0, 10),
                             grid   = NULL,
                             KK     = 0.1,
                             Ss     = 1e-4,
                             h0     = NULL,
                             Qinf   = data.frame(Qp = 10, x = 20.5, y = 20.5, z = 5.5),
                             times  = seq(0, 100, by = 10),
                             lrw    = 20000000) {

  require('deSolve')
  require('ReacTran')

  if (is.null(grid)) grid <- GenGrid3D(domain)
  n  <- grid$n
  nx <- grid$nx; ny <- grid$ny; nz <- grid$nz
  dx <- grid$dx; dy <- grid$dy; dz <- grid$dz

  Kp  <- if (length(KK) == 1L) rep(KK, n) else KK
  Ssp <- if (length(Ss) == 1L) rep(Ss, n) else Ss
  if (is.null(h0)) h0 <- rep(0, n)

  # ---- detect well-screen pumping (z_top / z_bottom) vs point pumping (z) ----
  if (all(c('z_top', 'z_bottom') %in% names(Qinf))) {
    # well-screen pumping: map to z-index vectors
    Nxp <- integer(nrow(Qinf)); Nyp <- integer(nrow(Qinf))
    Nzp <- vector('list', nrow(Qinf))
    for (i in seq_len(nrow(Qinf))) {
      Nxp[i] <- which.min(abs(Qinf$x[i] - grid$xmid))
      Nyp[i] <- which.min(abs(Qinf$y[i] - grid$ymid))
      z_lo   <- min(Qinf$z_top[i], Qinf$z_bottom[i])
      z_hi   <- max(Qinf$z_top[i], Qinf$z_bottom[i])
      iz_vec <- which(grid$zmid >= z_lo & grid$zmid <= z_hi)
      if (length(iz_vec) == 0) {
        iz_vec <- which.min(abs(grid$zmid - mean(c(z_lo, z_hi))))
      }
      Nzp[[i]] <- iz_vec
    }
  } else {
    Qseq <- getQseq3D(grid = grid, Qinf = Qinf)
    Nxp  <- Qseq$Nxp; Nyp <- Qseq$Nyp; Nzp <- Qseq$Nzp
  }
  Qp <- Qinf$Qp

  para <- list(dx = dx, dy = dy, dz = dz,
               nx = nx, ny = ny, nz = nz,
               Kp = Kp, Ssp = Ssp,
               Qp = Qp, Nxp = Nxp, Nyp = Nyp, Nzp = Nzp)

  out <- ode.3D(y      = h0,
                times  = times,
                func   = diffusion3D_tr,
                parms  = para,
                dimens = c(nx, ny, nz),
                lrw    = lrw)

  list(out   = unclass(out),
       grid  = grid,
       times = out[, 1])
}


#' Sample transient 3D head values at observation locations and times
#'
#' Extracts simulated head values from the output of
#' \code{\link{Ftransient3dsim}} at the spatial locations and observation
#' times listed in \code{Oinf}.  The nearest available time step is used when
#' an observation time does not exactly match the ODE output.
#'
#' @param Oinf      data frame with columns \code{data}, \code{x}, \code{y},
#'   \code{z}, \code{time}.
#' @param grid      grid from \code{GenGrid3D()}.
#' @param result_tr list returned by \code{\link{Ftransient3dsim}}.
#' @return \code{Oinf} with the \code{data} column filled.
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15, 15, 5, 0, 15, 0, 15, 0, 5))
#' times  <- c(0, 10, 50, 100)
#' res3d  <- Ftransient3dsim(grid = grid3d, times = times)
#' Oinf   <- data.frame(data = NA,
#'                      x    = c(5.5, 10.5),
#'                      y    = c(7.5,  7.5),
#'                      z    = c(2.5,  2.5),
#'                      time = c(50,  100))
#' Oinf   <- samData3DTr(Oinf = Oinf, grid = grid3d, result_tr = res3d)
samData3DTr <- function(Oinf,
                        grid,
                        result_tr) {

  out_mat   <- result_tr$out      # nt × (n+1); col 1 = time
  sim_times <- out_mat[, 1]
  Oinf_elem <- getOelem3D(grid = grid, Oinf = Oinf)

  Oinf$data <- mapply(
    function(nelem, t_obs) {
      it <- which.min(abs(sim_times - t_obs))
      out_mat[it, nelem + 1L]     # +1 offset: col 1 is time
    },
    Oinf_elem$nelem,
    Oinf$time
  )
  Oinf
}


#' Sample transient 3D head values at well-screen intervals
#'
#' Transient counterpart of \code{\link{samData3DScreen}}.  For each
#' observation well with a vertical screen defined by \code{z_top} and
#' \code{z_bottom}, computes the arithmetic mean of simulated heads at the
#' requested observation time over all grid cells within the screen interval.
#'
#' @param Oinf      data frame with columns \code{data}, \code{x}, \code{y},
#'   \code{z_top}, \code{z_bottom}, \code{time}.
#' @param grid      grid from \code{\link{GenGrid3D}}.
#' @param result_tr list returned by \code{\link{Ftransient3dsim}}.
#' @return \code{Oinf} with the \code{data} column filled with
#'   interval-averaged heads.
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15, 15, 5, 0, 15, 0, 15, 0, 5))
#' times  <- c(0, 10, 50, 100)
#' res3d  <- Ftransient3dsim(grid = grid3d, times = times)
#' Oinf   <- data.frame(data = NA,
#'                      x      = c(5.5, 10.5),
#'                      y      = c(7.5,  7.5),
#'                      z_top  = c(1.0,  1.0),
#'                      z_bottom = c(4.0, 4.0),
#'                      time   = c(50, 100))
#' Oinf   <- samData3DTrScreen(Oinf = Oinf, grid = grid3d, result_tr = res3d)
samData3DTrScreen <- function(Oinf,
                              grid,
                              result_tr) {

  out_mat   <- result_tr$out      # nt × (n+1); col 1 = time
  sim_times <- out_mat[, 1]
  Oinf_elem <- getOelem3DScreen(grid = grid, Oinf = Oinf)

  Oinf$data <- mapply(
    function(nelem_vec, t_obs) {
      it <- which.min(abs(sim_times - t_obs))
      mean(out_mat[it, nelem_vec + 1L])
    },
    Oinf_elem$nelem_list,
    Oinf$time
  )
  Oinf
}
