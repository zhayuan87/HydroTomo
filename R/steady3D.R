# =============================================================================
# 3D Steady-State Groundwater Flow Simulation
# Mirrors the structure of steady2D.R, extended to three dimensions.
#
# Key differences from 2D:
#   - domain has 9 elements: c(nx, ny, nz, x1, x2, y1, y2, z1, z2)
#   - parameter field is hydraulic conductivity K \code{[L/T]} (not transmissivity T)
#   - pump/obs locations require (x, y, z) coordinates
#   - PDE solved with ReacTran::tran.3D() + rootSolve::steady.3D()
#   - element numbering follows expand.grid(x, y, z): x varies fastest
# =============================================================================

#' Solve 3D steady-state groundwater flow equation
#'
#' Solves \eqn{\nabla \cdot (K \nabla h) = Q} in a 3D confined aquifer domain
#' using finite differences (\code{ReacTran::tran.3D}) and a steady-state solver
#' (\code{rootSolve::steady.3D}).  Boundary conditions are zero drawdown on all
#' six faces.
#'
#' @param domain a 9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2, z1, z2)}.
#' @param grid   grid from \code{GenGrid3D()}; generated from \code{domain} if
#'   \code{NULL}.
#' @param KK     hydraulic conductivity \eqn{K} \code{[L/T]} — either a scalar or a
#'   length-\code{n} vector (one value per cell, ordered as
#'   \code{expand.grid(x, y, z)}).
#' @param Qinf   a data frame with columns \code{Qp} L\eqn{^3}/T, \code{x},
#'   \code{y}, \code{z} giving pumping rate and location of each pumping well.
#'   Positive \code{Qp} = extraction.
#' @param lrw    length of the real work array passed to the solver.  Increase
#'   for larger grids (rule of thumb: \eqn{\geq 100 \times n}).
#' @return A data frame with columns \code{x}, \code{y}, \code{z}, and
#'   \code{solution} (steady-state head/drawdown at every cell centre).
#' @export
#' @examples
#' # Homogeneous K, single pumping well at domain centre
#' s3 <- Fsteady3dsim()
#'
#' # Heterogeneous K field
#' grid3d <- GenGrid3D()
#' # KK <- random3d(nsim = 1, grid = grid3d)$Kp   # once random3d() is available
#' Qinf  <- data.frame(Qp = 10, x = 20.5, y = 20.5, z = 5.5)
#' s3    <- Fsteady3dsim(grid = grid3d, Qinf = Qinf)
Fsteady3dsim <- function(domain = c(40, 40, 10, 0, 40, 0, 40, 0, 10),
                         grid   = NULL,
                         KK     = 0.1,
                         Qinf   = data.frame(Qp = 10, x = 20.5, y = 20.5, z = 5.5),
                         lrw    = 20000000) {

  require('rootSolve')
  require('ReacTran')

  if (is.null(grid)) grid <- GenGrid3D(domain)

  n  <- grid$n
  nx <- grid$nx; ny <- grid$ny; nz <- grid$nz
  dx <- grid$dx; dy <- grid$dy; dz <- grid$dz

  Kp <- if (length(KK) == 1L) rep(KK, n) else KK
  h0 <- rep(0, n)   # initial guess (all zero drawdown)

  Qseq <- getQseq3D(grid = grid, Qinf = Qinf)
  Nxp  <- Qseq$Nxp; Nyp <- Qseq$Nyp; Nzp <- Qseq$Nzp
  Qp   <- Qinf$Qp

  para <- list(dx = dx, dy = dy, dz = dz,
               nx = nx, ny = ny, nz = nz,
               Kp = Kp, Qp = Qp,
               Nxp = Nxp, Nyp = Nyp, Nzp = Nzp)

  s_steady <- steady.3D(y      = h0,
                        parms  = para,
                        func   = diffusion3D_GW,
                        dimens = c(nx, ny, nz),
                        lrw    = lrw)$y

  data.frame(grid$grid, solution = s_steady)
}


#' 3D groundwater diffusion residual function
#'
#' Computes the finite-difference residual of the steady-state groundwater flow
#' equation in three dimensions.  Called internally by \code{Fsteady3dsim()} via
#' \code{rootSolve::steady.3D()}.
#'
#' @param t   time argument (unused; required by \code{steady.3D} interface).
#' @param h   current head vector (length \code{nx * ny * nz}).
#' @param par named list of parameters produced inside \code{Fsteady3dsim()}.
#' @return A list containing the residual vector (length \code{n}).
#' @export
diffusion3D_GW <- function(t, h, par) {
  with(par, {

    Km <- array(Kp, dim = c(nx, ny, nz))
    y  <- array(h,  dim = c(nx, ny, nz))

    # Zero-drawdown boundary conditions on all six faces
    BNDx <- matrix(0, nrow = ny, ncol = nz)   # x = x1 and x = x2 faces
    BNDy <- matrix(0, nrow = nx, ncol = nz)   # y = y1 and y = y2 faces
    BNDz <- matrix(0, nrow = nx, ncol = ny)   # z = z1 and z = z2 faces

    # Arithmetic-mean K at cell interfaces (same convention as diffusion2D)
    # x-direction: dim (nx+1, ny, nz)
    Kx <- array(0, dim = c(nx + 1L, ny, nz))
    Kx[1L, , ]        <- Km[1L, , ]
    Kx[nx + 1L, , ]   <- Km[nx, , ]
    if (nx > 1L) Kx[2L:nx, , ] <- (Km[1L:(nx - 1L), , ] + Km[2L:nx, , ]) * 0.5

    # y-direction: dim (nx, ny+1, nz)
    Ky <- array(0, dim = c(nx, ny + 1L, nz))
    Ky[, 1L, ]        <- Km[, 1L, ]
    Ky[, ny + 1L, ]   <- Km[, ny, ]
    if (ny > 1L) Ky[, 2L:ny, ] <- (Km[, 1L:(ny - 1L), ] + Km[, 2L:ny, ]) * 0.5

    # z-direction: dim (nx, ny, nz+1)
    Kz <- array(0, dim = c(nx, ny, nz + 1L))
    Kz[, , 1L]        <- Km[, , 1L]
    Kz[, , nz + 1L]   <- Km[, , nz]
    if (nz > 1L) Kz[, , 2L:nz] <- (Km[, , 1L:(nz - 1L)] + Km[, , 2L:nz]) * 0.5

    # Flux divergence: ∇·(K ∇h)
    dY <- tran.3D(C        = y,
                  C.x.up   = BNDx, C.x.down = BNDx,
                  C.y.up   = BNDy, C.y.down = BNDy,
                  C.z.up   = BNDz, C.z.down = BNDz,
                  D.x      = Kx, D.y = Ky, D.z = Kz,
                  dx       = dx, dy = dy, dz = dz)$dC

    # Pumping source terms: Qp [L³/T] distributed over one cell volume
    for (p in seq_along(Qp)) {
      dY[Nxp[p], Nyp[p], Nzp[p]] <-
        dY[Nxp[p], Nyp[p], Nzp[p]] - Qp[p] / dx / dy / dz
    }

    list(as.vector(dY))
  })
}
