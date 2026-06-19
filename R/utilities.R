#' Generate 2-D random field of transmissivity
#'
#' Based on package \pkg{gstat}, generates a spatially correlated 2D random field
#' of log-transmissivity \eqn{\ln T} using unconditional Gaussian simulation,
#' then returns \eqn{T = \exp(\ln T)}.  The variogram model supports exponential
#' (\code{"Exp"}), Gaussian (\code{"Gau"}), spherical (\code{"Sph"}), and other
#' \code{gstat} model types.  Anisotropy is specified by the azimuth angle
#' (degrees clockwise from north) and the anisotropy ratio.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.  \code{nx}
#'   and \code{ny} are cell counts; \code{x1:x2} and \code{y1:y2} are domain
#'   extents \code{[L]}.
#' @param grid grid list from \code{GenGrid()}; generated from \code{domain} if
#'   \code{NULL}.
#' @param geo geostatistical parameters as a named list:
#'   \describe{
#'     \item{me}{mean of ln(T), default 0}
#'     \item{var}{variance of ln(T), default 1}
#'     \item{geomod}{variogram model type, e.g. \code{"Exp"}, \code{"Gau"},
#'       \code{"Sph"}; default \code{"Exp"}}
#'     \item{anis}{2D anisotropy: 2-element vector \code{c(azimuth, ratio)}.
#'       \code{azimuth} = anisotropy direction (degrees clockwise from north),
#'       \code{ratio} = range ratio in the minor direction (1 = isotropic).
#'       Default \code{c(90, 1)}.}
#'     \item{range}{correlation range \code{[L]}, default 30}
#'     \item{nugget}{nugget variance, default 0}
#'   }
#' @param nsim number of realisations to generate (default 1).
#' @return A data frame with columns \code{x}, \code{y} and \code{nsim} columns
#'   of \eqn{T} \code{[L\eqn{^2}/T]} values named \code{sim1}, \code{sim2}, … (or a
#'   single column \code{Tp} when \code{nsim = 1}).
#' @export
#' @examples
#' TT= random2d()
#' Plotparameter2d(TT,iflog=T)
#'  TT= random2d(geo=list(me=0,var=1,geomod="Exp",anis=c(90,0.2),range=30,nugget=0))
#' Plotparameter2d(TT,iflog=T)

random2d <- function(domain=c(40,40,0,40,0,40),
                     grid = NULL,
                     geo=list(me=0,var=1,geomod="Exp",anis=c(90,1),range=30,nugget=0),
                     nsim=1){

  require('gstat')
  if(is.null(grid))grid = GenGrid(domain)
  # Nx <- domain[1]
  # Ny <- domain[2]
  # # input
  # x <- seq(domain[3],domain[4],length.out=Nx)
  # y <- seq(domain[5],domain[6],length.out=Ny)
  # xy <- expand.grid(x=x, y=y)
  x = grid$xmid
  y = grid$ymid
  xy = grid$grid
  me <- geo$me #mean value. lnK
  var <- geo$var # variance.
  geomod <- geo$geomod
  range <- geo$range
  nugget <- geo$nugget
  anis <- geo$anis
  m1 <- vgm(psill=var,model=geomod,range=range,anis=anis,nugget=nugget)
  g.dummy <- gstat(formula=z~1, locations=~x+y,
                   dummy=T, beta=me,
                   model=m1, nmax=20)

  # make four simulations based on the gstat object
  yy <- predict(g.dummy, newdata=xy, nsim=nsim)
  Tp <- exp(yy[,-c(1,2)])
  return(data.frame(xy,Tp))
}

#' Generate 3-D random field of hydraulic conductivity
#'
#' Based on package gstat, generates a spatially correlated 3D random field of
#' log-hydraulic conductivity using unconditional Gaussian simulation, then
#' returns \eqn{K = \exp(\ln K)}.
#'
#' @param domain a 9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2, z1, z2)}.
#' @param grid   grid from \code{GenGrid3D()}; generated from \code{domain} if
#'   \code{NULL}.
#' @param geo    geostatistical properties as a named list:
#'   \describe{
#'     \item{me}{mean of ln(K), default 0}
#'     \item{var}{variance of ln(K), default 1}
#'     \item{geomod}{variogram model type, e.g. \code{"Exp"}, \code{"Gau"},
#'       default \code{"Exp"}}
#'     \item{range}{correlation range, default 30}
#'     \item{nugget}{nugget variance, default 0}
#'     \item{anis}{3D anisotropy: 5-element vector
#'       \code{c(alpha, beta, phi, ratio_y, ratio_z)}.
#'       \code{alpha} = azimuth (deg, rotation around z), \code{beta} = dip
#'       (deg, rotation around y), \code{phi} = plunge (deg, rotation around x),
#'       \code{ratio_y} and \code{ratio_z} are range ratios in the minor
#'       directions (1 = isotropic).  Default \code{c(0, 0, 0, 1, 1)}.}
#'   }
#' @param nsim   number of realisations to generate.
#' @return A data frame with columns \code{x}, \code{y}, \code{z} and
#'   \code{nsim} columns of K values named \code{sim1}, \code{sim2}, …
#'   (or a single unnamed column when \code{nsim = 1}).
#' @export
#' @examples
#' KK <- random3d(nsim = 1)
#' # Anisotropic: longer range in x (ratio_y = 0.3, ratio_z = 0.1)
#' KK <- random3d(geo = list(me = 0, var = 1, geomod = "Exp",
#'                           range = 30, nugget = 0,
#'                           anis = c(0, 0, 0, 0.3, 0.1)))
random3d <- function(domain = c(40, 40, 10, 0, 40, 0, 40, 0, 10),
                     grid   = NULL,
                     geo    = list(me = 0, var = 1, geomod = "Exp",
                                   range = 30, nugget = 0,
                                   anis = c(0, 0, 0, 1, 1)),
                     nsim   = 1) {

  require('gstat')
  if (is.null(grid)) grid <- GenGrid3D(domain)

  xyz    <- grid$grid
  me     <- geo$me
  var    <- geo$var
  geomod <- geo$geomod
  range  <- geo$range
  nugget <- geo$nugget
  anis   <- geo$anis

  m1 <- vgm(psill = var, model = geomod, range = range,
             anis = anis, nugget = nugget)

  # 'value~1' avoids naming conflict with the spatial coordinate column 'z'
  g.dummy <- gstat(formula   = value ~ 1,
                   locations = ~ x + y + z,
                   dummy     = TRUE,
                   beta      = me,
                   model     = m1,
                   nmax      = 20)

  yy <- predict(g.dummy, newdata = xyz, nsim = nsim)
  # predict() returns xyz columns first, then sim1, sim2, ...
  Kp <- exp(yy[, -c(1, 2, 3)])
  data.frame(xyz, Kp)
}


#' Generate a 2D finite difference grid
#'
#' Creates a uniform 2D grid for finite difference simulations.  The grid cells
#' are centred at half-integer node positions and numbered by
#' \code{expand.grid(x, y)} (x varies fastest).  The grid list can be passed to
#' all forward and inverse functions to avoid recomputing grid geometry.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.
#'   \code{nx} and \code{ny} are cell counts; \code{x1:x2} and \code{y1:y2}
#'   are domain extents \code{[L]}.
#' @return A list with elements:
#'   \describe{
#'     \item{n}{total number of cells (\code{nx * ny}).}
#'     \item{grid}{data frame of cell-centre coordinates (\code{x}, \code{y}).}
#'     \item{xmid, ymid}{vectors of cell-centre coordinates along each axis \code{[L]}.}
#'     \item{nx, ny}{number of cells in each direction.}
#'     \item{dx, dy}{cell sizes \code{[L]}.}
#'     \item{x1, x2, y1, y2}{domain bounds \code{[L]}.}
#'   }
#' @export
#' @examples
#' grid <- GenGrid()
#' str(grid)
GenGrid <- function(domain=c(40,40,0,40,0,40))
{
  nx <- domain[1]
  ny <- domain[2]
  x1 <- domain[3]
  x2 <- domain[4]
  y1 <- domain[5]
  y2 <- domain[6]
  dx <- (x2-x1)/nx
  dy <- (y2-y1)/ny
  n <- nx * ny
  xmid <- seq(0.5*dx,x2-0.5*dx,by=dx)
  ymid <- seq(0.5*dy,y2-0.5*dy,by=dy)
  grid <- expand.grid(x=xmid,y=ymid)
  return(list(n=n,grid=grid,xmid = xmid, ymid = ymid,
              nx=nx,ny=ny,dx=dx,dy=dy,x1=x1,x2=x2,y1=y1,y2=y2))
}

#' Map pumping well (x, y) coordinates to grid indices
#'
#' Maps each pumping well's physical coordinates to the nearest grid cell indices
#' \code{Nxp} (x-direction) and \code{Nyp} (y-direction).  These indices are used
#' to apply the source term in the finite difference scheme.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.  Ignored
#'   when \code{grid} is supplied.
#' @param grid grid list from \code{GenGrid()}.
#' @param Qinf a data frame with columns \code{Qp} \code{[L\eqn{^3}/T]}, \code{x} \code{[L]},
#'   \code{y} \code{[L]} giving pumping rate and location.
#' @return \code{Qinf} with additional columns \code{Nxp} (x-index) and
#'   \code{Nyp} (y-index).
#' @export
#' @examples
#' grid <- GenGrid()
#' Qseq <- getQseq(grid = grid)
#'
getQseq <- function(domain=c(40,40,0,40,0,40),
                    grid=NULL,
                    Qinf=data.frame(Qp=10,x=20.5,y=20.5)){
  Np <- length(Qinf$x)
  xmid = grid$xmid
  ymid = grid$ymid
  Nxp <- integer(length= Np)
  for(i in 1:Np){
    Nxp[i] <- which.min(abs(Qinf$x[i]-xmid))
  }
  Nyp <- integer(length= Np)
  for(i in 1:Np){
    Nyp[i] <- which.min(abs(Qinf$y[i]-ymid))
  }
  Qinf$Nxp = Nxp
  Qinf$Nyp = Nyp
  return(Qinf)
}

#' Map observation well (x, y) coordinates to flat element numbers
#'
#' Element numbering follows \code{expand.grid(x, y)}: x varies fastest.
#' For indices (ix, iy): \code{nelem = (iy - 1) * nx + ix}.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.  Ignored
#'   when \code{grid} is supplied.
#' @param grid grid list from \code{GenGrid()}.
#' @param Oinf a data frame with columns \code{data} \code{[L]}, \code{x} \code{[L]},
#'   \code{y} \code{[L]} giving observation drawdown and location.  \code{data} may
#'   be \code{NA} for synthetic sampling.
#' @return \code{Oinf} with an additional column \code{nelem} (flat element
#'   number).
#' @export
#' @examples
#' grid <- GenGrid()
#' Oinf=data.frame(data=NA,x=30.5,y=30.5)
#' Oinf <- getOelem(grid = grid, Oinf=Oinf)
#'
getOelem <- function(domain=c(40,40,0,40,0,40),
                    grid =NULL,
                    Oinf=data.frame(data=NA,x=30.5,y=30.5)){
  No <- length(Oinf$x)
  xmid = grid$xmid
  ymid = grid$ymid
  Nxo <- integer(length= No)
  for(i in 1:No){
    Nxo[i] <- which.min(abs(Oinf$x[i]-xmid))
  }
  Nyo <- integer(length= No)
  for(i in 1:No){
    Nyo[i] <- which.min(abs(Oinf$y[i]-ymid))
  }
 ### return the actual element number... numbering from lowerleft to upperright.
  nelem = (Nyo - 1) * grid$nx +Nxo
  Oinf$nelem = nelem
  return(Oinf)
}

#' Generate a 2D finite difference grid with pump and observation locations
#'
#' Convenience function that creates a grid and simultaneously maps a single
#' pumping well and a single observation well to their grid indices.  Useful
#' for quick tests and demonstrations.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.
#'   \code{nx} and \code{ny} are cell counts; \code{x1:x2} and \code{y1:y2}
#'   are domain extents \code{[L]}.
#' @param pump a list with elements \code{xp} \code{[L]} and \code{yp} \code{[L]} giving the
#'   pumping well coordinates.
#' @param obs a list with elements \code{xo} \code{[L]} and \code{yo} \code{[L]} giving the
#'   observation well coordinates.
#' @return A list with elements \code{n}, \code{grid}, \code{pump} (c(Nxp, Nyp)),
#'   and \code{obs} (c(Nxo, Nyo)).
#' @export
#' @examples
#' grid <- GenGrid2()
GenGrid2 <- function(domain=c(40,40,0,40,0,40),
                    pump=list(xp=20.5,yp=20.5),
                    obs=list(xo=30,yo=30))
{
  nx <- domain[1]
  ny <- domain[2]
  x1 <- domain[3]
  x2 <- domain[4]
  y1 <- domain[5]
  y2 <- domain[6]
  dx <- (x2-x1)/nx
  dy <- (y2-y1)/ny
  n <- nx * ny
  xmid <- seq(0.5*dx,x2-0.5*dx,by=dx)
  ymid <- seq(0.5*dy,y2-0.5*dy,by=dy)
  grid <- expand.grid(x=xmid,y=ymid)
  Nxp <- which.min(abs(pump$xp-xmid))
  Nyp <- which.min(abs(pump$yp-ymid))
  Nxo <- which.min(abs(obs$xo-xmid))
  Nyo <- which.min(abs(obs$yo-ymid))
  return(list(n=n,grid=grid,pump=c(Nxp,Nyp),obs=c(Nxo,Nyo)))
}


# ---- 3D Grid utilities -------------------------------------------------------

#' Generate a 3D finite difference grid
#'
#' @param domain a 9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2, z1, z2)}.
#'   \code{nx/ny/nz} are cell counts; \code{x1:x2}, \code{y1:y2}, \code{z1:z2}
#'   are domain extents.
#' @return A list with elements \code{n}, \code{grid}, \code{xmid}, \code{ymid},
#'   \code{zmid}, \code{nx}, \code{ny}, \code{nz}, \code{dx}, \code{dy},
#'   \code{dz}, and domain bounds.
#' @export
#' @examples
#' grid3d <- GenGrid3D()
#' head(grid3d$grid)
GenGrid3D <- function(domain = c(40, 40, 10, 0, 40, 0, 40, 0, 10)) {
  nx <- domain[1]; ny <- domain[2]; nz <- domain[3]
  x1 <- domain[4]; x2 <- domain[5]
  y1 <- domain[6]; y2 <- domain[7]
  z1 <- domain[8]; z2 <- domain[9]

  dx <- (x2 - x1) / nx
  dy <- (y2 - y1) / ny
  dz <- (z2 - z1) / nz
  n  <- nx * ny * nz

  xmid <- seq(x1 + 0.5 * dx, x2 - 0.5 * dx, by = dx)
  ymid <- seq(y1 + 0.5 * dy, y2 - 0.5 * dy, by = dy)
  zmid <- seq(z1 + 0.5 * dz, z2 - 0.5 * dz, by = dz)

  # expand.grid varies x fastest, then y, then z — matches array(, dim=c(nx,ny,nz))
  grid <- expand.grid(x = xmid, y = ymid, z = zmid)

  list(n = n, grid = grid,
       xmid = xmid, ymid = ymid, zmid = zmid,
       nx = nx, ny = ny, nz = nz,
       dx = dx, dy = dy, dz = dz,
       x1 = x1, x2 = x2, y1 = y1, y2 = y2, z1 = z1, z2 = z2)
}


#' Map pumping well (x, y, z) coordinates to grid indices (3D)
#'
#' @param grid grid generated by \code{GenGrid3D}.
#' @param Qinf a data frame with columns \code{Qp}, \code{x}, \code{y}, \code{z}.
#' @return \code{Qinf} with additional columns \code{Nxp}, \code{Nyp}, \code{Nzp}.
#' @export
#' @examples
#' grid3d <- GenGrid3D()
#' Qinf <- data.frame(Qp = 10, x = 20.5, y = 20.5, z = 5.5)
#' getQseq3D(grid = grid3d, Qinf = Qinf)
getQseq3D <- function(grid,
                      Qinf = data.frame(Qp = 10, x = 20.5, y = 20.5, z = 5.5)) {
  Np  <- nrow(Qinf)
  Nxp <- integer(Np); Nyp <- integer(Np); Nzp <- integer(Np)
  for (i in seq_len(Np)) {
    Nxp[i] <- which.min(abs(Qinf$x[i] - grid$xmid))
    Nyp[i] <- which.min(abs(Qinf$y[i] - grid$ymid))
    Nzp[i] <- which.min(abs(Qinf$z[i] - grid$zmid))
  }
  Qinf$Nxp <- Nxp; Qinf$Nyp <- Nyp; Qinf$Nzp <- Nzp
  Qinf
}


#' Map observation well (x, y, z) coordinates to flat element numbers (3D)
#'
#' Element numbering follows \code{expand.grid(x, y, z)}: x varies fastest.
#' For indices (ix, iy, iz): \code{nelem = (iz-1)*nx*ny + (iy-1)*nx + ix}.
#'
#' @param grid grid generated by \code{GenGrid3D}.
#' @param Oinf a data frame with columns \code{data}, \code{x}, \code{y}, \code{z}.
#' @return \code{Oinf} with an additional column \code{nelem}.
#' @export
#' @examples
#' grid3d <- GenGrid3D()
#' Oinf <- data.frame(data = NA, x = 30.5, y = 30.5, z = 5.5)
#' getOelem3D(grid = grid3d, Oinf = Oinf)
getOelem3D <- function(grid,
                       Oinf = data.frame(data = NA, x = 30.5, y = 30.5, z = 5.5)) {
  No  <- nrow(Oinf)
  Nxo <- integer(No); Nyo <- integer(No); Nzo <- integer(No)
  for (i in seq_len(No)) {
    Nxo[i] <- which.min(abs(Oinf$x[i] - grid$xmid))
    Nyo[i] <- which.min(abs(Oinf$y[i] - grid$ymid))
    Nzo[i] <- which.min(abs(Oinf$z[i] - grid$zmid))
  }
  Oinf$nelem <- (Nzo - 1L) * grid$nx * grid$ny +
                (Nyo - 1L) * grid$nx +
                Nxo
  Oinf
}


#' Sample simulated head values at observation well locations (3D)
#'
#' @param Oinf a data frame with columns \code{data}, \code{x}, \code{y}, \code{z}.
#' @param grid grid generated by \code{GenGrid3D}.
#' @param h the head vector from \code{Fsteady3dsim()$solution} (length \code{n}).
#' @return \code{Oinf} with the \code{data} column filled.
#' @export
#' @examples
#' grid3d <- GenGrid3D()
#' Oinf <- data.frame(data = NA, x = 10.5, y = 10.5, z = 5.5)
#' h <- rnorm(grid3d$n)
#' samData3D(Oinf = Oinf, grid = grid3d, h = h)
samData3D <- function(Oinf = data.frame(data = NA, x = 10.5, y = 10.5, z = 5.5),
                      grid,
                      h) {
  Oinf       <- getOelem3D(grid = grid, Oinf = Oinf)
  Oinf$data  <- h[Oinf$nelem]
  Oinf
}


#' Sample simulated head values at observation well locations (2D)
#'
#' Extracts simulated head values from a forward simulation result at the
#' spatial locations listed in \code{Oinf}.  Used to generate synthetic
#' observation data for inversion testing.
#'
#' @param Oinf a data frame with columns \code{data} \code{[L]}, \code{x} \code{[L]},
#'   \code{y} \code{[L]}.  The \code{data} column is filled with sampled head values.
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.  Ignored
#'   when \code{grid} is supplied.
#' @param grid grid list from \code{GenGrid()}.
#' @param h the head vector \code{[L]} from \code{Fsteady2dsim()$solution} (length
#'   \code{n}).
#' @return \code{Oinf} with the \code{data} column filled.
#' @export
#' @examples
#' grid <- GenGrid()
#' Oinf=data.frame(data=NA,x=10.5,y=10.5)
#' h = rnorm(1600)
#' Oinf <- samData(Oinf=Oinf,grid=grid,h=h)

samData <- function(Oinf=data.frame(data=NA,x=10.5,y=10.5),
                    domain=c(40,40,0,40,0,40),
                    grid=NULL,
                    h=rnorm(1600)){
  if(is.null(grid)) grid = GenGrid(domain = domain)
  Oinf <- getOelem(grid = grid,Oinf = Oinf)
  nelem = Oinf$nelem
  obs <- h[nelem]
  Oinf$data = obs
  return(Oinf)
}

#' Calculate goodness-of-fit statistics between observed and simulated values
#'
#' Computes four commonly used metrics:
#' \describe{
#'   \item{L1}{Mean absolute error: \eqn{\text{L1} = \frac{1}{N} \sum |y_i - \hat{y}_i|}}
#'   \item{RMSE}{Root mean square error: \eqn{\text{RMSE} = \sqrt{\frac{1}{N} \sum (y_i - \hat{y}_i)^2}}}
#'   \item{r}{Pearson correlation coefficient between true and simulated values.}
#'   \item{NSE}{Nash–Sutcliffe efficiency: \eqn{\text{NSE} = 1 - \frac{\sum (y_i - \hat{y}_i)^2}{\sum (y_i - \bar{y})^2}}.
#'     NSE = 1 indicates perfect fit; NSE \eqn{\leq} 0 indicates the mean is a
#'     better predictor.}
#' }
#'
#' @param df a data frame with two columns: the first column holds the true
#'   (observed) values, the second column holds the simulated values.
#' @return A data frame with columns \code{L1}, \code{RMSE}, \code{r}, and
#'   \code{NSE}.
#' @export

statData <- function(df){
  true = df[,1]
  simu = df[,2]
  rmse = sqrt(mean((true - simu)^2))
  r = cor(true, simu)
  l1 <- mean(abs(true - simu))
  y_ = mean(true)
  r2 =1 - sum((true - simu)^2)/sum((true - y_)^2)
  # https://towardsdatascience.com/explaining-negative-r-squared-17894ca26321/
  # we can treat it as nash coefficient, NSE.
  return(data.frame(L1=l1,RMSE=rmse,r=r,NSE=r2))
}
