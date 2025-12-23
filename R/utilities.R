#' Generate 2-D random field of hydraulic parameters
#'
#' Based on package gstat, it generate correlated random field using Gaussian simulation.

#' @param domain c(Nx,Ny,x1,x2,y1,y2) domain of the grid.
#' @param geo geostatistical features: list(mean,variance,covarirance_function,anisotropy,range,nugget)
#' @param nsim number of simulations.
#' @return A data frame of three elements with grid showing generated grid (in x and y) and
#'    parameter values.
#' @export
#' @examples
#' TT= random2d()
#' library('ggplot2')
#' ggplot(TT) +
#'   aes(x = x, y = y, fill = Tp) +
#'   geom_tile() +
#'   scale_fill_viridis_c(option = "viridis", direction = 1,trans="log") +
#'   theme_minimal()
#'  TT= random2d(geo=list(me=0,var=1,geomod="Exp",anis=c(90,0.2),range=30,nugget=0))
#' ggplot(TT) +
#'   aes(x = x, y = y, fill = Tp) +
#'   geom_tile() +
#'   scale_fill_viridis_c(option = "viridis", direction = 1,trans="log") +
#'   theme_minimal()
random2d <- function(domain=c(40,40,0,40,0,40),
                     geo=list(me=0,var=1,geomod="Exp",anis=c(90,1),range=30,nugget=0),
                     nsim=1){

  require('gstat')
  Nx <- domain[1]
  Ny <- domain[2]
  # input
  x <- seq(domain[3],domain[4],length.out=Nx)
  y <- seq(domain[5],domain[6],length.out=Ny)
  xy <- expand.grid(x=x, y=y)

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


#' generate forward running grid and output its grid and observation/source grid.
#'  input domain, pump and obs location
#' @param domain c(Nx,Ny,x1,x2,y1,y2) domain of the grid.
#' @param pump c(xp,yp) pump location.
#' @param obs c(xo,yo) observation location.
#' @return output grid, pump and obs element number
#' @export
#' @examples
#' grid <- GenGrid()
GenGrid <- function(domain=c(40,40,0,40,0,40),
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




#' sample obs for synthetic data.
#' @export
#' @param loc_obs c(elem1,elem2) elem number of the observation
samData <- function(loc_obs=c(1,2),h=rnorm(10)){
  obs <- h[loc_obs]
  return(obs)
}
