#' Solve 2-D steady-state groundwater flow equation
#'
#' Solve 2-D steady-state groundwater flow equation. The boundary conditions at
#' both sides, the source term, and the parameters should be specified. A grid to
#' run the simulation should be given and parameters and source term rely on the
#' definition of grid.
#'
#' @return A data frame of two elements with grid showing generated grid and
#'    steady-state solution of head distribution.
#' @param domain a 6 * 1 vector giving (nx,ny,x1,x2,y1,y2)
#' @param TT a 1 * n vector giving transimisivity (L2/T) at each grid point.
#' @param Qinf a list or data.frame giving pumping rates (Qp) and locations (xp,yp).
#' @export
#' @examples
#' s <- Fsteady2dsim()
#' library('ggplot2')
#' ggplot(s) +
#'   aes(x = x, y = y, fill = solution) +
#'   geom_tile() +
#'   scale_fill_viridis_c(option = "viridis", direction = 1) +
#'   theme_minimal()
Fsteady2dsim <- function(domain=c(40,40,0,40,0,40),
                         TT=0.1,
                         Qinf=list(Qp=10,xp=20.5,yp=20.5)
){

  require('rootSolve')
  nx <- domain[1]
  ny <- domain[2]
  x1 <- domain[3]
  x2 <- domain[4]
  y1 <- domain[5]
  y2 <- domain[6]
  dx <- (x2-x1)/nx
  dy <- (y2-y1)/ny
  n <- nx * ny
  Tp <- if(length(TT)==1) rep(TT,n) else TT
  Sp <- rep(1,n) # useless in steady-state solution.
  h <- rep(0,n)
  h0 <- matrix(0,nx,ny)
  Qp <- Qinf$Qp
  xmid <- seq(0.5*dx,x2-0.5*dx,by=dx)
  ymid <- seq(0.5*dy,y2-0.5*dy,by=dy)
  Nxp <- which.min(abs(Qinf$xp-xmid))
  Nyp <- which.min(abs(Qinf$yp-ymid))
  # solve the equation
  y <- rep(0,n)
  para =list(dx=dx,dy=dy,nx=nx,ny=ny,Tp=Tp,Sp=Sp,Qp=Qp,Nxp=Nxp,Nyp=Nyp)
  s_steady <- steady.2D(y = y, parms=para,func = diffusion2D, dimens = c(nx,ny),lrw = 160000)$y
  xy <- expand.grid(x=xmid, y = ymid)
  return(data.frame(xy, solution = s_steady))
}

#' 2-D diffusion function solver
#' @param t time
#' @param h initial condition (not important in steady-state solution)
#' @param par parameters
#' @return list of state variables (e.g., head values)
#' @export

diffusion2D <- function(t, h, par)   {
  with(par,{
    Tpm <- matrix(Tp,nx,ny) # impose at the interface.
    Spm <- matrix(Sp,nx,ny)
    #Dx <- (Tpm[-nx,] + Tpm[-1,]) * 0.5
    #Dy <- (Tpm[,-ny] + Tpm[,-1]) * 0.5
    Dx <- (rbind(Tpm[1,],Tpm) + rbind(Tpm,Tpm[nx,])) * 0.5
    Dy <- (cbind(Tpm[,1],Tpm) + cbind(Tpm,Tpm[,ny])) * 0.5
    y    <- matrix(nr=nx,nc=ny,data=h)  # vector to 2-D matrix
    dY   <- matrix(nr=nx,nc=ny,data=0)  # initial condition
    BNDx1 <- BNDx2 <- rep(0,nx)   # boundary drawdown == 0
    BNDy1 <- BNDy2 <- rep(0,ny)   # boundary drawdown == 0
    #diffusion in X-direction; boundaries=imposed concentration: nx+1, ny
    Flux <- -Dx * rbind(y[1,]-BNDx1,(y[2:nx,]-y[1:(nx-1),]),BNDx2-y[nx,])/dx
    dY   <- dY - (Flux[2:(nx+1),]-Flux[1:nx,])/dx

    #diffusion in Y-direction
    Flux <- -Dy * cbind(y[,1]-BNDy1,(y[,2:ny]-y[,1:(ny-1)]),BNDy2-y[,ny])/dy
    dY    <- dY - (Flux[,2:(ny+1)]-Flux[,1:ny])/dy
    dY <- dY / Spm
    # add pumping rate locations.
    dY[Nxp,Nyp] <- dY[Nxp,Nyp] - Qp/dx/dy/Spm[Nxp,Nyp]

    return(list(as.vector(dY)))
  })

}
