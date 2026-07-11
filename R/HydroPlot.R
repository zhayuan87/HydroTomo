#' Plot 2D steady-state head or drawdown field
#'
#' Creates a filled contour plot (using \code{ggplot2}) of the steady-state head
#' or drawdown field from \code{\link{Fsteady2dsim}}.  Pumping well and
#' observation well locations can be overlaid.
#'
#' @param s a data frame returned by \code{Fsteady2dsim()} (columns \code{x},
#'   \code{y}, \code{solution}).
#' @param plotshow a logical value; if \code{TRUE} (default), the plot is
#'   displayed.
#' @param plotfile optional file path; if supplied the plot is saved as a PDF.
#' @param plotwidth,plotheight plot dimensions in inches.
#' @param ifdrawdown logical; if \code{TRUE} (default), negate the solution to
#'   display drawdown (positive = water-level decline).
#' @param palette a string, the \code{viridis} palette option to use:
#'   \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"viridis"}
#'   (default), \code{"cividis"}, \code{"rocket"}, \code{"mako"}, \code{"turbo"}.
#' @param Qinf data frame with columns \code{Qp}, \code{x}, \code{y} for pumping
#'   wells (plotted as red circles).  \code{NULL} to omit.
#' @param Oinf data frame with columns \code{data}, \code{x}, \code{y} for
#'   observation wells (plotted as yellow squares).  \code{NULL} to omit.
#' @param title plot title.  Automatically set if \code{NULL}.
#' @param z colour-bar label.  Automatically set if \code{NULL}.
#' @param xlab,ylab axis labels (default \code{"X/m"}, \code{"Y/m"}).
#' @return A \code{ggplot} object (invisibly when \code{plotshow = TRUE}).
#' @export
#' @examples
#' Qinf = data.frame(Qp=c(10,10),x=c(10.5,30.5),y=c(20.5,20.5))
#' s <- Fsteady2dsim(Qinf=Qinf)
#' Plotsteady2d(s,Qinf=Qinf,Oinf=data.frame(data=NA,x=c(15.5,25.5),y=c(20.5,20.5)))

Plotsteady2d <- function(s, palette="viridis",plotshow = TRUE,
                         plotfile = NULL, plotwidth = 10, plotheight = 10, ifdrawdown = TRUE,
                         Qinf = NULL, Oinf = NULL,
                         title =NULL, z = NULL,xlab="X/m",ylab="Y/m") {

  if (ifdrawdown) {
    s$solution = - s$solution
    if(!is.null(title)) title = "2D Steady-State Drawdown Distribution"
    if(!is.null(z))z = "Drawdown"
  } else {
    if(!is.null(title))title = "2D Steady-State Head Distribution"
    if(!is.null(z))z = "Head"
  }

  require('ggplot2')
  p = ggplot() +
    geom_tile(data = s, aes(x, y, fill = solution)) +
    scale_fill_viridis_c(option = palette)+
    labs(title = title, x = xlab, y = ylab, fill = z)
# add the location of pumping well and observation well if Qinf and Oinf are given.
  if(is.null(Qinf)==FALSE){
    p = p + geom_point(data=Qinf,aes(x=x,y=y),color='red',fill="white",shape=21,size=4)+
      geom_text(data=Qinf,aes(x=x,y=y,label=paste0("Q=",Qp)),color='red',vjust=-1)
  }
  if(is.null(Oinf)==FALSE){
    p = p + geom_point(data=Oinf,aes(x=x,y=y),color='yellow',fill="white",shape=22,size=3)
  }


  if(!is.null(plotfile)){
    ggsave(plotfile, width = plotwidth, height = plotheight)
  }
  if(plotshow) return(p)
}

#' Plot 2D transmissivity field
#'
#' Creates a filled contour plot (using \code{ggplot2}) of the transmissivity
#' field from \code{\link{random2d}}.  Supports both raw \eqn{T} and
#' log-transformed \eqn{\ln T} display.
#'
#' @param TT a data frame returned by \code{random2d()} (columns \code{x},
#'   \code{y}, \code{Tp} or \code{sim1}, …).
#' @param plotshow a logical value; if \code{TRUE} (default), the plot is
#'   displayed.
#' @param plotfile optional file path for saving as PDF.
#' @param plotwidth,plotheight plot dimensions in inches.
#' @param iflog logical; if \code{TRUE}, plot \eqn{\ln T} instead of \eqn{T}.
#' @param palette a string, the \code{viridis} palette option to use:
#'   \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"viridis"}
#'   (default), \code{"cividis"}, \code{"rocket"}, \code{"mako"}, \code{"turbo"}.
#' @param title plot title.  Automatically set if \code{NULL}.
#' @param z colour-bar label.  Automatically set if \code{NULL}.
#' @param xlab,ylab axis labels (default \code{"X/m"}, \code{"Y/m"}).
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' TT= random2d()
#' Plotparameter2d(TT,iflog=T)
#' TT= random2d(geo=list(me=0,var=1,geomod="Exp",anis=c(90,0.2),range=30,nugget=0))
#' Plotparameter2d(TT,iflog=T)

Plotparameter2d <- function(TT, palette="viridis",plotshow = TRUE,
                         plotfile = NULL, plotwidth = 10, plotheight = 10, iflog = FALSE,
                         title =NULL, z = NULL,xlab="X/m",ylab="Y/m") {

  if (iflog) {
    TT$Tp = log(TT$Tp)
    if(!is.null(title)) title = "2D log-transformed Transimisivity"
    if(!is.null(z))z = expression("log(T) ["*L^2*"/T]")
  } else {
    if(!is.null(title)) title = "2D Transimisivity"
    if(!is.null(z)) z = expression("T ["*L^2*"/T]")
  }

  require('ggplot2')
  p = ggplot(TT, aes(x, y, fill = Tp)) +
    geom_tile() +
    scale_fill_viridis_c(option = palette)+
    labs(title = title, x = xlab, y = ylab, fill = z)



  if(!is.null(plotfile)){
    ggsave(plotfile, width = plotwidth, height = plotheight)
  }
  if(plotshow) return(p)
}

#' Visualise 2D inversion results (mean/variance maps + scatter diagnostics)
#'
#' Creates a multi-panel diagnostic plot for a given iteration of the
#' ensemble-based inversion.  Panels 1 and 2 show the estimated ln(T) mean and
#' variance fields; panel 3 shows observed vs simulated drawdown; panel 4
#' (if \code{trueK} is provided) shows true vs estimated ln(T).
#'
#' @param niterm iteration number to visualise (index into \code{iterdf}).
#' @param grid grid list from \code{GenGrid()}.
#' @param iterdf list returned by \code{Finverse()}, \code{Finverse2()}, or
#'   \code{Finverse3()}.
#' @param oHT list of observation data frames used in the inversion (columns
#'   \code{data, x, y}).
#' @param trueK (optional) the true transmissivity field from
#'   \code{random2d()} — if supplied, a true-vs-estimated ln(T) scatter is
#'   added as panel 4.
#' @param p1title,p1z,p1xlab,p1ylab labels for panel 1 (mean ln(T) map).
#' @param p2title,p2z,p2xlab,p2ylab labels for panel 2 (variance map).
#' @param p3title,p3xlab,p3ylab labels for panel 3 (drawdown scatter).
#' @param p4title,p4xlab,p4ylab labels for panel 4 (ln(T) scatter).
#' @return A named list of ggplot objects: \code{lnk}, \code{varlnk},
#'   \code{headscatter}, and (if \code{trueK} is provided) \code{lnkscatter}.
#' @export
#' @examples
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT <- trueK[,-c(1,2)]
#' domain=c(40,40,0,40,0,40)
#' grid = GenGrid(domain)
#' Qinf1=data.frame(Qp=10,xp=20.5,yp=20.5)
#' qHT <- list(test1 = Qinf1)
#' trueh <- Fsteady2dsim(TT=TT,Qinf=Qinf1,grid=grid)
#' trueh <- trueh$solution
#' locx = c(15,18,22,25,30)
#' locy = c(15,18,22,25,30)
#' loc= expand.grid(x=locx,y=locy)
#' Oinf <- data.frame(data=NA,x=loc$x,y=loc$y)
#' Oinf <- getOelem(grid = grid)
#' Oinf$data <- trueh[Oinf$nelem]
#' oHT <- list(test1 = Oinf)
#' result <- Finverse(grid =grid, qHT = qHT, oHT = oHT)
#' inversePlot(niterm=5,grid=grid,iterdf = result,oHT = oHT,trueK=trueK)
inversePlot <- function(niterm=1,grid,iterdf,oHT,trueK=NULL,
                        p1title= "Est. lnT value",p1z= "lnT",p1xlab="X/m",p1ylab="Y/m",
                        p2title= 'Est. var. of lnT',p2z= "var. of lnT",p2xlab="X/m",p2ylab="Y/m",
                        p3title= "Drawdown scatter",p3xlab="Obs. Drawdown/m",p3ylab="Sim. Drawdown/m",
                        p4title= "lnT scatter",p4xlab="True lnT /(m2/d)",p4ylab="Est. lnT /(m2/d)"){
  # niterm = 5
  # note that trueh is actually a list contains many pumping tests.
  require(dplyr)
  oHTdf = bind_rows(oHT,.id='id')
  # np <- sapply(trueh,length) # get the number of observations for each test.
  # ntest <- rep("test1",np[1])
  # if(length(np)>=2){
  #   for (i in 2:length(np)){
  #     ntest <- append(ntest,rep(paste0("test",i),np[i]))
  #   }
  # }

  trueh <- oHTdf$data # the observation data in a vector.
  library(ggplot2)
  # build df first, before constructing p1
  if ("varT" %in% names(iterdf[[niterm]]) && !is.null(iterdf[[niterm]]$varT)) {
    df <- data.frame(x=grid$grid$x,y=grid$grid$y,v=iterdf[[niterm]]$meanT,var=iterdf[[niterm]]$varT)
    p2 <- ggplot(df) +
      aes(x = x, y = y, fill = var) +
      geom_tile()+
      scale_fill_viridis_c(option = "viridis", direction = 1)+
      labs(title = p2title, x = p2xlab, y = p2ylab, fill = p2z)
  } else {
    df <- data.frame(x=grid$grid$x,y=grid$grid$y,v=iterdf[[niterm]]$meanT)
    p2 <- NULL
  }
  p1 <- ggplot(df) +
    aes(x = x, y = y, fill = v) +
    geom_tile()+
    scale_fill_viridis_c(option = "viridis", direction = 1)+
    labs(title = p1title, x = p1xlab, y = p1ylab, fill = p1z)
  # change it to drawdown not head.
  has_varobsh <- "varobsh" %in% names(iterdf[[niterm]]) && !is.null(iterdf[[niterm]]$varobsh)
  if (has_varobsh) {
    dfh <- data.frame(trueh=-trueh,simh=-iterdf[[niterm]]$meanobsh,varh=iterdf[[niterm]]$varobsh,ntest=oHTdf$id)
  } else {
    dfh <- data.frame(trueh=-trueh,simh=-iterdf[[niterm]]$meanobsh,ntest=oHTdf$id)
  }
  #https://www.roelpeters.be/how-to-add-a-regression-equation-and-r-squared-in-ggplot2/
  #https://rpkgs.datanovia.com/ggpubr/reference/stat_regline_equation.html
  require(ggpubr)
  require(ggpp)
  require(tibble)
  xm = max(dfh[,1])
  xn = min(dfh[,1])
  ym = max(dfh[,2])
  yn = min(dfh[,2])
  x = xm- 0.15 *(xm - xn)
  y = yn + 0.01*(ym - yn)
  tb = tibble(round(statData(dfh),3))
  df = tibble(x=x,y=y,tb=list(tb))
  p3 <- ggplot(dfh) +
    geom_point(size=4,aes(x = trueh, y = simh,colour=ntest))+
    geom_smooth(method = "lm", aes(x = trueh, y = simh),se=FALSE) +
    stat_regline_equation(aes(x = trueh, y = simh),label.x = xn + 0.01*(xm - xn),label.y = ym- 0.01 *(ym - yn)) +
    geom_table(data=df,aes(x=x,y=y,label=tb)) +
    labs(title = p3title, x = p3xlab, y = p3ylab)
  if (has_varobsh) {
    p3 <- p3 + geom_errorbar(aes(x = trueh,ymin=simh-varh^0.5,ymax=simh+varh^0.5),width=0.02,color="gray")
  }

  p1
  p2
  p3

  if(!is.null(trueK)){
    # scatterks <- data.frame(x=log(trueK$Tp),y=log(iterdf[[niterm]]$meanT))
    scatterks <- data.frame(x=log(trueK$Tp),y=iterdf[[niterm]]$meanT) # inverse directly store lnT. 2026.01.01
    xm = max(scatterks[,1])
    xn = min(scatterks[,1])
    ym = max(scatterks[,2])
    yn = min(scatterks[,2])
    x = xm- 0.15 *(xm - xn)
    y = yn + 0.01*(ym - yn)
    tb = tibble(round(statData(scatterks),3))
    df = tibble(x=x,y=y,tb=list(tb))
    p4 <- ggplot(scatterks) +
      aes(x,y) +
      geom_point(size=1)+
      geom_smooth(method = "lm", se=T) +
      stat_regline_equation(label.x = xn + 0.01*(xm - xn),label.y = ym- 0.01 *(ym - yn))+
      geom_table(data=df,aes(x=x,y=y,label=tb)) +
      labs(title = p4title, x = p4xlab, y = p4ylab)
  }
  if(!is.null(trueK) && !is.null(p2))
    return(list(lnk=p1,varlnk=p2,lnkscatter=p4,headscatter=p3))
  else if(!is.null(trueK))
    return(list(lnk=p1,lnkscatter=p4,headscatter=p3))
  else if(!is.null(p2))
    return(list(lnk=p1,varlnk=p2,headscatter=p3))
  else
    return(list(lnk=p1,headscatter=p3))
}

# ---- 3D Visualization --------------------------------------------------------

#' Plot 3D hydraulic conductivity field from random3d()
#'
#' Visualises the 3D hydraulic conductivity (or log-K) field generated by
#' \code{\link{random3d}} using the \pkg{plot3D} package.  Three display types
#' are supported: mutually-perpendicular slices, isosurfaces, and a voxel plot.
#'
#' @param KK    data frame returned by \code{random3d()} (columns
#'   \code{x, y, z, Kp} or \code{sim1, …}).  If \code{nsim > 1} was used,
#'   the first simulation column is plotted.
#' @param grid  grid object from \code{GenGrid3D()}.
#' @param iflog logical; if \code{TRUE} (default) plot \eqn{\ln K}, otherwise
#'   plot \eqn{K} directly.
#' @param type  character, one of \code{"slice"} (default), \code{"isosurf"},
#'   or \code{"voxel"}.
#' @param xs,ys,zs  positions of the slice planes (passed to \code{slice3D}).
#'   Default: midpoint of each axis.  Ignored for \code{type = "isosurf"} /
#'   \code{"voxel"}.
#' @param level  isosurface level for \code{type = "isosurf"}.  Default: median
#'   of the plotted variable.
#' @param palette character; colour palette: \code{"jet"} (default),
#'   \code{"jet2"}, \code{"gg"}, or \code{"heat"}.
#' @param theta,phi  viewing angles (degrees) passed to \pkg{plot3D}.
#' @param title,clab  plot title and colour-key label.  Sensible defaults are
#'   set automatically.
#' @param plotfile optional file path; if supplied the plot is saved as a PDF.
#' @param plotwidth,plotheight  PDF dimensions in inches.
#' @param ...  additional arguments forwarded to the underlying \pkg{plot3D}
#'   function.
#' @return Invisibly returns \code{NULL}; the plot is drawn as a side-effect.
#' @export
#' @examples
#' grid3d <- GenGrid3D()
#' KK     <- random3d(nsim = 1, grid = grid3d)
#' Plotparameter3d(KK, grid3d)
#' Plotparameter3d(KK, grid3d, type = "isosurf")
Plotparameter3d <- function(KK,
                             grid,
                             iflog      = TRUE,
                             type       = "slice",
                             xs         = NULL,
                             ys         = NULL,
                             zs         = NULL,
                             level      = NULL,
                             palette    = "jet",
                             theta      = 40,
                             phi        = 25,
                             title      = NULL,
                             clab       = NULL,
                             plotfile   = NULL,
                             plotwidth  = 8,
                             plotheight = 7,
                             ...) {

  require('plot3D')

  # ---- extract first K column (4th column of the data frame) ----------------
  k_vec <- KK[[4]]
  V     <- if (iflog) log(k_vec) else k_vec
  V3    <- array(V, dim = c(grid$nx, grid$ny, grid$nz))

  xx <- grid$xmid; yy <- grid$ymid; zz <- grid$zmid

  # ---- defaults ---------------------------------------------------------------
  if (is.null(xs)) xs <- xx[round(grid$nx / 2)]
  if (is.null(ys)) ys <- yy[round(grid$ny / 2)]
  if (is.null(zs)) zs <- zz[round(grid$nz / 2)]
  if (is.null(title))
    title <- if (iflog) "3D ln(K) Field" else "3D Hydraulic Conductivity Field"
  if (is.null(clab))
    clab  <- if (iflog) "ln(K)" else expression("K ["*L*"/T]")

  col <- .hydro_col(palette, 100)

  # ---- open file device if requested -----------------------------------------
  if (!is.null(plotfile)) { pdf(plotfile, width = plotwidth, height = plotheight); on.exit(dev.off()) }

  # ---- dispatch on type -------------------------------------------------------
  if (type == "slice") {
    slice3D(xx, yy, zz, colvar = V3,
            xs = xs, ys = ys, zs = zs,
            col = col, theta = theta, phi = phi,
            main = title, clab = clab,
            xlab = "X (m)", ylab = "Y (m)", zlab = "Z (m)", ...)

  } else if (type == "isosurf") {
    # default: three levels at Q25 / Q50 / Q75 to show low/mid/high zones
    if (is.null(level)) level <- quantile(V, c(0.25, 0.5, 0.75))
    iso_col <- .hydro_col(palette, length(level))
    isosurf3D(xx, yy, zz, colvar = V3,
              level = level,
              col   = iso_col,
              shade = 0.3, alpha = 0.5,
              theta = theta, phi = phi,
              main = title,
              xlab = "X (m)", ylab = "Y (m)", zlab = "Z (m)", ...)

  } else if (type == "voxel") {
    # show cells whose value exceeds the median (highlights high-K zone)
    threshold <- if (is.null(level)) median(V) else level
    voxel3D(xx, yy, zz, colvar = V3,
            threshold = threshold,
            col = col,
            theta = theta, phi = phi,
            main = title, clab = clab,
            xlab = "X (m)", ylab = "Y (m)", zlab = "Z (m)", ...)

  } else {
    stop("'type' must be one of: \"slice\", \"isosurf\", \"voxel\"")
  }

  invisible(NULL)
}


#' Plot 3D steady-state head or drawdown field from Fsteady3dsim()
#'
#' Visualises the 3D head/drawdown field produced by
#' \code{\link{Fsteady3dsim}} using the \pkg{plot3D} package.  Pumping well
#' and observation well locations can be overlaid.
#'
#' @param s3   data frame returned by \code{Fsteady3dsim()} (columns
#'   \code{x, y, z, solution}).
#' @param grid grid object from \code{GenGrid3D()}.
#' @param ifdrawdown logical; if \code{TRUE} (default) negate the solution to
#'   display drawdown (positive values = water-level decline).
#' @param type  character, one of \code{"slice"} (default) or
#'   \code{"isosurf"}.
#' @param xs,ys,zs  slice-plane positions for \code{type = "slice"}.
#'   Default: midpoint of each axis.
#' @param level  isosurface level for \code{type = "isosurf"}.  Default: 25th
#'   percentile of non-zero values.
#' @param palette character; colour palette — \code{"jet"} (default),
#'   \code{"jet2"}, \code{"gg"}, or \code{"heat"}.
#' @param Qinf  data frame with columns \code{x, y, z} for pumping wells
#'   (plotted as red filled circles).  \code{NULL} to omit.
#' @param Oinf  data frame with columns \code{x, y, z} for observation wells
#'   (plotted as yellow filled squares).  \code{NULL} to omit.
#' @param theta,phi  viewing angles (degrees).
#' @param title,clab  plot title and colour-key label.
#' @param plotfile  optional PDF output path.
#' @param plotwidth,plotheight  PDF dimensions in inches.
#' @param ...  additional arguments forwarded to the \pkg{plot3D} function.
#' @return Invisibly returns \code{NULL}.
#' @export
#' @examples
#' grid3d <- GenGrid3D()
#' s3     <- Fsteady3dsim(grid = grid3d)
#' Plotsteady3d(s3, grid3d)
#' Qinf   <- data.frame(Qp = 10, x = 20.5, y = 20.5, z = 5.5)
#' Plotsteady3d(s3, grid3d, Qinf = Qinf, type = "isosurf")
Plotsteady3d <- function(s3,
                          grid,
                          ifdrawdown = TRUE,
                          type       = "slice",
                          xs         = NULL,
                          ys         = NULL,
                          zs         = NULL,
                          level      = NULL,
                          palette    = "jet",
                          Qinf       = NULL,
                          Oinf       = NULL,
                          theta      = 40,
                          phi        = 25,
                          title      = NULL,
                          clab       = NULL,
                          plotfile   = NULL,
                          plotwidth  = 8,
                          plotheight = 7,
                          ...) {

  require('plot3D')

  # ---- reshape solution to 3D array ------------------------------------------
  vals <- if (ifdrawdown) -s3$solution else s3$solution
  V3   <- array(vals, dim = c(grid$nx, grid$ny, grid$nz))

  xx <- grid$xmid; yy <- grid$ymid; zz <- grid$zmid

  # ---- defaults ---------------------------------------------------------------
  if (is.null(xs)) xs <- xx[round(grid$nx / 2)]
  if (is.null(ys)) ys <- yy[round(grid$ny / 2)]
  if (is.null(zs)) zs <- zz[round(grid$nz / 2)]
  if (is.null(title))
    title <- if (ifdrawdown) "3D Steady-State Drawdown" else "3D Steady-State Head"
  if (is.null(clab))
    clab  <- if (ifdrawdown) "Drawdown (m)" else "Head (m)"

  col <- .hydro_col(palette, 100)

  # ---- open file device if requested -----------------------------------------
  if (!is.null(plotfile)) { pdf(plotfile, width = plotwidth, height = plotheight); on.exit(dev.off()) }

  # ---- dispatch on type -------------------------------------------------------
  if (type == "slice") {
    slice3D(xx, yy, zz, colvar = V3,
            xs = xs, ys = ys, zs = zs,
            col = col, theta = theta, phi = phi,
            main = title, clab = clab,
            xlab = "X (m)", ylab = "Y (m)", zlab = "Z (m)", ...)

  } else if (type == "isosurf") {
    # default: one level at Q25 (captures drawdown cone near the well)
    nz_vals <- vals[vals != 0]
    if (is.null(level)) level <- quantile(nz_vals, 0.25)
    iso_col <- .hydro_col(palette, length(level))
    isosurf3D(xx, yy, zz, colvar = V3,
              level = level,
              col   = iso_col,
              shade = 0.3, alpha = 0.5,
              theta = theta, phi = phi,
              main = title,
              xlab = "X (m)", ylab = "Y (m)", zlab = "Z (m)", ...)

  } else {
    stop("'type' must be one of: \"slice\", \"isosurf\"")
  }

  # ---- overlay well locations ------------------------------------------------
  if (!is.null(Qinf)) {
    scatter3D(Qinf$x, Qinf$y, Qinf$z,
              add = TRUE, pch = 21, cex = 2,
              col = "red", bg = "red", colkey = FALSE)
  }
  if (!is.null(Oinf)) {
    scatter3D(Oinf$x, Oinf$y, Oinf$z,
              add = TRUE, pch = 22, cex = 1.5,
              col = "yellow", bg = "yellow", colkey = FALSE)
  }

  invisible(NULL)
}


# ---- internal helper: translate palette name to colour vector ----------------
.hydro_col <- function(palette, n) {
  switch(palette,
    "jet"  = jet.col(n),
    "jet2" = jet2.col(n),
    "gg"   = gg.col(n),
    "heat" = heat.colors(n),
    jet.col(n)  # fallback
  )
}


#' Scatter plot of observed vs simulated drawdown
#'
#' Creates a scatter plot of observed vs simulated drawdown with a linear
#' regression line, regression equation annotation, and a summary statistics
#' table (L1, RMSE, r, NSE from \code{\link{statData}}).
#'
#' @param df a data frame with columns \code{observed} and \code{simulated}.
#' @param title plot title (default \code{"Drawdown scatter"}).
#' @param xlab,ylab axis labels.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' df = data.frame(observed = rnorm(100), simulated = rnorm(100))
#' predictPlot(df)

predictPlot <- function(df,
                        title= "Drawdown scatter",
                        xlab="Obs. Drawdown/m",
                        ylab="Sim. Drawdown/m"
    )
{

library(ggplot2)
  require(ggpubr)
  require(ggpp)
  require(tibble)
  xm = max(df$observed)
  xn = min(df$observed)
  ym = max(df$simulated)
  yn = min(df$simulated)
  x = xm- 0.15 *(xm - xn)
  y = yn + 0.01*(ym - yn)
  tb = tibble(round(statData(df[,-1]),3))
  df1 = tibble(x=x,y=y,tb=list(tb))
  ggplot(df, aes(x=observed, y=simulated)) +
    geom_point(size=4) +
    labs(title=title, x=xlab, y=ylab) +
    geom_smooth(method = "lm", aes(x = observed, y = simulated),se=FALSE) +
    stat_regline_equation(aes(x = observed, y = simulated),label.x = xn + 0.01*(xm - xn),label.y = ym- 0.01 *(ym - yn)) +
    geom_table(data=df1,aes(x=x,y=y,label=tb)) +
    theme_minimal()
}


# ==============================================================================
# 3D Inversion Result Visualisation
# ==============================================================================

#' Plot 3D inversion results (horizontal z-slices + scatter diagnostics)
#'
#' 3D counterpart of \code{\link{inversePlot}}.  Panels 1 and 2 show
#' horizontal slices of the estimated ln(K) mean and variance fields at a
#' chosen z level; panels 3 and 4 show the same drawdown and ln(K) scatter
#' diagnostics as the 2D version.
#'
#' @param niterm  iteration number to visualise (index into \code{iterdf}).
#' @param grid    grid from \code{GenGrid3D()}.
#' @param iterdf  list returned by \code{Finverse3D()}.
#' @param oHT     list of observation data frames used in the inversion
#'   (columns \code{data, x, y, z}).
#' @param trueK   (optional) true K field from \code{random3d()} — if
#'   supplied, a true-vs-estimated ln(K) scatter is added as panel 4.
#' @param zslice  z coordinate of the horizontal slice to display in panels
#'   1 and 2.  Defaults to the middle layer of the grid.
#' @param p1title,p1z,p1xlab,p1ylab  labels for panel 1 (mean ln K map).
#' @param p2title,p2z,p2xlab,p2ylab  labels for panel 2 (variance map).
#' @param p3title,p3xlab,p3ylab  labels for panel 3 (drawdown scatter).
#' @param p4title,p4xlab,p4ylab  labels for panel 4 (ln K scatter).
#' @return a named list of ggplot objects:
#'   \code{lnk}, \code{varlnk}, \code{headscatter}, and (if \code{trueK}
#'   is provided) \code{lnkscatter}.
#' @export
#' @examples
#' # see Finverse3D() example for how to generate result3d
#' # figs <- inversePlot3D(niterm=1, grid=grid3d, iterdf=result3d,
#' #                       oHT=oHT3d, trueK=trueK3d)
#' # ggpubr::ggarrange(plotlist = figs, ncol=2, nrow=2)
inversePlot3D <- function(niterm  = 1,
                          grid,
                          iterdf,
                          oHT,
                          trueK   = NULL,
                          zslice  = NULL,
                          p1title = "Est. ln(K) — horizontal slice",
                          p1z     = "ln(K)",
                          p1xlab  = "X/m", p1ylab = "Y/m",
                          p2title = "Est. var. of ln(K) — horizontal slice",
                          p2z     = "var. of ln(K)",
                          p2xlab  = "X/m", p2ylab = "Y/m",
                          p3title = "Drawdown scatter",
                          p3xlab  = "Obs. Drawdown/m",
                          p3ylab  = "Sim. Drawdown/m",
                          p4title = "ln(K) scatter",
                          p4xlab  = "True ln(K)",
                          p4ylab  = "Est. ln(K)") {

  require(ggplot2)
  require(ggpubr)
  require(ggpp)
  require(tibble)
  require(dplyr)

  # ---- choose z slice -------------------------------------------------------
  if (is.null(zslice)) zslice <- grid$zmid[round(grid$nz / 2)]
  iz  <- which.min(abs(grid$zmid - zslice))
  # elements in this z-layer (x varies fastest in expand.grid)
  idx <- seq((iz - 1L) * grid$nx * grid$ny + 1L,
             iz        * grid$nx * grid$ny)

  # ---- panel 1: mean ln(K) slice (always present) ---------------------------
  if ("varT" %in% names(iterdf[[niterm]]) && !is.null(iterdf[[niterm]]$varT)) {
    df_slice <- data.frame(
      x   = grid$grid$x[idx],
      y   = grid$grid$y[idx],
      v   = iterdf[[niterm]]$meanT[idx],
      var = iterdf[[niterm]]$varT[idx]
    )
  } else {
    df_slice <- data.frame(
      x = grid$grid$x[idx],
      y = grid$grid$y[idx],
      v = iterdf[[niterm]]$meanT[idx]
    )
  }

  p1 <- ggplot(df_slice, aes(x = x, y = y, fill = v)) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis") +
    labs(title = paste0(p1title, " (z = ", round(zslice, 2), " m)"),
         x = p1xlab, y = p1ylab, fill = p1z)

  # ---- panel 2: variance slice (only if varT exists) ------------------------
  if ("varT" %in% names(iterdf[[niterm]]) && !is.null(iterdf[[niterm]]$varT)) {
    p2 <- ggplot(df_slice, aes(x = x, y = y, fill = var)) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma") +
      labs(title = paste0(p2title, " (z = ", round(zslice, 2), " m)"),
           x = p2xlab, y = p2ylab, fill = p2z)
  } else {
    p2 <- NULL
  }

  # ---- panel 3: drawdown scatter --------------------------------------------
  oHTdf  <- bind_rows(oHT, .id = 'id')
  trueh  <- oHTdf$data
  has_varobsh <- "varobsh" %in% names(iterdf[[niterm]]) && !is.null(iterdf[[niterm]]$varobsh)
  if (has_varobsh) {
    dfh <- data.frame(trueh  = -trueh,
                      simh   = -iterdf[[niterm]]$meanobsh,
                      varh   =  iterdf[[niterm]]$varobsh,
                      ntest  =  oHTdf$id)
  } else {
    dfh <- data.frame(trueh  = -trueh,
                      simh   = -iterdf[[niterm]]$meanobsh,
                      ntest  =  oHTdf$id)
  }
  xm <- max(dfh$trueh); xn <- min(dfh$trueh)
  ym <- max(dfh$simh);  yn <- min(dfh$simh)
  tb  <- tibble(round(statData(dfh[, 1:2]), 3))
  dft <- tibble(x = xm - 0.15 * (xm - xn),
                y = yn + 0.01  * (ym - yn),
                tb = list(tb))
  p3 <- ggplot(dfh) +
    geom_point(aes(x = trueh, y = simh, colour = ntest), size = 3) +
    geom_smooth(aes(x = trueh, y = simh), method = "lm", se = FALSE) +
    stat_regline_equation(aes(x = trueh, y = simh),
                          label.x = xn + 0.01 * (xm - xn),
                          label.y = ym - 0.01 * (ym - yn)) +
    geom_table(data = dft, aes(x = x, y = y, label = tb)) +
    labs(title = p3title, x = p3xlab, y = p3ylab)
  if (has_varobsh) {
    p3 <- p3 + geom_errorbar(aes(x = trueh,
                                 ymin = simh - varh^0.5,
                                 ymax = simh + varh^0.5),
                             width = 0.02, colour = "gray")
  }

  # ---- panel 4 (optional): true vs estimated ln(K) -------------------------
  if (!is.null(trueK)) {
    # trueK from random3d() — K values in 4th column
    scatterks <- data.frame(x = log(trueK[[4]]),
                            y = iterdf[[niterm]]$meanT)
    xm <- max(scatterks$x); xn <- min(scatterks$x)
    ym <- max(scatterks$y); yn <- min(scatterks$y)
    tb  <- tibble(round(statData(scatterks), 3))
    dft <- tibble(x = xm - 0.15 * (xm - xn),
                  y = yn + 0.01  * (ym - yn),
                  tb = list(tb))
    p4 <- ggplot(scatterks, aes(x, y)) +
      geom_point(size = 1) +
      geom_smooth(method = "lm", se = TRUE) +
      stat_regline_equation(label.x = xn + 0.01 * (xm - xn),
                            label.y = ym - 0.01 * (ym - yn)) +
      geom_table(data = dft, aes(x = x, y = y, label = tb)) +
      labs(title = p4title, x = p4xlab, y = p4ylab)

    if (!is.null(p2))
      return(list(lnk = p1, varlnk = p2, headscatter = p3, lnkscatter = p4))
    else
      return(list(lnk = p1, headscatter = p3, lnkscatter = p4))
  }

  if (!is.null(p2))
    return(list(lnk = p1, varlnk = p2, headscatter = p3))
  else
    return(list(lnk = p1, headscatter = p3))
}


# ==============================================================================
# Transient time-drawdown plots (2D / 3D)
# ==============================================================================

#' Plot time-drawdown curves for 2D transient flow simulation
#'
#' Given the output of \code{\link{Ftransient2dsim}} and one or more
#' observation locations \code{(x, y)}, extracts the head time-series at each
#' location and plots the drawdown (or head) evolution over time.
#'
#' @param result_tr list returned by \code{\link{Ftransient2dsim}}.
#' @param grid      grid list from \code{\link{GenGrid}}.
#' @param loc       a data frame with columns \code{x} and \code{y} giving the
#'   coordinates of observation points.  One curve is drawn per row.
#' @param label     optional character vector of legend labels, one per
#'   location.  If \code{NULL}, labels are auto-generated as
#'   \code{"(x, y)"}.
#' @param ifdrawdown logical; if \code{TRUE} (default), negate head to show
#'   drawdown (positive = water-level decline).
#' @param title     plot title.  Auto-set if \code{NULL}.
#' @param xlab,ylab axis labels (default \code{"Time"} and \code{"Drawdown / m"}
#'   or \code{"Head / m"}).
#' @param palette   character; \code{viridis} palette option for colour mapping.
#' @param plotfile  optional file path; if supplied the plot is saved as a PDF.
#' @param plotwidth,plotheight  PDF dimensions in inches.
#' @param linewidth line width for the curves (default 1).
#' @param iflogx  logical; if \code{TRUE}, use log10 scale for the time axis.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' grid  <- GenGrid()
#' res   <- Ftransient2dsim(grid = grid, times = seq(0, 1, by = 0.05))
#' loc   <- data.frame(x = c(10.5, 20.5, 30.5), y = c(20.5, 15.5, 25.5))
#' PlotTr2d(res, grid, loc)
#' PlotTr2d(res, grid, loc, label = c("Well A", "Well B", "Well C"), iflogx = TRUE)
PlotTr2d <- function(result_tr,
                     grid,
                     loc,
                     label       = NULL,
                     ifdrawdown  = TRUE,
                     iflogx      = FALSE,
                     title       = NULL,
                     xlab        = "Time",
                     ylab        = NULL,
                     palette     = "viridis",
                     plotfile    = NULL,
                     plotwidth   = 8,
                     plotheight  = 5,
                     linewidth   = 1) {

  require('ggplot2')

  # ---- map (x, y) → element numbers -----------------------------------------
  Oinf <- data.frame(data = NA, x = loc$x, y = loc$y)
  Oinf <- getOelem(grid = grid, Oinf = Oinf)
  nelems <- Oinf$nelem

  # ---- auto-generate labels --------------------------------------------------
  if (is.null(label)) {
    label <- paste0("(", round(loc$x, 2), ", ", round(loc$y, 2), ")")
  }

  # ---- extract time series from ODE output matrix ----------------------------
  sim_times <- result_tr$out[, 1]
  nloc      <- nrow(loc)
  nt        <- length(sim_times)

  df_list <- vector("list", nloc)
  for (i in seq_len(nloc)) {
    h   <- result_tr$out[, nelems[i] + 1L]  # +1 offset: col 1 is time
    val <- if (ifdrawdown) -h else h
    df_list[[i]] <- data.frame(
      time  = sim_times,
      value = val,
      label = label[i],
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, df_list)

  # ---- defaults ---------------------------------------------------------------
  if (is.null(title))
    title <- if (ifdrawdown) "2D Transient Drawdown" else "2D Transient Head"
  if (is.null(ylab))
    ylab  <- if (ifdrawdown) "Drawdown / m" else "Head / m"

  # ---- plot -------------------------------------------------------------------
  p <- ggplot(df, aes(x = time, y = value, colour = label)) +
    geom_line(linewidth = linewidth) +
    scale_colour_viridis_d(option = palette, end = 0.85) +
    labs(title = title, x = xlab, y = ylab, colour = "Location") +
    theme_minimal()

  if (iflogx) {
    p <- p + scale_x_log10()
  }

  if (!is.null(plotfile)) {
    ggsave(plotfile, width = plotwidth, height = plotheight)
  }

  return(p)
}


#' Plot time-drawdown curves for 3D transient flow simulation
#'
#' Given the output of \code{\link{Ftransient3dsim}} and one or more
#' observation locations \code{(x, y, z)}, extracts the head time-series at
#' each location and plots the drawdown (or head) evolution over time.
#'
#' @param result_tr list returned by \code{\link{Ftransient3dsim}}.
#' @param grid      grid list from \code{\link{GenGrid3D}}.
#' @param loc       a data frame with columns \code{x}, \code{y}, \code{z}
#'   giving the coordinates of observation points.  One curve per row.
#' @param label     optional character vector of legend labels, one per
#'   location.  If \code{NULL}, labels are auto-generated as
#'   \code{"(x, y, z)"}.
#' @param ifdrawdown logical; if \code{TRUE} (default), negate head to show
#'   drawdown (positive = water-level decline).
#' @param title     plot title.  Auto-set if \code{NULL}.
#' @param xlab,ylab axis labels (default \code{"Time"} and \code{"Drawdown / m"}
#'   or \code{"Head / m"}).
#' @param palette   character; \code{viridis} palette option for colour mapping.
#' @param plotfile  optional file path; if supplied the plot is saved as a PDF.
#' @param plotwidth,plotheight  PDF dimensions in inches.
#' @param linewidth line width for the curves (default 1).
#' @param iflogx  logical; if \code{TRUE}, use log10 scale for the time axis.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' grid3d <- GenGrid3D()
#' res3d  <- Ftransient3dsim(grid = grid3d, times = seq(0, 1, by = 0.05))
#' loc    <- data.frame(x = c(10.5, 25.5), y = c(20.5, 20.5), z = c(5.5, 5.5))
#' PlotTr3d(res3d, grid3d, loc)
#' PlotTr3d(res3d, grid3d, loc, label = c("Shallow", "Deep"), iflogx = TRUE)
PlotTr3d <- function(result_tr,
                     grid,
                     loc,
                     label       = NULL,
                     ifdrawdown  = TRUE,
                     iflogx      = FALSE,
                     title       = NULL,
                     xlab        = "Time",
                     ylab        = NULL,
                     palette     = "viridis",
                     plotfile    = NULL,
                     plotwidth   = 8,
                     plotheight  = 5,
                     linewidth   = 1) {

  require('ggplot2')

  # ---- map (x, y, z) → element numbers ---------------------------------------
  Oinf <- data.frame(data = NA, x = loc$x, y = loc$y, z = loc$z)
  Oinf <- getOelem3D(grid = grid, Oinf = Oinf)
  nelems <- Oinf$nelem

  # ---- auto-generate labels --------------------------------------------------
  if (is.null(label)) {
    label <- paste0("(", round(loc$x, 2), ", ",
                    round(loc$y, 2), ", ",
                    round(loc$z, 2), ")")
  }

  # ---- extract time series from ODE output matrix ----------------------------
  sim_times <- result_tr$out[, 1]
  nloc      <- nrow(loc)
  nt        <- length(sim_times)

  df_list <- vector("list", nloc)
  for (i in seq_len(nloc)) {
    h   <- result_tr$out[, nelems[i] + 1L]  # +1 offset: col 1 is time
    val <- if (ifdrawdown) -h else h
    df_list[[i]] <- data.frame(
      time  = sim_times,
      value = val,
      label = label[i],
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, df_list)

  # ---- defaults ---------------------------------------------------------------
  if (is.null(title))
    title <- if (ifdrawdown) "3D Transient Drawdown" else "3D Transient Head"
  if (is.null(ylab))
    ylab  <- if (ifdrawdown) "Drawdown / m" else "Head / m"

  # ---- plot -------------------------------------------------------------------
  p <- ggplot(df, aes(x = time, y = value, colour = label)) +
    geom_line(linewidth = linewidth) +
    scale_colour_viridis_d(option = palette, end = 0.85) +
    labs(title = title, x = xlab, y = ylab, colour = "Location") +
    theme_minimal()

  if (iflogx) {
    p <- p + scale_x_log10()
  }

  if (!is.null(plotfile)) {
    ggsave(plotfile, width = plotwidth, height = plotheight)
  }

  return(p)
}


# ==============================================================================
# Cross-covariance (covhk) visualisation — 2D and 3D
# ==============================================================================

#' Plot cross-covariance (sensitivity) map for a single 2D observation
#'
#' Given a \code{covhk} matrix returned by any parallel 2D inversion function
#' (e.g. \code{Finverse3}, \code{Finverse3Tr}) with \code{ifcor = TRUE}, this
#' function plots the cross-covariance between one selected observation's head
#' and the log-transmissivity at every grid cell.  Pumping-well and
#' observation-well locations are overlaid.
#'
#' @param covhk  \code{nobs × n} cross-covariance matrix (heads × ln-T),
#'   returned by a \code{Finverse*} function with \code{ifcor = TRUE}.
#' @param grid   grid list from \code{GenGrid()}.
#' @param iobs   integer; which observation row of \code{covhk} to plot
#'   (1-based index into the concatenated observation vector).
#' @param qHT    list of pumping-test data frames (columns \code{Qp, x, y}),
#'   used to overlay well locations and auto-label.
#' @param oHT    list of observation data frames (columns \code{data, x, y}),
#'   used to overlay well locations and auto-label.
#' @param labels  optional character vector of custom labels for observations.
#'   If \code{NULL} (default), labels are auto-generated as \code{"P1 (x,y)"},
#'   \code{"P2 (x,y)"}, … for pumping wells and \code{"O (x,y)"} for
#'   observation wells.
#' @param palette  \code{viridis} palette option (default \code{"viridis"}).
#' @param title    plot title; auto-generated if \code{NULL}.
#' @param zlab     colour-bar label; default \code{"Cov(h, lnT)"}.
#' @param xlab,ylab axis labels.
#' @param plotfile optional PDF file path.
#' @param plotwidth,plotheight PDF dimensions in inches.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#' covhk <- Finverse3(grid = grid, qHT = qHT, oHT = oHT, ifcor = TRUE)
#' covhkPlot2D(covhk, grid, iobs = 1, qHT = qHT, oHT = oHT)
#' }
covhkPlot2D <- function(covhk,
                         grid,
                         iobs       = 1L,
                         qHT        = NULL,
                         oHT        = NULL,
                         labels     = NULL,
                         palette    = "viridis",
                         title      = NULL,
                         zlab       = "Cov(h, lnT)",
                         xlab       = "X/m",
                         ylab       = "Y/m",
                         plotfile   = NULL,
                         plotwidth  = 8,
                         plotheight = 7) {

  require('ggplot2')

  n <- grid$nx * grid$ny
  if (ncol(covhk) != n) {
    stop("ncol(covhk) must equal grid$nx * grid$ny")
  }
  if (iobs < 1L || iobs > nrow(covhk)) {
    stop("iobs out of range [1, ", nrow(covhk), "]")
  }

  # ---- build data frame for the sensitivity field ----------------------------
  vals <- as.vector(covhk[iobs, ])
  df   <- data.frame(x = grid$grid$x, y = grid$grid$y, cov = vals)

  # ---- build observation-index lookup ----------------------------------------
  # Each observation in oHT is a row; concatenated across all tests.
  obs_info <- .build_obs_info_2d(oHT, qHT, labels)

  # ---- auto title ------------------------------------------------------------
  if (is.null(title)) {
    title <- paste0("Cross-covariance: obs ", iobs,
                    " [", obs_info$label[iobs], "]")
  }

  # ---- plot ------------------------------------------------------------------
  p <- ggplot(df, aes(x = x, y = y, fill = cov)) +
    geom_tile() +
    scale_fill_viridis_c(option = palette) +
    labs(title = title, x = xlab, y = ylab, fill = zlab) +
    coord_fixed()

  # ---- overlay pumping wells -------------------------------------------------
  if (!is.null(qHT)) {
    qdf <- .flatten_qHT_2d(qHT, labels)
    p <- p +
      geom_point(data = qdf, aes(x = x, y = y),
                 color = "red", fill = "white", shape = 21, size = 3,
                 inherit.aes = FALSE) +
      geom_text(data = qdf, aes(x = x, y = y, label = label),
                color = "red", vjust = -1, size = 3, inherit.aes = FALSE)
  }

  # ---- overlay observation wells ---------------------------------------------
  if (!is.null(oHT)) {
    odf <- obs_info
    # highlight the selected observation
    odf$highlight <- seq_len(nrow(odf)) == iobs
    p <- p +
      geom_point(data = odf, aes(x = x, y = y,
                                  color = highlight, shape = highlight),
                 fill = "yellow", size = 2.5, inherit.aes = FALSE) +
      scale_color_manual(values = c("FALSE" = "yellow", "TRUE" = "cyan"),
                         guide = "none") +
      scale_shape_manual(values = c("FALSE" = 22, "TRUE" = 23),
                         guide = "none") +
      geom_text(data = odf, aes(x = x, y = y, label = label),
                color = "yellow3", vjust = -1, size = 2.5, inherit.aes = FALSE)
  }

  if (!is.null(plotfile)) {
    ggsave(plotfile, width = plotwidth, height = plotheight)
  }

  return(p)
}


# ---- Internal: build observation info for 2D ----------------------------------
.build_obs_info_2d <- function(oHT, qHT, labels) {
  if (is.null(oHT)) {
    return(data.frame(x = numeric(0), y = numeric(0),
                      label = character(0), stringsAsFactors = FALSE))
  }
  require(dplyr)
  odf_all <- bind_rows(oHT, .id = "test_id")

  # determine x/y column names
  xy_cols <- if ("xp" %in% names(odf_all)) c("xp", "yp") else c("x", "y")
  # for Qinf data frames the columns are Qp, xp, yp (steady) or Qp, x, y (transient)
  q_xy <- if (!is.null(qHT)) {
    q1 <- qHT[[1]]
    if ("xp" %in% names(q1)) c("xp", "yp") else c("x", "y")
  } else {
    c("x", "y")
  }

  # observation well positions
  x <- odf_all[[xy_cols[1]]]
  y <- odf_all[[xy_cols[2]]]

  if (is.null(labels)) {
    label <- paste0("O(", round(x, 1), ",", round(y, 1), ")")
  } else {
    label <- labels
  }

  data.frame(x = x, y = y, label = label, stringsAsFactors = FALSE)
}


# ---- Internal: flatten qHT to a data frame with labels -----------------------
.flatten_qHT_2d <- function(qHT, labels) {
  if (is.null(qHT)) return(data.frame())
  ntest <- length(qHT)
  xy_cols <- if ("xp" %in% names(qHT[[1]])) c("xp", "yp") else c("x", "y")

  x <- numeric(ntest)
  y <- numeric(ntest)
  for (i in seq_len(ntest)) {
    x[i] <- qHT[[i]][[xy_cols[1]]][1]
    y[i] <- qHT[[i]][[xy_cols[2]]][1]
  }

  if (is.null(labels)) {
    label <- paste0("P", seq_len(ntest), "(", round(x, 1), ",", round(y, 1), ")")
  } else {
    label <- labels[seq_len(ntest)]
  }

  data.frame(x = x, y = y, label = label, stringsAsFactors = FALSE)
}


# ==============================================================================
# 3D covhk plotting
# ==============================================================================

#' Plot cross-covariance (sensitivity) slice for a single 3D observation
#'
#' 3D counterpart of \code{\link{covhkPlot2D}}.  Given a \code{covhk} matrix
#' from \code{Finverse3D}, \code{Finverse3DTr}, or \code{Finverse3DTrScreen}
#' (all with \code{ifcor = TRUE}), this function plots a horizontal z-slice of
#' the cross-covariance between one observation's head and log-K at every grid
#' cell.  Pumping and observation well locations (projected onto the slice) are
#' overlaid.
#'
#' @param covhk   \code{nobs × n} cross-covariance matrix.
#' @param grid    grid from \code{GenGrid3D()}.
#' @param iobs    integer; which observation row to plot.
#' @param zslice  z coordinate \code{[L]} of the horizontal slice to display.
#'   Default: middle layer of the grid.
#' @param qHT     list of pumping-test data frames (columns
#'   \code{Qp, x, y, z} or \code{Qp, x, y, z_top, z_bottom}).
#' @param oHT     list of observation data frames (columns
#'   \code{data, x, y, z} or \code{data, x, y, z_top, z_bottom}).
#' @param labels  optional custom labels.  Auto-generated if \code{NULL}.
#' @param palette  \code{viridis} palette option.
#' @param title    plot title.
#' @param zlab     colour-bar label.
#' @param xlab,ylab axis labels.
#' @param plotfile optional PDF path.
#' @param plotwidth,plotheight PDF dimensions.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#' covhk <- Finverse3D(grid = grid3d, qHT = qHT, oHT = oHT, ifcor = TRUE)
#' covhkPlot3D(covhk, grid3d, iobs = 1, qHT = qHT, oHT = oHT)
#' }
covhkPlot3D <- function(covhk,
                         grid,
                         iobs       = 1L,
                         zslice     = NULL,
                         qHT        = NULL,
                         oHT        = NULL,
                         labels     = NULL,
                         palette    = "viridis",
                         title      = NULL,
                         zlab       = "Cov(h, lnK)",
                         xlab       = "X/m",
                         ylab       = "Y/m",
                         plotfile   = NULL,
                         plotwidth  = 8,
                         plotheight = 7) {

  require('ggplot2')

  n <- grid$nx * grid$ny * grid$nz
  if (ncol(covhk) != n) {
    stop("ncol(covhk) must equal grid$nx * grid$ny * grid$nz")
  }
  if (iobs < 1L || iobs > nrow(covhk)) {
    stop("iobs out of range [1, ", nrow(covhk), "]")
  }

  # ---- choose z slice --------------------------------------------------------
  if (is.null(zslice)) zslice <- grid$zmid[round(grid$nz / 2)]
  iz <- which.min(abs(grid$zmid - zslice))
  # elements: x varies fastest, then y, then z
  idx <- seq((iz - 1L) * grid$nx * grid$ny + 1L,
             iz        * grid$nx * grid$ny)

  vals <- as.vector(covhk[iobs, idx])
  df   <- data.frame(x   = grid$grid$x[idx],
                     y   = grid$grid$y[idx],
                     cov = vals)

  # ---- build observation info ------------------------------------------------
  obs_info <- .build_obs_info_3d(oHT, qHT, labels)

  # ---- auto title ------------------------------------------------------------
  if (is.null(title)) {
    title <- paste0("Cross-covariance: obs ", iobs,
                    " [", obs_info$label[iobs], "]",
                    "  (z = ", round(zslice, 2), " m)")
  }

  # ---- plot ------------------------------------------------------------------
  p <- ggplot(df, aes(x = x, y = y, fill = cov)) +
    geom_tile() +
    scale_fill_viridis_c(option = palette) +
    labs(title = title, x = xlab, y = ylab, fill = zlab) +
    coord_fixed()

  # ---- overlay pumping wells -------------------------------------------------
  if (!is.null(qHT)) {
    qdf <- .flatten_qHT_3d(qHT, labels)
    p <- p +
      geom_point(data = qdf, aes(x = x, y = y),
                 color = "red", fill = "white", shape = 21, size = 3,
                 inherit.aes = FALSE) +
      geom_text(data = qdf, aes(x = x, y = y, label = label),
                color = "red", vjust = -1, size = 3, inherit.aes = FALSE)
  }

  # ---- overlay observation wells ---------------------------------------------
  if (!is.null(oHT)) {
    odf <- obs_info
    odf$highlight <- seq_len(nrow(odf)) == iobs
    p <- p +
      geom_point(data = odf, aes(x = x, y = y,
                                  color = highlight, shape = highlight),
                 fill = "yellow", size = 2.5, inherit.aes = FALSE) +
      scale_color_manual(values = c("FALSE" = "yellow", "TRUE" = "cyan"),
                         guide = "none") +
      scale_shape_manual(values = c("FALSE" = 22, "TRUE" = 23),
                         guide = "none") +
      geom_text(data = odf, aes(x = x, y = y, label = label),
                color = "yellow3", vjust = -1, size = 2.5, inherit.aes = FALSE)
  }

  if (!is.null(plotfile)) {
    ggsave(plotfile, width = plotwidth, height = plotheight)
  }

  return(p)
}


# ---- Internal: build observation info for 3D ---------------------------------
.build_obs_info_3d <- function(oHT, qHT, labels) {
  if (is.null(oHT)) {
    return(data.frame(x = numeric(0), y = numeric(0),
                      label = character(0), stringsAsFactors = FALSE))
  }
  require(dplyr)
  odf_all <- bind_rows(oHT, .id = "test_id")

  # determine x/y column names (robust to steady vs transient vs screen)
  xcol <- if ("xp" %in% names(odf_all)) "xp" else "x"
  ycol <- if ("yp" %in% names(odf_all)) "yp" else "y"

  x <- odf_all[[xcol]]
  y <- odf_all[[ycol]]

  if (is.null(labels)) {
    label <- paste0("O(", round(x, 1), ",", round(y, 1), ")")
  } else {
    label <- labels
  }

  data.frame(x = x, y = y, label = label, stringsAsFactors = FALSE)
}


# ==============================================================================
# Jacobian (sensitivity) visualisation — 2D steady & transient
# ==============================================================================

#' Plot Jacobian sensitivity maps for 2D hydraulic tomography
#'
#' Visualises the sensitivity of each observation's head to the
#' log-transmissivity \eqn{\ln T} at every grid cell.  Each observation
#' is shown as a separate tile-map panel, faceted so that many observations
#' can be inspected at once.  Works with the output of
#' \code{\link{jacobian2D}} (steady-state) and \code{\link{jacobian2DTr}}
#' (transient).
#'
#' Positive sensitivity (red) means increasing \eqn{T} at that cell raises
#' the observed head; negative sensitivity (blue) means it lowers it.
#'
#' @param J      \code{nobs × n} Jacobian matrix from
#'   \code{jacobian2D()} or \code{jacobian2DTr()}.
#' @param grid   grid list from \code{GenGrid()}.
#' @param oHT    list of observation data frames used to compute \code{J}
#'   (columns \code{data, x, y} for steady, plus \code{time} for transient).
#'   If supplied, observation labels include coordinates and times.
#' @param maxobs maximum number of observation panels to display (default
#'   \code{16}).  Observations beyond \code{maxobs} are omitted to keep the
#'   plot readable.  Set to \code{Inf} to show all.
#' @param ncol   number of facet columns (default \code{4}).
#' @param palette \code{viridis} palette option (default \code{"RdBu"} —
#'   uses \code{ggplot2}'s built-in diverging palette, not \code{viridis}).
#' @param title  plot title.  Auto-generated if \code{NULL}.
#' @param zlab   colour-bar label.  Default \code{expression(partial(h) /
#'   partial(lnT))}.
#' @param xlab,ylab axis labels (default \code{"X/m"}, \code{"Y/m"}).
#' @param plotfile optional PDF file path.
#' @param plotwidth,plotheight PDF dimensions in inches.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#' grid <- GenGrid(c(20, 20, 0, 20, 0, 20))
#' set.seed(123)
#' TT <- random2d(nsim = 1, grid = grid)$Tp
#' qHT <- list(data.frame(Qp = 10, x = 10.5, y = 10.5))
#' oHT <- list(data.frame(data = -1, x = 5.5, y = 5.5))
#' J   <- jacobian2D(grid = grid, TT = TT, qHT = qHT, oHT = oHT)
#' PlotJ2d(J, grid, oHT)
#'
#' # transient example
#' oHTt <- list(data.frame(data = -1, x = 5.5, y = 5.5, time = 50))
#' Jt   <- jacobian2DTr(grid = grid, TT = TT, SS = 1e-4,
#'                      qHT = qHT, oHT = oHTt)
#' PlotJ2d(Jt, grid, oHTt)
#' }
PlotJ2d <- function(J,
                     grid,
                     oHT       = NULL,
                     maxobs    = 16,
                     ncol      = 4,
                     palette   = "RdBu",
                     title     = NULL,
                     zlab      = NULL,
                     xlab      = "X/m",
                     ylab      = "Y/m",
                     plotfile  = NULL,
                     plotwidth = 10,
                     plotheight = 8) {

  require('ggplot2')

  n   <- grid$n
  nx  <- grid$nx
  ny  <- grid$ny
  nobs <- nrow(J)

  if (ncol(J) != n) {
    stop("ncol(J) = ", ncol(J), " does not match grid$n = ", n)
  }

  # ---- auto-generate observation labels ---------------------------------------
  obs_label <- paste0("Obs", seq_len(nobs))
  if (!is.null(oHT)) {
    require(dplyr)
    odf <- dplyr::bind_rows(oHT, .id = "test")
    if (nrow(odf) != nobs) {
      warning("nrow(flattened oHT) != nrow(J); using generic observation labels")
    } else {
      if ("time" %in% names(odf)) {
        # transient
        obs_label <- paste0("(", round(odf$x, 1), ", ",
                            round(odf$y, 1), ", t=", round(odf$time, 1), ")")
      } else {
        # steady
        obs_label <- paste0("(", round(odf$x, 1), ", ", round(odf$y, 1), ")")
      }
      # if multiple tests, prefix with test ID
      if (length(unique(odf$test)) > 1) {
        obs_label <- paste0(odf$test, " ", obs_label)
      }
    }
  }

  # ---- limit number of panels -------------------------------------------------
  if (nobs > maxobs) {
    message("Showing ", maxobs, " of ", nobs, " observations. Set maxobs = Inf to show all.")
    nobs_show <- maxobs
  } else {
    nobs_show <- nobs
  }

  # ---- build long-format data frame -------------------------------------------
  df_list <- vector("list", nobs_show)
  for (i in seq_len(nobs_show)) {
    df_list[[i]] <- data.frame(
      x        = grid$grid$x,
      y        = grid$grid$y,
      sens     = as.vector(J[i, ]),
      obs      = factor(obs_label[i], levels = obs_label[seq_len(nobs_show)]),
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, df_list)

  # ---- defaults ---------------------------------------------------------------
  if (is.null(title)) {
    title <- if ("time" %in% names(oHT[[1]])) {
      "2D Transient Jacobian Sensitivity"
    } else {
      "2D Steady-State Jacobian Sensitivity"
    }
  }
  if (is.null(zlab)) {
    zlab <- expression(partial(h) / partial(lnT))
  }

  # ---- diverging colour scale -------------------------------------------------
  # find symmetric limits
  abs_max <- max(abs(df$sens), na.rm = TRUE)
  if (abs_max < 1e-12) abs_max <- 1  # avoid zero range

  # ---- plot -------------------------------------------------------------------
  p <- ggplot(df, aes(x = x, y = y, fill = sens)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, limits = c(-abs_max, abs_max)) +
    facet_wrap(~ obs, ncol = ncol) +
    labs(title = title, x = xlab, y = ylab, fill = zlab) +
    coord_fixed() +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(size = 7))

  if (!is.null(plotfile)) {
    ggsave(plotfile, width = plotwidth, height = plotheight)
  }

  return(p)
}


# ---- Internal: flatten 3D qHT to a data frame with labels --------------------
.flatten_qHT_3d <- function(qHT, labels) {
  if (is.null(qHT)) return(data.frame())
  ntest <- length(qHT)
  xy_cols <- if ("xp" %in% names(qHT[[1]])) c("xp", "yp") else c("x", "y")

  x <- numeric(ntest)
  y <- numeric(ntest)
  for (i in seq_len(ntest)) {
    x[i] <- qHT[[i]][[xy_cols[1]]][1]
    y[i] <- qHT[[i]][[xy_cols[2]]][1]
  }

  if (is.null(labels)) {
    label <- paste0("P", seq_len(ntest), "(", round(x, 1), ",", round(y, 1), ")")
  } else {
    label <- labels[seq_len(ntest)]
  }

  data.frame(x = x, y = y, label = label, stringsAsFactors = FALSE)
}
