#' 2D Ensemble-based Bayesian inversion for hydraulic tomography (serial)
#'
#' Performs ensemble-smoother-based Bayesian inversion to estimate the 2D
#' log-transmissivity field \eqn{\ln T(\mathbf{x})} from hydraulic tomography
#' (HT) data.  The algorithm iteratively updates an ensemble of \eqn{T} fields
#' using the ensemble Kalman-like update equation:
#' \deqn{\mathbf{K}_i^{new} = \mathbf{K}_i^{old} \odot
#' \exp\left( \mathbf{a}^\top \cdot (\mathbf{h}_{obs} - \mathbf{h}_i^{sim}) \right)}
#' where the gain matrix is \eqn{\mathbf{a} = \mathbf{C}_{hh}^{-1} \mathbf{C}_{hk}},
#' \eqn{\mathbf{C}_{hh}} is the ensemble covariance of simulated heads at
#' observation locations, and \eqn{\mathbf{C}_{hk}} is the cross-covariance
#' between simulated heads and log-transmissivity.  A stabiliser \code{mul} is
#' applied to the diagonal of \eqn{\mathbf{C}_{hh}} and decays by a factor
#' \code{decay} each iteration.
#'
#' This version uses serial \code{lapply} for forward simulations.  For
#' parallel execution see \code{\link{Finverse3}}.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.  \code{nx}
#'   and \code{ny} are cell counts; \code{x1:x2} and \code{y1:y2} are domain
#'   extents \code{[L]}.
#' @param grid grid list from \code{GenGrid()}; generated from \code{domain} if
#'   \code{NULL}.
#' @param qHT list of pumping test data frames, each with columns \code{Qp}
#'   \code{[L\eqn{^3}/T]}, \code{x} \code{[L]}, \code{y} \code{[L]}.
#' @param nsim ensemble size (number of realisations).  Default 50; use smaller
#'   values for testing.
#' @param itermax maximum number of iterations.
#' @param varmeanTmax convergence threshold on the spatial variance of mean
#'   ln(T).  Iteration stops when this variance exceeds the threshold.
#' @param rmsemin minimum weighted RMSE stopping criterion.  Iteration stops
#'   when RMSE falls below this value.
#' @param mul stabiliser coefficient (default 1.0).  Added to the diagonal of
#'   the head covariance matrix to improve numerical stability.  The value
#'   decays by \code{decay} each iteration.
#' @param decay stabiliser decay factor per iteration (default 1.05).
#'   \code{mul <- mul / decay} at the end of each iteration.
#' @param oHT list of observation data frames, each with columns \code{data}
#'   \code{[L]}, \code{x} \code{[L]}, \code{y} \code{[L]}.
#' @param lrw integer, the length of the real work array for \code{steady.2D}.
#'   Default 160000; increase for larger grids.
#' @param geo prior geostatistical parameters — same list structure as
#'   \code{random2d()}: \code{list(me, var, geomod, anis, range, nugget)}.
#'   \code{me} and \code{var} are the mean and variance of ln(T);
#'   \code{geomod} is the variogram model type (e.g. \code{"Exp"});
#'   \code{anis} is \code{c(azimuth_deg, ratio)}; \code{range} is the
#'   correlation range \code{[L]}; \code{nugget} is the nugget variance.
#' @param ifcor    logical; if \code{TRUE}, returns the ensemble results after
#'   the first iteration without performing the update. Useful for checking#'   forward simulations before running the full inversion.
#' @return A list of per-iteration results, each element a list with
#'   \code{meanT} (mean ln(T) vector), \code{varT} (variance of ln(T)),
#'   \code{meanobsh} (mean simulated head), \code{varobsh} (variance of
#'   simulated head).  When \code{ifcor = TRUE}, returns list with
#'   \code{obsh} (simulated heads), \code{Tnew} (ensemble), and
#'   \code{meanobsh}.
#' @export
#' @examples
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT <- trueK[,-c(1,2)]
#' domain=c(40,40,0,40,0,40)
#' grid = GenGrid(domain)
#' Qinf1=data.frame(Qp=10,xp=20.5,yp=20.5)
#' qHT <- list(test1 = Qinf1)
#' trueh <- Fsteady2dsim(TT=TT,Qinf=Qinf1,domain=domain)
#' locx = c(15,18,22,25,30)
#' locy = c(15,18,22,25,30)
#' loc= expand.grid(x=locx,y=locy)
#' Oinf1 <- data.frame(data=NA,x=loc$x,y=loc$y)
#' Oinf1 <- samData(grid = grid,Oinf = Oinf1,h=trueh$solution)
#' oHT <- list(test1 = Oinf1)
#' result <- Finverse(grid =grid, qHT = qHT, oHT = oHT)
#' inversePlot(niterm=5,iterdf = result,oHT = oHT,trueK=trueK,grid=grid)
#' Qinf2=data.frame(Qp=10,xp=10.5,yp=10.5)
#' trueh2 <- Fsteady2dsim(TT=TT,Qinf=Qinf2,domain=domain)
#' Oinf2 <- Oinf1
#' Oinf2 <- samData(grid = grid,Oinf = Oinf2,h=trueh2$solution)
#' oHT <- list(test1 = Oinf1,test2 = Oinf2)
#' qHT <- list(test1 = Qinf1,test2 = Qinf2)
#' result <- Finverse(grid =grid, qHT = qHT, oHT = oHT)
#' inversePlot(niterm=5,iterdf = result,oHT = oHT,trueK=trueK,grid=grid)

Finverse <- function(
    domain=c(40,40,0,40,0,40),
    grid = NULL,
    qHT=list(data.frame(Qp=10,x=20.5,y=20.5)), # should be list since it is HT.
    nsim=50,
    itermax=5,
    varmeanTmax =5,
    rmsemin = 0,
    mul=1.0,
    decay = 1.05,
    oHT= list(data.frame(data=-1,x=11,y=11)),
    lrw=160000,
    geo=list(me=0,var=1,geomod="Exp",anis=c(90,1),range=30,nugget=0), # should be list since multiple pumping test.
    ifcor = FALSE) {

  #1. ----- before iteration (generating ensemble.....)----------
  set.seed(200)
  ### record the time so to see how long it takes.
  startTime <- Sys.time()
  if(is.null(grid))grid = GenGrid(domain)

  #nsim = 50
  nHT <- length(qHT) # number of pumping test.also list length of oHT.
  trueobshHT <- list()
  loc_obsHT <- list()
  for (i in 1:nHT){
    oinf <- oHT[[i]]
    oelemdf <- getOelem(domain = domain,grid = grid,Oinf = oinf)
    loc_obsHT[[i]] <- oelemdf$nelem
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT) # the data format in HT is a list.
  nobs <- length(trueobsh)

  yy <- random2d(nsim=nsim,grid = grid, geo = geo)
  # initial ensemble.
  Tnew <- yy[,-c(1,2)]
  # nelem*nsim
  ### get the variance map of Tnew.
  ### this should goes with iterations.
  niter <- 1
  varmeanT <- 0
  rmse <- 1e10
  #itermax <- 100
  #varmeanTmax <- 5
  #rmsemin <- 0
  #mul <- 1.0 # stablizer.
  msgdf <- data.frame(niter = niter,varmeanT=varmeanT, rmse=rmse)
  iterdf <- list()
  #1. ----- before iteration----------0.94 0.00 0.94 (blh_synthetic)
  # start of the itertion loop.
  while(niter<=itermax & varmeanT<varmeanTmax & rmse>rmsemin){
  #2. -----------get the obsh......------
     varT <- apply(log(Tnew),1,var)
    meanT <- apply(log(Tnew),1,mean)
    varmeanT = var(meanT)
    # from Tnew to h.

    # hHT <- list()
    # for (j in 1:nHT){
    #   Qinf <- qHT[[j]]  # the information of jth pumping test.
    #   h <- Fsteady2dsim(TT=Tnew[,1],Qinf=Qinf)$solution
    #   for(i in 2:nsim){
    #   h <- rbind(h,Fsteady2dsim(TT=Tnew[,i],Qinf=Qinf)$solution)
    #   }
    #   # h: nsim*nelem
    #   obsh<- h[,loc_obsHT[[j]]]
    #   hHT[[j]] <- obsh # store the j th pumping test observations.
    # }
    ### upgrade to parallel computing...........
     # require(parallel)
     # n = parallel::detectCores() - 1
     # print(paste("using cores of ", n))
     # cl = parallel::makeCluster(n)

    hHT <- lapply(1:nHT, function(j) {
      Qinf <- qHT[[j]]
      # 使用 lapply 生成所有模拟结果
       h_list <- lapply(1:nsim, function(i){
         tmp = Fsteady2dsim(grid = grid, TT = Tnew[, i], Qinf = Qinf,lrw=lrw)$solution
         index = loc_obsHT[[j]]
         tmp[index]
       })


      # https://hansekbrand.se/code/clusterExport.html
      ## it seems parallel package has some problem for local variable broadcasting.
      #  clusterExport(cl,varlist = c("Fsteady2dsim"))
      #  clusterExport(cl,varlist = c("nsim","grid","Qinf","Tnew",'lrw'))
      #   h_list <- parLapply(cl,1:nsim, function(i){
      #     tmp = Fsteady2dsim(grid = grid, TT = Tnew[, i], Qinf = Qinf,lrw=lrw)$solution
      #     index = loc_obsHT[[j]]
      #     tmp[index]
      #   }
      #   )
      # 一次性合并所有结果
      h_matrix <- do.call(rbind, h_list)
    })

    # Stop the cluster
    # stopCluster(cl)
    # obsh: nsim*nobs
    if(nHT>1) obsh <- do.call("cbind",hHT) else obsh <- hHT[[1]]
  #2. -----------get the obsh......------10.51  0.41 10.97

    # get the head variance for each observation.
  #3. ----------- get the covh and covhk..............
    varobsh <- apply(obsh,2,var)
    meanobsh <- apply(obsh,2,mean)

    # get the misfit.
    weigs <- 1/varobsh/sum(1/varobsh)
    rmse <- mean((trueobsh - meanobsh)^2*weigs)^0.5
    l2 <- mean((trueobsh - meanobsh)^2)^0.5
    l1 <- mean(abs(trueobsh - meanobsh))

    # get the covariance of h and h-T
    # covh <- cov(obsh)
    # data = list(obsh,t(log(Tnew)))
    # # notice it is lnK rather than K.
    # cc <- function(obsh,TT){
    #   covhk = cov(obsh,TT)
    # }
    # covhk = suppressWarnings(multiApply::Apply(data,target_dims = list(1,1),cc))
    #3. ----------- get the covh and covhk.............. blh 37.38  0.05 37.46
    ## we now update the covariance calculation method.
    # obsh --- n_ensemble * n_obs.
    covh = cov(obsh) # the default is y = x
    # covh n_obs*n_obs
    tmpy = t(log(Tnew))
    covhk = cov(obsh,tmpy)
    if (ifcor) {
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }
    # add stablizer term.
    #4. ----------- solve covh ..............
    covh1 <- covh
    #diag(covh1) <-  (1+mul)*diag(covh)
    diag(covh1) <-  rep((1+mul)*max(diag(covh)),nobs)
    #### now we need to invert the covariance function and get the covh-1* %*% covhf
    #a = solve(covh1,covhk[[1]])
    a = solve(covh1,covhk)
    # a in nobs*nelem

    ### now we update the k value based on the difference between obs. and sim.
    for(i in 1:nsim){
      Tnew[,i] <- Tnew[,i] *  exp( t(a) %*% (trueobsh - obsh[i,]))
    }

    msg <- paste('niter =', niter,
                 'varmeanT=', round(varmeanT,4),
                 'rmse=',round(rmse,4),
                 "l2=",round(l2,4),
                 "l1=",round(l1,4))
    msgdf <- rbind(msgdf,c(niter,varmeanT,rmse))
    #4. ----------- solve covh ..............0.42 0.07 0.49
    print(msg)
    ### we need to store the iteration data.
    iterdf[[niter]] <- list(meanT = as.vector(meanT), varT = as.vector(varT), meanobsh = meanobsh, varobsh = varobsh)
    if(niter == 1){
      et = Sys.time()
      print(" -------------------The Computational time for one iteration--------- ")
      print( difftime(et,startTime))
      print(" -------------------The Computational time for one iteration--------- ")
    }
    niter <- niter + 1
    mul <- mul/decay

    ### iteration over time.
  }
  endTime <- Sys.time()
  # get the time in seconds or mins.
  print("-------------------The Computational--------- ")
  print(difftime(endTime,startTime))
  print("-------------------The Computational--------- ")
  return(iterdf)
}


#' 2D Ensemble-based Bayesian inversion for hydraulic tomography (foreach serial)
#'
#' Variant of \code{\link{Finverse}} that uses the \pkg{foreach} framework
#' (\code{\%do\%}) for forward simulations instead of \code{lapply}.  The
#' Bayesian updating algorithm is identical to \code{Finverse}.  See
#' \code{\link{Finverse}} for the full algorithm description and parameter
#' documentation.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)} \code{[L]}.
#' @param grid grid list from \code{GenGrid()}.
#' @param qHT list of pumping test data frames (\code{Qp, x, y}).
#' @param nsim ensemble size (default 50).
#' @param itermax maximum number of iterations.
#' @param varmeanTmax convergence threshold on variance of mean ln(T).
#' @param rmsemin minimum RMSE stopping criterion.
#' @param mul stabiliser coefficient (default 1.0).
#' @param decay stabiliser decay factor (default 1.05).
#' @param oHT list of observation data frames (\code{data, x, y}).
#' @param lrw real work array length for \code{steady.2D} (default 160000).
#' @param geo prior geostatistical parameters: \code{list(me, var, geomod, anis, range, nugget)}.
#' @return A list of per-iteration results (see \code{\link{Finverse}}).
#' @export
#' @examples
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT <- trueK[,-c(1,2)]
#' domain=c(40,40,0,40,0,40)
#' grid = GenGrid(domain)
#' Qinf1=data.frame(Qp=10,xp=20.5,yp=20.5)
#' qHT <- list(test1 = Qinf1)
#' trueh <- Fsteady2dsim(TT=TT,Qinf=Qinf1,domain=domain)
#' locx = c(15,18,22,25,30)
#' locy = c(15,18,22,25,30)
#' loc= expand.grid(x=locx,y=locy)
#' Oinf1 <- data.frame(data=NA,x=loc$x,y=loc$y)
#' Oinf1 <- samData(grid = grid,Oinf = Oinf1,h=trueh$solution)
#' oHT <- list(test1 = Oinf1)
#' result <- Finverse2(grid =grid, qHT = qHT, oHT = oHT)
#' inversePlot(niterm=5,iterdf = result,oHT = oHT,trueK=trueK,grid=grid)
#' Qinf2=data.frame(Qp=10,xp=10.5,yp=10.5)
#' trueh2 <- Fsteady2dsim(TT=TT,Qinf=Qinf2,domain=domain)
#' Oinf2 <- Oinf1
#' Oinf2 <- samData(grid = grid,Oinf = Oinf2,h=trueh2$solution)
#' oHT <- list(test1 = Oinf1,test2 = Oinf2)
#' qHT <- list(test1 = Qinf1,test2 = Qinf2)
#' result <- Finverse2(grid =grid, qHT = qHT, oHT = oHT)
#' inversePlot(niterm=5,iterdf = result,oHT = oHT,trueK=trueK,grid=grid)

Finverse2 <- function(
    domain=c(40,40,0,40,0,40),
    grid = NULL,
    qHT=list(data.frame(Qp=10,x=20.5,y=20.5)), # should be list since it is HT.
    nsim=50,
    itermax=5,
    varmeanTmax =5,
    rmsemin = 0,
    mul=1.0,
    decay = 1.05,
    oHT= list(data.frame(data=-1,x=11,y=11)),
    lrw=160000,
    geo=list(me=0,var=1,geomod="Exp",anis=c(90,1),range=30,nugget=0), # should be list since multiple pumping test.
    ifcor = FALSE)
{
  #1. ----- before iteration (generating ensemble.....)----------
  set.seed(200)
  ### record the time so to see how long it takes.
  startTime <- Sys.time()
  if(is.null(grid))grid = GenGrid(domain)

  #nsim = 50
  nHT <- length(qHT) # number of pumping test.also list length of oHT.
  # create the ij grid.
  library('foreach')
  #isim = 1:nsim
  #jHT = 1: nHT
  #ij = expand.grid(isim = isim, jHT = jHT)
  trueobshHT <- list()
  loc_obsHT <- list()
  for (i in 1:nHT){
    oinf <- oHT[[i]]
    oelemdf <- getOelem(domain = domain,grid = grid,Oinf = oinf)
    loc_obsHT[[i]] <- oelemdf$nelem
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT) # the data format in HT is a list.
  nobs <- length(trueobsh)

  yy <- random2d(nsim=nsim,grid = grid, geo = geo)
  # initial ensemble.
  Tnew <- yy[,-c(1,2)]
  # nelem*nsim
  ### get the variance map of Tnew.
  ### this should goes with iterations.
  niter <- 1
  varmeanT <- 0
  rmse <- 1e10
  #itermax <- 100
  #varmeanTmax <- 5
  #rmsemin <- 0
  #mul <- 1.0 # stablizer.
  msgdf <- data.frame(niter = niter,varmeanT=varmeanT, rmse=rmse)
  iterdf <- list()
  #1. ----- before iteration----------0.94 0.00 0.94 (blh_synthetic)
  # start of the itertion loop.
  while(niter<=itermax & varmeanT<varmeanTmax & rmse>rmsemin){
    #2. -----------get the obsh......------
    varT <- apply(log(Tnew),1,var)
    meanT <- apply(log(Tnew),1,mean)
    varmeanT = var(meanT)

# for each test and each nsim, simulate.

    hHT <- lapply(1:nHT, function(j) {
      Qinf <- qHT[[j]]
      # 使用 lapply 生成所有模拟结果
      h_mat <- foreach(i= 1:nsim,.combine=rbind) %do% {
        tmp = Fsteady2dsim(grid = grid, TT = Tnew[, i], Qinf = Qinf,lrw=lrw)$solution
        index = loc_obsHT[[j]]
        tmp[index]
      }
    })

    # Stop the cluster
    # stopCluster(cl)
    # obsh: nsim*nobs
    if(nHT>1) obsh <- do.call("cbind",hHT) else obsh <- hHT[[1]]
    #2. -----------get the obsh......------10.51  0.41 10.97

    # get the head variance for each observation.
    #3. ----------- get the covh and covhk..............
    varobsh <- apply(obsh,2,var)
    meanobsh <- apply(obsh,2,mean)
    # get the misfit.
    weigs <- 1/varobsh/sum(1/varobsh)
    rmse <- mean((trueobsh - meanobsh)^2*weigs)^0.5
    l2 <- mean((trueobsh - meanobsh)^2)^0.5
    l1 <- mean(abs(trueobsh - meanobsh))

    ## we now update the covariance calculation method.
    # obsh --- n_ensemble * n_obs.
    covh = cov(obsh) # the default is y = x
    # covh n_obs*n_obs
    tmpy = t(log(Tnew))
    covhk = cov(obsh,tmpy)
    if (ifcor) {
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }
    # add stablizer term.
    #4. ----------- solve covh ..............
    covh1 <- covh
    #diag(covh1) <-  (1+mul)*diag(covh)
    diag(covh1) <-  rep((1+mul)*max(diag(covh)),nobs)
    #### now we need to invert the covariance function and get the covh-1* %*% covhf
    #a = solve(covh1,covhk[[1]])
    a = solve(covh1,covhk)
    # a in nobs*nelem

    ### now we update the k value based on the difference between obs. and sim.
    for(i in 1:nsim){
      Tnew[,i] <- Tnew[,i] *  exp( t(a) %*% (trueobsh - obsh[i,]))
    }

    msg <- paste('niter =', niter,
                 'varmeanT=', round(varmeanT,4),
                 'rmse=',round(rmse,4),
                 "l2=",round(l2,4),
                 "l1=",round(l1,4))
    msgdf <- rbind(msgdf,c(niter,varmeanT,rmse))
    #4. ----------- solve covh ..............0.42 0.07 0.49
    print(msg)
    ### we need to store the iteration data.
    iterdf[[niter]] <- list(meanT = as.vector(meanT), varT = as.vector(varT), meanobsh = meanobsh, varobsh = varobsh)
    if(niter == 1){
      et = Sys.time()
      print(" -------------------The Computational time for one iteration--------- ")
      print( difftime(et,startTime))
      print(" -------------------The Computational time for one iteration--------- ")
    }
    niter <- niter + 1
    mul <- mul/decay

    ### iteration over time.
  }
  endTime <- Sys.time()
  # get the time in seconds or mins.
  print("-------------------The Computational--------- ")
  print(difftime(endTime,startTime))
  print("-------------------The Computational--------- ")
  return(iterdf)
}





#' 2D Ensemble-based Bayesian inversion for hydraulic tomography (parallel)
#'
#' Variant of \code{\link{Finverse}} that uses the \pkg{foreach} +
#' \pkg{doParallel} framework (\code{\%dopar\%}) for parallel forward
#' simulations across ensemble members.  The Bayesian updating algorithm is
#' identical to \code{Finverse}.  See \code{\link{Finverse}} for the full
#' algorithm description.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)} \code{[L]}.
#' @param grid grid list from \code{GenGrid()}.
#' @param qHT list of pumping test data frames (\code{Qp, x, y}).
#' @param nsim ensemble size (default 50).
#' @param itermax maximum number of iterations.
#' @param varmeanTmax convergence threshold on variance of mean ln(T).
#' @param rmsemin minimum RMSE stopping criterion.
#' @param mul stabiliser coefficient (default 1.0).
#' @param decay stabiliser decay factor (default 1.05).
#' @param oHT list of observation data frames (\code{data, x, y}).
#' @param lrw real work array length for \code{steady.2D} (default 160000).
#' @param ncore integer, the number of CPU cores for parallel execution
#'   (default 10).
#' @param geo prior geostatistical parameters: \code{list(me, var, geomod, anis, range, nugget)}.
#' @return A list of per-iteration results (see \code{\link{Finverse}}).
#' @export
#' @examples
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT <- trueK[,-c(1,2)]
#' domain=c(40,40,0,40,0,40)
#' grid = GenGrid(domain)
#' Qinf1=data.frame(Qp=10,xp=20.5,yp=20.5)
#' qHT <- list(test1 = Qinf1)
#' trueh <- Fsteady2dsim(TT=TT,Qinf=Qinf1,domain=domain)
#' locx = c(15,18,22,25,30)
#' locy = c(15,18,22,25,30)
#' loc= expand.grid(x=locx,y=locy)
#' Oinf1 <- data.frame(data=NA,x=loc$x,y=loc$y)
#' Oinf1 <- samData(grid = grid,Oinf = Oinf1,h=trueh$solution)
#' oHT <- list(test1 = Oinf1)
#' result <- Finverse3(grid =grid, qHT = qHT, oHT = oHT)
#' inversePlot(niterm=5,iterdf = result,oHT = oHT,trueK=trueK,grid=grid)
#' Qinf2=data.frame(Qp=10,xp=10.5,yp=10.5)
#' trueh2 <- Fsteady2dsim(TT=TT,Qinf=Qinf2,domain=domain)
#' Oinf2 <- Oinf1
#' Oinf2 <- samData(grid = grid,Oinf = Oinf2,h=trueh2$solution)
#' oHT <- list(test1 = Oinf1,test2 = Oinf2)
#' qHT <- list(test1 = Qinf1,test2 = Qinf2)
#' result <- Finverse3(grid =grid, qHT = qHT, oHT = oHT)
#' inversePlot(niterm=5,iterdf = result,oHT = oHT,trueK=trueK,grid=grid)

Finverse3 <- function(
    domain=c(40,40,0,40,0,40),
    grid = NULL,
    qHT=list(data.frame(Qp=10,x=20.5,y=20.5)), # should be list since it is HT.
    nsim=50,
    itermax=5,
    varmeanTmax =5,
    rmsemin = 0,
    mul=1.0,
    decay = 1.05,
    oHT= list(data.frame(data=-1,x=11,y=11)),
    lrw=160000,
    ncore =10,
    geo=list(me=0,var=1,geomod="Exp",anis=c(90,1),range=30,nugget=0), # should be list since multiple pumping test.
    ifcor = FALSE)
{
  #1. ----- before iteration (generating ensemble.....)----------
  set.seed(200)
  ### record the time so to see how long it takes.
  startTime <- Sys.time()
  if(is.null(grid))grid = GenGrid(domain)

  #nsim = 50
  nHT <- length(qHT) # number of pumping test.also list length of oHT.
  # create the ij grid.
  library('foreach')
  library("doParallel")
  #isim = 1:nsim
  #jHT = 1: nHT
  #ij = expand.grid(isim = isim, jHT = jHT)
  trueobshHT <- list()
  loc_obsHT <- list()
  for (i in 1:nHT){
    oinf <- oHT[[i]]
    oelemdf <- getOelem(domain = domain,grid = grid,Oinf = oinf)
    loc_obsHT[[i]] <- oelemdf$nelem
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT) # the data format in HT is a list.
  nobs <- length(trueobsh)

  yy <- random2d(nsim=nsim,grid = grid, geo = geo)
  # initial ensemble.
  Tnew <- yy[,-c(1,2)]
  # nelem*nsim
  ### get the variance map of Tnew.
  ### this should goes with iterations.
  niter <- 1
  varmeanT <- 0
  rmse <- 1e10
  #itermax <- 100
  #varmeanTmax <- 5
  #rmsemin <- 0
  #mul <- 1.0 # stablizer.
  msgdf <- data.frame(niter = niter,varmeanT=varmeanT, rmse=rmse)
  iterdf <- list()
  #1. ----- before iteration----------0.94 0.00 0.94 (blh_synthetic)
  # start of the itertion loop.
  while(niter<=itermax & varmeanT<varmeanTmax & rmse>rmsemin){
    #2. -----------get the obsh......------
    varT <- apply(log(Tnew),1,var)
    meanT <- apply(log(Tnew),1,mean)
    varmeanT = var(meanT)

    # for each test and each nsim, simulate.
    obsh = samHTmcPar(grid =grid, TT=Tnew, qHT= qHT, oHT=oHT,lrw=lrw,ncore =ncore)

    #2. -----------get the obsh......------10.51  0.41 10.97

    # get the head variance for each observation.
    #3. ----------- get the covh and covhk..............
    varobsh <- apply(obsh,2,var)
    meanobsh <- apply(obsh,2,mean)
    # get the misfit.
    weigs <- 1/varobsh/sum(1/varobsh)
    rmse <- mean((trueobsh - meanobsh)^2*weigs)^0.5
    l2 <- mean((trueobsh - meanobsh)^2)^0.5
    l1 <- mean(abs(trueobsh - meanobsh))

    ## we now update the covariance calculation method.
    # obsh --- n_ensemble * n_obs.
    covh = cov(obsh) # the default is y = x
    # covh n_obs*n_obs
    tmpy = t(log(Tnew))
    covhk = cov(obsh,tmpy)
    if (ifcor) {
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }
    # add stablizer term.
    #4. ----------- solve covh ..............
    covh1 <- covh
    #diag(covh1) <-  (1+mul)*diag(covh)
    diag(covh1) <-  rep((1+mul)*max(diag(covh)),nobs)
    #### now we need to invert the covariance function and get the covh-1* %*% covhf
    #a = solve(covh1,covhk[[1]])
    a = solve(covh1,covhk)
    # a in nobs*nelem

    ### now we update the k value based on the difference between obs. and sim.
    for(i in 1:nsim){
      Tnew[,i] <- Tnew[,i] *  exp( t(a) %*% (trueobsh - obsh[i,]))
    }

    msg <- paste('niter =', niter,
                 'varmeanT=', round(varmeanT,4),
                 'rmse=',round(rmse,4),
                 "l2=",round(l2,4),
                 "l1=",round(l1,4))
    msgdf <- rbind(msgdf,c(niter,varmeanT,rmse))
    #4. ----------- solve covh ..............0.42 0.07 0.49
    print(msg)
    ### we need to store the iteration data.
    iterdf[[niter]] <- list(meanT = as.vector(meanT), varT = as.vector(varT), meanobsh = meanobsh, varobsh = varobsh)
    if(niter == 1){
      et = Sys.time()
      print(" -------------------The Computational time for one iteration--------- ")
      print( difftime(et,startTime))
      print(" -------------------The Computational time for one iteration--------- ")
    }
    niter <- niter + 1
    mul <- mul/decay

    ### iteration over time.
  }
  endTime <- Sys.time()
  # get the time in seconds or mins.
  print("-------------------The Computational--------- ")
  print(difftime(endTime,startTime))
  print("-------------------The Computational--------- ")
  return(iterdf)
}




#' Run forward simulations for one T field across all 2D pumping tests
#'
#' For a single transmissivity field, runs \code{Fsteady2dsim} for each pumping
#' test in \code{qHT} and samples the simulated heads at the observation
#' locations given in \code{oHT}.  Used internally by the inversion functions.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)} \code{[L]}.
#' @param grid grid list from \code{GenGrid()}.
#' @param TT transmissivity \eqn{T} \code{[L\eqn{^2}/T]} — either a scalar or a
#'   length-\code{n} vector (one value per cell).
#' @param qHT list of pumping test data frames, each with columns \code{Qp}
#'   \code{[L\eqn{^3}/T]}, \code{x} \code{[L]}, \code{y} \code{[L]}.
#' @param oHT list of observation data frames, each with columns \code{data}
#'   \code{[L]}, \code{x} \code{[L]}, \code{y} \code{[L]}.
#' @param lrw real work array length for \code{steady.2D} (default 160000).
#' @param simplify if \code{TRUE} (default) return a plain vector of sampled
#'   heads concatenated across all tests; if \code{FALSE} return the updated
#'   \code{oHT} list.
#' @return A numeric vector of sampled drawdown values (when
#'   \code{simplify = TRUE}), or the updated \code{oHT} list (when
#'   \code{simplify = FALSE}).
#' @export
#' @examples
#' domain=c(80,80,0,80,0,80)
#' grid <- GenGrid(domain=domain)
#' TT =0.1
#' data("oHT")
#' data("qHT")
#' da = samHT(grid =grid, TT=TT, qHT= qHT, oHT=oHT,lrw=320000)
samHT <- function(domain=c(40,40,0,40,0,40),
                   grid = NULL,
                   TT=0.1,
                   qHT=list(data.frame(Qp=10,x=20.5,y=20.5)),
                   oHT= list(data.frame(data=-1,x=11,y=11)),
                   lrw = 160000,
                   simplify =TRUE
 ){
  if(is.null(grid))grid = GenGrid(domain)
  nHT <- length(qHT) # number of pumping test.also list length of oHT.
  for (i in 1:nHT){
    qinf <- qHT[[i]]
    oinf <- oHT[[i]]
    trueh <- Fsteady2dsim(TT=TT,Qinf=qinf,grid=grid,lrw=lrw)$solution
    oinf = samData(grid = grid,Oinf = oinf,h=trueh)
    oHT[[i]] <- oinf[,c('data','x','y')]
  }
  if(simplify){
    oHTdf = dplyr::bind_rows(oHT,.id='id')
    return(oHTdf$data)
  }else {return(oHT)}
 }

#' Serial Monte Carlo forward runs for 2D HT (ensemble)
#'
#' Runs \code{samHT} for each ensemble member using \pkg{foreach}
#' (\code{\%do\%}, serial).  Returns a matrix where each row corresponds to
#' one ensemble member and columns are observed heads concatenated across all
#' pumping tests.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)} \code{[L]}.
#' @param grid grid list from \code{GenGrid()}.
#' @param TT \code{n × nsim} matrix of transmissivity realisations \code{[L\eqn{^2}/T]}.
#' @param qHT list of pumping test data frames (\code{Qp, x, y}).
#' @param oHT list of observation data frames (\code{data, x, y}).
#' @param lrw real work array length for \code{steady.2D} (default 160000).
#' @return An \code{nsim × nobs} matrix of simulated heads at observation wells.
#' @export
#' @examples
#' domain=c(80,80,0,80,0,80)
#' grid <- GenGrid(domain=domain)
#' TT= random2d(nsim=10,grid=grid)
#' TT = TT[,-c(1,2)]
#' data("oHT")
#' data("qHT")
#' da = samHTmc(grid =grid, TT=TT, qHT= qHT, oHT=oHT,lrw=320000)
samHTmc <- function(domain=c(40,40,0,40,0,40),
                  grid = NULL,
                  TT=0.1,
                  qHT=list(data.frame(Qp=10,x=20.5,y=20.5)),
                  oHT= list(data.frame(data=-1,x=11,y=11)),
                  lrw = 160000

){
  nsim = ncol(TT)
  library('foreach')
  foreach(i =1:nsim,.combine = rbind) %do%{
    samHT(grid =grid, TT=TT[,i], qHT= qHT, oHT=oHT,lrw=lrw)
  }
}


#' Parallel Monte Carlo forward runs for 2D HT (ensemble)
#'
#' Runs \code{samHT} for each ensemble member in parallel using
#' \pkg{foreach} + \pkg{doParallel} (\code{\%dopar\%}).  Returns a matrix
#' where each row corresponds to one ensemble member and columns are observed
#' heads concatenated across all pumping tests.  Used internally by
#' \code{\link{Finverse3}}.
#'
#' @param domain a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)} \code{[L]}.
#' @param grid grid list from \code{GenGrid()}.
#' @param TT \code{n × nsim} matrix of transmissivity realisations \code{[L\eqn{^2}/T]}.
#' @param qHT list of pumping test data frames (\code{Qp, x, y}).
#' @param oHT list of observation data frames (\code{data, x, y}).
#' @param lrw real work array length for \code{steady.2D} (default 160000).
#' @param ncore number of CPU cores for parallel execution (default 7).
#' @return An \code{nsim × nobs} matrix of simulated heads at observation wells.
#' @export
#' @examples
#' domain=c(80,80,0,80,0,80)
#' grid <- GenGrid(domain=domain)
#' TT= random2d(nsim=50,grid=grid)
#' TT = TT[,-c(1,2)]
#' data("oHT")
#' data("qHT")
#' da = samHTmcPar(grid =grid, TT=TT, qHT= qHT, oHT=oHT,lrw=320000,ncore=10)
samHTmcPar <- function(domain=c(40,40,0,40,0,40),
                    grid = NULL,
                    TT=0.1,
                    qHT=list(data.frame(Qp=10,x=20.5,y=20.5)),
                    oHT= list(data.frame(data=-1,x=11,y=11)),
                    lrw = 160000,
                    ncore = 7

){
  library('foreach')
  library("doParallel")
  registerDoParallel(cores=ncore)
  nsim = ncol(TT)
  x <-foreach(i =1:nsim,.combine = rbind) %dopar%{
    HydroTomo::samHT(grid =grid, TT=TT[,i], qHT= qHT, oHT=oHT,lrw=lrw)
  }


  stopImplicitCluster()
  return(x)
}


# ==============================================================================
# 3D Hydraulic Tomography — Inverse Functions
# ==============================================================================

#' Run forward simulations for one K field across all 3D pumping tests
#'
#' 3D counterpart of \code{samHT}. Calls \code{Fsteady3dsim} for each pumping
#' test and samples simulated heads at observation well locations.
#'
#' @param grid  grid from \code{GenGrid3D()}.
#' @param TT    length-\code{n} vector of hydraulic conductivity K \code{[L/T]}.
#' @param qHT   list of pumping test data frames, each with columns
#'   \code{Qp, x, y, z}.
#' @param oHT   list of observation data frames, each with columns
#'   \code{data, x, y, z}.
#' @param lrw   real work array length for \code{steady.3D}; default 20000000.
#' @param simplify if \code{TRUE} return a plain vector of sampled heads
#'   (all tests concatenated); if \code{FALSE} return updated \code{oHT}.
#' @return vector of sampled drawdown values (when \code{simplify = TRUE}).
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15,15,5,0,15,0,15,0,5))
#' Qinf   <- data.frame(Qp=10, x=7.5, y=7.5, z=2.5)
#' Oinf   <- data.frame(data=NA, x=c(10,12), y=c(7.5,7.5), z=c(2.5,2.5))
#' da     <- samHT3D(grid=grid3d, TT=0.1,
#'                   qHT=list(Qinf), oHT=list(Oinf))
samHT3D <- function(grid,
                    TT      = 0.1,
                    qHT     = list(data.frame(Qp=10, x=10.5, y=10.5, z=2.5)),
                    oHT     = list(data.frame(data=NA, x=5, y=5, z=2.5)),
                    lrw     = 20000000,
                    simplify = TRUE) {

  nHT <- length(qHT)
  for (i in seq_len(nHT)) {
    qinf     <- qHT[[i]]
    oinf     <- oHT[[i]]
    h        <- Fsteady3dsim(KK = TT, Qinf = qinf, grid = grid, lrw = lrw)$solution
    oinf     <- samData3D(Oinf = oinf, grid = grid, h = h)
    oHT[[i]] <- oinf[, c('data', 'x', 'y', 'z')]
  }
  if (simplify) {
    oHTdf <- dplyr::bind_rows(oHT, .id = 'id')
    return(oHTdf$data)
  } else {
    return(oHT)
  }
}


#' Parallel Monte Carlo forward runs for 3D HT (ensemble)
#'
#' 3D counterpart of \code{samHTmcPar}. Runs \code{samHT3D} in parallel across
#' all ensemble members.
#'
#' @param grid   grid from \code{GenGrid3D()}.
#' @param TT     \code{n × nsim} matrix of K realisations.
#' @param qHT    list of pumping test data frames (columns \code{Qp, x, y, z}).
#' @param oHT    list of observation data frames (columns \code{data, x, y, z}).
#' @param lrw    real work array length; default 20000000.
#' @param ncore  number of parallel cores.
#' @return \code{nsim × nobs} matrix of simulated heads at observation wells.
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15,15,5,0,15,0,15,0,5))
#' KK     <- random3d(nsim=5, grid=grid3d)
#' TT     <- as.matrix(KK[,-c(1,2,3)])
#' Qinf   <- data.frame(Qp=10, x=7.5, y=7.5, z=2.5)
#' Oinf   <- data.frame(data=NA, x=c(10,12), y=c(7.5,7.5), z=c(2.5,2.5))
#' da     <- samHTmcPar3D(grid=grid3d, TT=TT,
#'                        qHT=list(Qinf), oHT=list(Oinf), ncore=2)
samHTmcPar3D <- function(grid,
                         TT    = 0.1,
                         qHT   = list(data.frame(Qp=10, x=10.5, y=10.5, z=2.5)),
                         oHT   = list(data.frame(data=NA, x=5, y=5, z=2.5)),
                         lrw   = 20000000,
                         ncore = 4) {

  library('foreach')
  library('doParallel')
  registerDoParallel(cores = ncore)

  nsim <- ncol(TT)
  x <- foreach(i = 1:nsim, .combine = rbind) %dopar% {
    HydroTomo::samHT3D(grid = grid, TT = TT[, i],
                       qHT = qHT, oHT = oHT, lrw = lrw)
  }

  stopImplicitCluster()
  return(x)
}


# ==============================================================================
# 3D Steady-State HT — Well-Screen (Vertical Interval) Functions
# ==============================================================================

#' Run forward simulations for one K field — well-screen version (3D steady)
#'
#' Steady-state counterpart of \code{\link{samHT3DtrScreen}}.  Calls
#' \code{Fsteady3dsim} for each pumping test and averages simulated heads
#' over the vertical screen interval of each observation well.
#'
#' @param grid  grid from \code{GenGrid3D()}.
#' @param TT    length-\code{n} vector of hydraulic conductivity K \code{[L/T]}.
#' @param qHT   list of pumping test data frames, each with columns
#'   \code{Qp, x, y, z_top, z_bottom} (well-screen pumping) or
#'   \code{Qp, x, y, z} (point pumping).
#' @param oHT   list of observation data frames, each with columns
#'   \code{data, x, y, z_top, z_bottom}.
#' @param lrw   real work array length for \code{steady.3D}; default 20000000.
#' @param simplify if \code{TRUE} return a plain vector of sampled heads
#'   (all tests concatenated); if \code{FALSE} return updated \code{oHT}.
#' @return vector of interval-averaged drawdown values (when \code{simplify = TRUE}).
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15,15,5,0,15,0,15,0,5))
#' Qinf   <- data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)
#' Oinf   <- data.frame(data=NA, x=c(5,10), y=c(7.5,7.5),
#'                      z_top=c(1,1), z_bottom=c(4,4))
#' da     <- samHT3DScreen(grid=grid3d, TT=0.1,
#'                         qHT=list(Qinf), oHT=list(Oinf))
samHT3DScreen <- function(grid,
                          TT       = 0.1,
                          qHT      = list(data.frame(Qp=10, x=7.5, y=7.5,
                                                     z_top=2, z_bottom=3)),
                          oHT      = list(data.frame(data=NA, x=5, y=5,
                                                     z_top=1, z_bottom=4)),
                          lrw      = 20000000,
                          simplify = TRUE) {

  nHT <- length(qHT)
  for (i in seq_len(nHT)) {
    qinf <- qHT[[i]]
    oinf <- oHT[[i]]
    h    <- Fsteady3dsim(KK = TT, Qinf = qinf, grid = grid, lrw = lrw)$solution
    oinf <- samData3DScreen(Oinf = oinf, grid = grid, h = h)
    oHT[[i]] <- oinf[, c('data', 'x', 'y', 'z_top', 'z_bottom')]
  }
  if (simplify) {
    oHTdf <- dplyr::bind_rows(oHT, .id = 'id')
    return(oHTdf$data)
  } else {
    return(oHT)
  }
}


#' Parallel MC forward runs for 3D steady HT — well-screen version (ensemble)
#'
#' Like \code{\link{samHTmcPar3D}} but uses \code{\link{samHT3DScreen}}
#' for interval-averaged well-screen sampling.
#'
#' @param grid   grid from \code{GenGrid3D()}.
#' @param TT     \code{n × nsim} matrix of K realisations.
#' @param qHT    list of pumping test data frames (columns
#'   \code{Qp, x, y, z_top, z_bottom} for well-screen pumping, or
#'   \code{Qp, x, y, z} for point pumping).
#' @param oHT    list of observation data frames (columns
#'   \code{data, x, y, z_top, z_bottom}).
#' @param lrw    real work array length; default 20000000.
#' @param ncore  number of parallel cores.
#' @return \code{nsim × nobs} matrix of interval-averaged heads.
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15,15,5,0,15,0,15,0,5))
#' KKmat  <- random3d(nsim=5, grid=grid3d)
#' KKmat  <- as.matrix(KKmat[,-c(1,2,3)])
#' Qinf   <- data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)
#' Oinf   <- data.frame(data=NA, x=c(5,10), y=c(7.5,7.5),
#'                      z_top=c(1,1), z_bottom=c(4,4))
#' da     <- samHTmcPar3DScreen(grid=grid3d, TT=KKmat,
#'                              qHT=list(Qinf), oHT=list(Oinf), ncore=2)
samHTmcPar3DScreen <- function(grid,
                               TT    = 0.1,
                               qHT   = list(data.frame(Qp=10, x=7.5, y=7.5,
                                                       z_top=2, z_bottom=3)),
                               oHT   = list(data.frame(data=NA, x=5, y=5,
                                                       z_top=1, z_bottom=4)),
                               lrw   = 20000000,
                               ncore = 4) {

  library('foreach')
  library('doParallel')
  registerDoParallel(cores = ncore)

  nsim <- ncol(TT)
  x <- foreach(i = 1:nsim, .combine = rbind) %dopar% {
    HydroTomo::samHT3DScreen(grid = grid, TT = TT[, i],
                             qHT = qHT, oHT = oHT, lrw = lrw)
  }

  stopImplicitCluster()
  return(x)
}


#' 3D Ensemble-based Bayesian inverse for hydraulic tomography (parallel)
#'
#' 3D extension of \code{Finverse3}, replacing all 2D components (grid,
#' random field, forward solver, observation sampling) with their 3D
#' counterparts.  The Bayesian updating equations are identical.
#'
#' @param domain  9-element vector \code{c(nx,ny,nz, x1,x2, y1,y2, z1,z2)}.
#'   **Keep small for testing** (e.g. \code{c(15,15,5,0,15,0,15,0,5)}).
#' @param grid    grid from \code{GenGrid3D()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT     list of pumping test data frames with columns
#'   \code{Qp, x, y, z_top, z_bottom} (well-screen pumping) or
#'   \code{Qp, x, y, z} (point pumping).
#' @param nsim    ensemble size (default 50; use ~20 for quick tests).
#' @param itermax maximum iterations (set to 1 for initial testing).
#' @param varmeanTmax  convergence threshold on variance of mean ln(K).
#' @param rmsemin minimum RMSE stopping criterion.
#' @param mul     stabiliser (default 1.0).
#' @param decay   stabiliser decay per iteration (default 1.05).
#' @param oHT     list of observation data frames with columns
#'   \code{data, x, y, z}.
#' @param lrw     real work array length for \code{steady.3D}
#'   (default 20000000; increase for larger grids).
#' @param ncore   number of parallel cores (default 4).
#' @param geo     prior geostatistical parameters — same list structure as
#'   \code{random3d()}: \code{list(me, var, geomod, range, nugget, anis)}.
#'   Note \code{anis} is a 5-element vector for 3D.
#' @param ifcor    logical; if \code{TRUE}, returns the ensemble results after
#'   the first iteration without performing the update. Useful for checking
#'   forward simulations before running the full inversion.
#' @return list of per-iteration results, each a list with
#'   \code{meanT} (mean ln K vector), \code{varT} (variance of ln K),
#'   \code{meanobsh}, \code{varobsh}.  When \code{ifcor = TRUE}, returns
#'   list with \code{obsh} (simulated heads), \code{Tnew} (ensemble), and
#'   \code{meanobsh}.
#' @export
#' @examples
#' # --- small synthetic test (15x15x5, 1 iteration) ---
#' domain3d <- c(15, 15, 5, 0, 15, 0, 15, 0, 5)
#' grid3d   <- GenGrid3D(domain3d)
#' set.seed(42)
#' trueK3d  <- random3d(nsim=1, grid=grid3d)
#' Qinf3d   <- data.frame(Qp=10, x=7.5, y=7.5, z=2.5)
#' qHT3d    <- list(test1 = Qinf3d)
#' trueh3d  <- Fsteady3dsim(grid=grid3d, KK=trueK3d$Kp, Qinf=Qinf3d)
#' loc      <- expand.grid(x=c(3,6,9,12), y=c(3,6,9,12))
#' Oinf3d   <- data.frame(data=NA, x=loc$x, y=loc$y, z=2.5)
#' Oinf3d   <- samData3D(Oinf=Oinf3d, grid=grid3d, h=trueh3d$solution)
#' oHT3d    <- list(test1 = Oinf3d)
#' result3d <- Finverse3D(grid=grid3d, qHT=qHT3d, oHT=oHT3d,
#'                        nsim=20, itermax=1, ncore=2)
Finverse3D <- function(
    domain      = c(15, 15, 5, 0, 15, 0, 15, 0, 5),
    grid        = NULL,
    qHT         = list(data.frame(Qp=10, x=7.5, y=7.5, z=2.5)),
    nsim        = 50,
    itermax     = 5,
    varmeanTmax = 5,
    rmsemin     = 0,
    mul         = 1.0,
    decay       = 1.05,
    oHT         = list(data.frame(data=-1, x=5, y=5, z=2.5)),
    lrw         = 20000000,
    ncore       = 4,
    geo         = list(me=0, var=1, geomod="Exp",
                       range=10, nugget=0,
                       anis=c(0, 0, 0, 1, 1)),
    ifcor       = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid3D(domain)

  nHT <- length(qHT)

  # ---- extract observation locations and true data ---------------------------
  trueobshHT <- list()
  loc_obsHT  <- list()
  for (i in seq_len(nHT)) {
    oinf           <- oHT[[i]]
    oelemdf        <- getOelem3D(grid = grid, Oinf = oinf)
    loc_obsHT[[i]] <- oelemdf$nelem
    trueobshHT[[i]]<- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs     <- length(trueobsh)

  # ---- initial ensemble ------------------------------------------------------
  yy   <- random3d(nsim = nsim, grid = grid, geo = geo)
  Tnew <- as.matrix(yy[, -c(1, 2, 3)])   # n × nsim matrix of K values

  # ---- iteration loop --------------------------------------------------------
  niter    <- 1
  varmeanT <- 0
  rmse     <- 1e10
  msgdf    <- data.frame(niter = niter, varmeanT = varmeanT, rmse = rmse)
  iterdf   <- list()

  while (niter <= itermax & varmeanT < varmeanTmax & rmse > rmsemin) {

    varT     <- apply(log(Tnew), 1, var)
    meanT    <- apply(log(Tnew), 1, mean)
    varmeanT <- var(meanT)

    # parallel forward runs for all ensemble members and pumping tests
    obsh <- samHTmcPar3D(grid  = grid,
                         TT    = Tnew,
                         qHT   = qHT,
                         oHT   = oHT,
                         lrw   = lrw,
                         ncore = ncore)

    # ---- statistics ----------------------------------------------------------
    varobsh  <- apply(obsh, 2, var)
    meanobsh <- apply(obsh, 2, mean)
    weigs    <- 1 / varobsh / sum(1 / varobsh)
    rmse     <- mean((trueobsh - meanobsh)^2 * weigs)^0.5
    l2       <- mean((trueobsh - meanobsh)^2)^0.5
    l1       <- mean(abs(trueobsh - meanobsh))

    # ---- Bayesian update -----------------------------------------------------
    covh  <- cov(obsh)
    covhk <- cov(obsh, t(log(Tnew)))
    if (ifcor) {
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }
    covh1 <- covh
    diag(covh1) <- rep((1 + mul) * max(diag(covh)), nobs)
    a <- solve(covh1, covhk)   # nobs × n

    for (i in seq_len(nsim)) {
      Tnew[, i] <- Tnew[, i] * exp(t(a) %*% (trueobsh - obsh[i, ]))
    }

    msg <- paste('niter =', niter,
                 'varmeanT =', round(varmeanT, 4),
                 'rmse =', round(rmse, 4),
                 'l2 =', round(l2, 4),
                 'l1 =', round(l1, 4))
    print(msg)

    iterdf[[niter]] <- list(meanT    = as.vector(meanT),
                            varT     = as.vector(varT),
                            meanobsh = meanobsh,
                            varobsh  = varobsh)

    if (niter == 1) {
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    mul   <- mul / decay
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}


# ==============================================================================
# 3D Steady-State HT — Well-Screen Inverse
# ==============================================================================

#' 3D Steady HT ensemble inverse — well-screen (vertical interval) version
#'
#' Like \code{\link{Finverse3D}} but for wells with vertical screens
#' (filter sections).  Observation wells are defined by \code{(x, y)} and a
#' screen interval \code{[z_top, z_bottom]}.  The forward solver
#' \code{\link{Fsteady3dsim}} computes full 3D steady heads, and
#' \code{\link{samData3DScreen}} averages heads over all grid cells within
#' the screen interval.
#'
#' @param domain  9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2, z1, z2)}.
#' @param grid    grid from \code{GenGrid3D()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT     list of pumping test data frames with columns
#'   \code{Qp, x, y, z_top, z_bottom} (well-screen pumping) or
#'   \code{Qp, x, y, z} (point pumping).
#' @param nsim    ensemble size (default 50).
#' @param itermax maximum iterations.
#' @param varmeanTmax  convergence threshold on variance of mean ln(K).
#' @param rmsemin minimum RMSE stopping criterion.
#' @param mul     stabiliser (default 1.0).
#' @param decay   stabiliser decay per iteration (default 1.05).
#' @param oHT     list of observation data frames with columns
#'   \code{data, x, y, z_top, z_bottom}.
#' @param lrw     real work array length for \code{steady.3D} (default 20000000).
#' @param ncore   number of parallel cores (default 4).
#' @param geo     prior geostatistical parameters (list).
#' @param ifcor    logical; if \code{TRUE}, returns the ensemble results after
#'   the first iteration without performing the update.
#' @return list of per-iteration results.  When \code{ifcor = TRUE}, returns
#'   list with \code{obsh}, \code{Tnew}, and \code{meanobsh}.
#' @export
#' @examples
#' # --- small synthetic test with well screens ---
#' domain3d <- c(15, 15, 5, 0, 15, 0, 15, 0, 5)
#' grid3d   <- GenGrid3D(domain3d)
#' set.seed(42)
#' trueK3d  <- random3d(nsim=1, grid=grid3d)
#' Qinf3d   <- data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)
#' qHT3d    <- list(test1 = Qinf3d)
#' res3d    <- Fsteady3dsim(grid=grid3d, KK=trueK3d$Kp, Qinf=Qinf3d)
#' loc      <- expand.grid(x=c(3,6,9,12), y=c(3,6,9,12))
#' Oinf3d   <- data.frame(data=NA, x=loc$x, y=loc$y,
#'                        z_top=1, z_bottom=4)
#' Oinf3d   <- samData3DScreen(Oinf=Oinf3d, grid=grid3d, h=res3d$solution)
#' oHT3d    <- list(test1 = Oinf3d)
#' result3d <- Finverse3DScreen(grid=grid3d, qHT=qHT3d, oHT=oHT3d,
#'                              nsim=20, itermax=1, ncore=2)
Finverse3DScreen <- function(
    domain      = c(15, 15, 5, 0, 15, 0, 15, 0, 5),
    grid        = NULL,
    qHT         = list(data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)),
    nsim        = 50,
    itermax     = 5,
    varmeanTmax = 5,
    rmsemin     = 0,
    mul         = 1.0,
    decay       = 1.05,
    oHT         = list(data.frame(data=-1, x=5, y=5,
                                 z_top=1, z_bottom=4)),
    lrw         = 20000000,
    ncore       = 4,
    geo         = list(me=0, var=1, geomod="Exp",
                       range=10, nugget=0,
                       anis=c(0, 0, 0, 1, 1)),
    ifcor       = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid3D(domain)

  nHT <- length(qHT)

  # ---- extract observation data ----------------------------------------------
  trueobshHT <- list()
  for (i in seq_len(nHT)) {
    oinf            <- oHT[[i]]
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs     <- length(trueobsh)

  # ---- initial ensemble ------------------------------------------------------
  yy   <- random3d(nsim = nsim, grid = grid, geo = geo)
  Knew <- as.matrix(yy[, -c(1, 2, 3)])

  # ---- iteration loop --------------------------------------------------------
  niter    <- 1
  varmeanT <- 0
  rmse     <- 1e10
  msgdf    <- data.frame(niter = niter, varmeanT = varmeanT, rmse = rmse)
  iterdf   <- list()

  while (niter <= itermax & varmeanT < varmeanTmax & rmse > rmsemin) {

    varT     <- apply(log(Knew), 1, var)
    meanT    <- apply(log(Knew), 1, mean)
    varmeanT <- var(meanT)

    # parallel steady forward runs — screen-averaged
    obsh <- samHTmcPar3DScreen(grid  = grid,
                               TT    = Knew,
                               qHT   = qHT,
                               oHT   = oHT,
                               lrw   = lrw,
                               ncore = ncore)

    # ---- statistics ----------------------------------------------------------
    varobsh  <- apply(obsh, 2, var)
    meanobsh <- apply(obsh, 2, mean)
    weigs    <- 1 / varobsh / sum(1 / varobsh)
    rmse     <- mean((trueobsh - meanobsh)^2 * weigs)^0.5
    l2       <- mean((trueobsh - meanobsh)^2)^0.5
    l1       <- mean(abs(trueobsh - meanobsh))

    # ---- Bayesian update -----------------------------------------------------
    covh  <- cov(obsh)
    covhk <- cov(obsh, t(log(Knew)))
    if (ifcor) {
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }
    covh1 <- covh
    diag(covh1) <- rep((1 + mul) * max(diag(covh)), nobs)
    a <- solve(covh1, covhk)

    for (i in seq_len(nsim)) {
      Knew[, i] <- Knew[, i] * exp(t(a) %*% (trueobsh - obsh[i, ]))
    }

    msg <- paste('niter =', niter,
                 'varmeanT =', round(varmeanT, 4),
                 'rmse =', round(rmse, 4),
                 'l2 =', round(l2, 4),
                 'l1 =', round(l1, 4))
    print(msg)

    iterdf[[niter]] <- list(meanT    = as.vector(meanT),
                            varT     = as.vector(varT),
                            meanobsh = meanobsh,
                            varobsh  = varobsh)

    if (niter == 1) {
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    mul   <- mul / decay
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}


# ==============================================================================
# 3D Adjoint-state sensitivity computation
# ==============================================================================

#' Solve the adjoint equation for a single observation point (3D steady-state)
#'
#' 3D counterpart of \code{\link{adjoint2D}}.  For the 3D steady-state
#' groundwater flow equation \eqn{\nabla \cdot (K \nabla h) = Q}, the adjoint
#' state \eqn{\lambda} for an observation at cell \code{iobs} satisfies:
#' \deqn{\nabla \cdot (K \nabla \lambda) = \delta(\mathbf{x} - \mathbf{x}_{iobs})}
#' with homogeneous (zero) boundary conditions.  The adjoint variable is used
#' to compute the Jacobian (sensitivity) matrix via:
#' \deqn{\frac{\partial h_{iobs}}{\partial (\ln K)_k} =
#'       -K_k \, \nabla h|_k \cdot \nabla \lambda|_k \, \Delta x \, \Delta y \, \Delta z}
#'
#' @param grid    grid list from \code{GenGrid3D()}.
#' @param KK      hydraulic conductivity field \eqn{K} \code{[L/T]}, length-\code{n} vector.
#' @param iobs    scalar integer, the flat element index of the observation point.
#' @param lrw     real work array length for \code{steady.3D} (default 20000000).
#' @return Numeric vector of adjoint variable \eqn{\lambda} (length \code{n}).
#' @keywords internal
adjoint3D <- function(grid, KK, iobs, lrw = 20000000) {
  require('rootSolve')

  n   <- grid$n
  nx  <- grid$nx
  ny  <- grid$ny
  nz  <- grid$nz
  dx  <- grid$dx
  dy  <- grid$dy
  dz  <- grid$dz

  Kp <- if (length(KK) == 1) rep(KK, n) else KK

  # Convert flat index iobs to (Nxp, Nyp, Nzp)
  Nxp <- ((iobs - 1) %% nx) + 1
  Nyp <- (((iobs - 1) %/% nx) %% ny) + 1
  Nzp <- ((iobs - 1) %/% (nx * ny)) + 1

  # Adjoint source: unit extraction at the observation cell
  Qp_adj <- (-dx * dy * dz)

  para <- list(dx = dx, dy = dy, dz = dz,
               nx = nx, ny = ny, nz = nz,
               Kp = Kp, Qp = Qp_adj,
               Nxp = Nxp, Nyp = Nyp, Nzp = Nzp)

  y_init <- rep(0, n)
  s_adj <- steady.3D(y = y_init, parms = para,
                     func = diffusion3D_GW, dimens = c(nx, ny, nz), lrw = lrw)$y

  return(s_adj)
}


#' Compute the Jacobian (sensitivity) matrix via the adjoint-state method (3D)
#'
#' 3D counterpart of \code{\link{jacobian2D}}.  For 3D steady-state hydraulic
#' tomography, computes the Jacobian matrix
#' \eqn{\mathbf{J} = \partial \mathbf{h} / \partial (\ln \mathbf{K})}
#' (dimension \code{nobs × nelem}) using the adjoint-state method.
#'
#' For each pumping test \eqn{j} and each observation \eqn{i}:
#' \enumerate{
#'   \item Solve the forward problem to obtain head field \eqn{h}
#'   \item Solve the adjoint problem to obtain \eqn{\lambda_i}
#'   \item Compute sensitivity:
#'         \eqn{J_{ik} = -K_k \, (\nabla h|_k \cdot \nabla \lambda_i|_k) \, \Delta x \, \Delta y \, \Delta z}
#' }
#'
#' @param grid    grid list from \code{GenGrid3D()}.
#' @param KK      hydraulic conductivity field \eqn{K} \code{[L/T]}, length-\code{n} vector.
#' @param qHT     list of pumping test data frames (columns \code{Qp, x, y, z}).
#' @param oHT     list of observation data frames (columns \code{data, x, y, z}).
#' @param lrw     real work array length for \code{steady.3D} (default 20000000).
#' @return A matrix of dimension \code{nobs × nelem}.
#' @export
jacobian3D <- function(grid,
                       KK,
                       qHT = list(data.frame(Qp = 10, x = 7.5, y = 7.5, z = 2.5)),
                       oHT = list(data.frame(data = -1, x = 5, y = 5, z = 2.5)),
                       lrw = 20000000) {

  require('rootSolve')

  nHT  <- length(qHT)
  n    <- grid$n
  nx   <- grid$nx
  ny   <- grid$ny
  nz   <- grid$nz
  dx   <- grid$dx
  dy   <- grid$dy
  dz   <- grid$dz
  dV   <- dx * dy * dz

  # ---- Collect observation element indices and counts ----
  loc_obsHT <- list()
  nobs_per_test <- integer(nHT)
  for (j in seq_len(nHT)) {
    oinf <- oHT[[j]]
    oelemdf <- getOelem3D(grid = grid, Oinf = oinf)
    loc_obsHT[[j]] <- oelemdf$nelem
    nobs_per_test[j] <- nrow(oinf)
  }
  nobs <- sum(nobs_per_test)

  Kp <- if (length(KK) == 1) rep(KK, n) else KK

  # ---- Pre-allocate Jacobian ----
  J <- matrix(0, nrow = nobs, ncol = n)

  # ---- Row offset tracker ----
  row_offset <- 0

  for (j in seq_len(nHT)) {
    qinf <- qHT[[j]]
    obs_elem <- loc_obsHT[[j]]
    nlocal <- length(obs_elem)

    # --- Forward solve for pumping test j ---
    h_forward <- Fsteady3dsim(grid = grid, KK = Kp, Qinf = qinf, lrw = lrw)$solution
    h_arr <- array(h_forward, dim = c(nx, ny, nz))

    # --- Compute forward head gradients (cell-centred) ---
    dhdx <- array(0, dim = c(nx, ny, nz))
    dhdy <- array(0, dim = c(nx, ny, nz))
    dhdz <- array(0, dim = c(nx, ny, nz))

    if (nx > 1) {
      dhdx[2:(nx - 1), , ] <- (h_arr[3:nx, , ] - h_arr[1:(nx - 2), , ]) / (2 * dx)
      dhdx[1, , ]  <- (h_arr[2, , ] - h_arr[1, , ]) / dx
      dhdx[nx, , ] <- (h_arr[nx, , ] - h_arr[nx - 1, , ]) / dx
    }
    if (ny > 1) {
      dhdy[, 2:(ny - 1), ] <- (h_arr[, 3:ny, ] - h_arr[, 1:(ny - 2), ]) / (2 * dy)
      dhdy[, 1, ]  <- (h_arr[, 2, ] - h_arr[, 1, ]) / dy
      dhdy[, ny, ] <- (h_arr[, ny, ] - h_arr[, ny - 1, ]) / dy
    }
    if (nz > 1) {
      dhdz[, , 2:(nz - 1)] <- (h_arr[, , 3:nz] - h_arr[, , 1:(nz - 2)]) / (2 * dz)
      dhdz[, , 1]  <- (h_arr[, , 2] - h_arr[, , 1]) / dz
      dhdz[, , nz] <- (h_arr[, , nz] - h_arr[, , nz - 1]) / dz
    }

    # --- For each observation in this test, solve adjoint and compute sensitivity ---
    for (i_local in seq_len(nlocal)) {
      iobs <- obs_elem[i_local]
      global_row <- row_offset + i_local

      # Solve adjoint equation
      lambda <- adjoint3D(grid = grid, KK = Kp, iobs = iobs, lrw = lrw)
      lambda_arr <- array(lambda, dim = c(nx, ny, nz))

      # Compute adjoint gradients
      dldx <- array(0, dim = c(nx, ny, nz))
      dldy <- array(0, dim = c(nx, ny, nz))
      dldz <- array(0, dim = c(nx, ny, nz))

      if (nx > 1) {
        dldx[2:(nx - 1), , ] <- (lambda_arr[3:nx, , ] - lambda_arr[1:(nx - 2), , ]) / (2 * dx)
        dldx[1, , ]  <- (lambda_arr[2, , ] - lambda_arr[1, , ]) / dx
        dldx[nx, , ] <- (lambda_arr[nx, , ] - lambda_arr[nx - 1, , ]) / dx
      }
      if (ny > 1) {
        dldy[, 2:(ny - 1), ] <- (lambda_arr[, 3:ny, ] - lambda_arr[, 1:(ny - 2), ]) / (2 * dy)
        dldy[, 1, ]  <- (lambda_arr[, 2, ] - lambda_arr[, 1, ]) / dy
        dldy[, ny, ] <- (lambda_arr[, ny, ] - lambda_arr[, ny - 1, ]) / dy
      }
      if (nz > 1) {
        dldz[, , 2:(nz - 1)] <- (lambda_arr[, , 3:nz] - lambda_arr[, , 1:(nz - 2)]) / (2 * dz)
        dldz[, , 1]  <- (lambda_arr[, , 2] - lambda_arr[, , 1]) / dz
        dldz[, , nz] <- (lambda_arr[, , nz] - lambda_arr[, , nz - 1]) / dz
      }

      # Sensitivity: J_{ik} = -K_k * (grad_h_k · grad_lambda_i_k) * dV
      grad_dot <- as.vector(dhdx * dldx + dhdy * dldy + dhdz * dldz)
      J[global_row, ] <- -Kp * grad_dot * dV
    }

    row_offset <- row_offset + nlocal
  }

  return(J)
}


#' Compute the Jacobian for 3D well-screen observations
#'
#' Like \code{\link{jacobian3D}} but for well-screen (vertical interval)
#' observations.  Each observation is the arithmetic average of heads over
#' a set of grid cells in the screen interval.  The sensitivity row for such
#' an observation is the weighted average of the point-sensitivities over
#' those cells.
#'
#' @param grid    grid list from \code{GenGrid3D()}.
#' @param KK      hydraulic conductivity field \eqn{K} \code{[L/T]}, length-\code{n} vector.
#' @param qHT     list of pumping test data frames with columns
#'   \code{Qp, x, y, z_top, z_bottom} (well-screen pumping) or
#'   \code{Qp, x, y, z} (point pumping).
#' @param oHT     list of observation data frames with columns
#'   \code{data, x, y, z_top, z_bottom}.
#' @param lrw     real work array length for \code{steady.3D} (default 20000000).
#' @return A matrix of dimension \code{nobs × nelem}.
#' @export
jacobian3DScreen <- function(grid,
                             KK,
                             qHT = list(data.frame(Qp = 10, x = 7.5, y = 7.5,
                                                   z_top = 2, z_bottom = 3)),
                             oHT = list(data.frame(data = -1, x = 5, y = 5,
                                                   z_top = 1, z_bottom = 4)),
                             lrw = 20000000) {

  require('rootSolve')

  nHT  <- length(qHT)
  n    <- grid$n
  nx   <- grid$nx
  ny   <- grid$ny
  nz   <- grid$nz
  dx   <- grid$dx
  dy   <- grid$dy
  dz   <- grid$dz
  dV   <- dx * dy * dz

  # ---- For each observation, determine which grid cells belong to the screen ----
  obs_elem_list <- list()   # list of cell-index vectors, one per observation
  nobs_per_test <- integer(nHT)
  for (j in seq_len(nHT)) {
    oinf <- oHT[[j]]
    nobs_per_test[j] <- nrow(oinf)
    for (r in seq_len(nrow(oinf))) {
      x_obs <- oinf$x[r]
      y_obs <- oinf$y[r]
      z_lo  <- min(oinf$z_top[r], oinf$z_bottom[r])
      z_hi  <- max(oinf$z_top[r], oinf$z_bottom[r])

      ix <- which.min(abs(grid$xmid - x_obs))
      iy <- which.min(abs(grid$ymid - y_obs))
      iz_vec <- which(grid$zmid >= z_lo & grid$zmid <= z_hi)
      if (length(iz_vec) == 0) {
        iz_vec <- which.min(abs(grid$zmid - mean(c(z_lo, z_hi))))
      }

      # flat element indices for cells in this screen
      elems <- integer(0)
      for (iz in iz_vec) {
        elems <- c(elems, (iz - 1L) * nx * ny + (iy - 1L) * nx + ix)
      }
      obs_elem_list[[length(obs_elem_list) + 1]] <- elems
    }
  }
  nobs <- sum(nobs_per_test)

  Kp <- if (length(KK) == 1) rep(KK, n) else KK

  # ---- Pre-allocate Jacobian ----
  J <- matrix(0, nrow = nobs, ncol = n)

  # ---- Global observation index tracker ----
  obs_global <- 0

  for (j in seq_len(nHT)) {
    qinf <- qHT[[j]]
    nlocal <- nobs_per_test[j]

    # --- Forward solve for pumping test j ---
    h_forward <- Fsteady3dsim(grid = grid, KK = Kp, Qinf = qinf, lrw = lrw)$solution
    h_arr <- array(h_forward, dim = c(nx, ny, nz))

    # --- Compute forward head gradients (cell-centred) ---
    dhdx <- array(0, dim = c(nx, ny, nz))
    dhdy <- array(0, dim = c(nx, ny, nz))
    dhdz <- array(0, dim = c(nx, ny, nz))

    if (nx > 1) {
      dhdx[2:(nx - 1), , ] <- (h_arr[3:nx, , ] - h_arr[1:(nx - 2), , ]) / (2 * dx)
      dhdx[1, , ]  <- (h_arr[2, , ] - h_arr[1, , ]) / dx
      dhdx[nx, , ] <- (h_arr[nx, , ] - h_arr[nx - 1, , ]) / dx
    }
    if (ny > 1) {
      dhdy[, 2:(ny - 1), ] <- (h_arr[, 3:ny, ] - h_arr[, 1:(ny - 2), ]) / (2 * dy)
      dhdy[, 1, ]  <- (h_arr[, 2, ] - h_arr[, 1, ]) / dy
      dhdy[, ny, ] <- (h_arr[, ny, ] - h_arr[, ny - 1, ]) / dy
    }
    if (nz > 1) {
      dhdz[, , 2:(nz - 1)] <- (h_arr[, , 3:nz] - h_arr[, , 1:(nz - 2)]) / (2 * dz)
      dhdz[, , 1]  <- (h_arr[, , 2] - h_arr[, , 1]) / dz
      dhdz[, , nz] <- (h_arr[, , nz] - h_arr[, , nz - 1]) / dz
    }

    # --- For each screen observation, average point-sensitivities ---
    for (i_local in seq_len(nlocal)) {
      obs_global <- obs_global + 1
      screen_elems <- obs_elem_list[[obs_global]]
      ns <- length(screen_elems)

      J_row <- rep(0, n)

      for (k in seq_len(ns)) {
        iobs <- screen_elems[k]

        # Solve adjoint for this cell
        lambda <- adjoint3D(grid = grid, KK = Kp, iobs = iobs, lrw = lrw)
        lambda_arr <- array(lambda, dim = c(nx, ny, nz))

        dldx <- array(0, dim = c(nx, ny, nz))
        dldy <- array(0, dim = c(nx, ny, nz))
        dldz <- array(0, dim = c(nx, ny, nz))

        if (nx > 1) {
          dldx[2:(nx - 1), , ] <- (lambda_arr[3:nx, , ] - lambda_arr[1:(nx - 2), , ]) / (2 * dx)
          dldx[1, , ]  <- (lambda_arr[2, , ] - lambda_arr[1, , ]) / dx
          dldx[nx, , ] <- (lambda_arr[nx, , ] - lambda_arr[nx - 1, , ]) / dx
        }
        if (ny > 1) {
          dldy[, 2:(ny - 1), ] <- (lambda_arr[, 3:ny, ] - lambda_arr[, 1:(ny - 2), ]) / (2 * dy)
          dldy[, 1, ]  <- (lambda_arr[, 2, ] - lambda_arr[, 1, ]) / dy
          dldy[, ny, ] <- (lambda_arr[, ny, ] - lambda_arr[, ny - 1, ]) / dy
        }
        if (nz > 1) {
          dldz[, , 2:(nz - 1)] <- (lambda_arr[, , 3:nz] - lambda_arr[, , 1:(nz - 2)]) / (2 * dz)
          dldz[, , 1]  <- (lambda_arr[, , 2] - lambda_arr[, , 1]) / dz
          dldz[, , nz] <- (lambda_arr[, , nz] - lambda_arr[, , nz - 1]) / dz
        }

        grad_dot <- as.vector(dhdx * dldx + dhdy * dldy + dhdz * dldz)
        J_row <- J_row + (-Kp * grad_dot * dV)
      }

      J[obs_global, ] <- J_row / ns
    }
  }

  return(J)
}


# ==============================================================================
# Adjoint-state sensitivity computation for 2D steady-state HT
# ==============================================================================

#' Solve the adjoint equation for a single observation point (2D steady-state)
#'
#' For the 2D steady-state groundwater flow equation
#' \eqn{\nabla \cdot (T \nabla h) = Q}, the adjoint state \eqn{\lambda} for an
#' observation at cell \code{iobs} satisfies:
#' \deqn{\nabla \cdot (T \nabla \lambda) = \delta(\mathbf{x} - \mathbf{x}_{iobs})}
#' with homogeneous (zero) boundary conditions.  This is the same PDE as the
#' forward problem but with a unit source at the observation location instead
#' of a pumping well.  The adjoint variable is used to compute the Jacobian
#' (sensitivity) matrix via:
#' \deqn{\frac{\partial h_{iobs}}{\partial (\ln T)_k} =
#'       -T_k \, \nabla h|_k \cdot \nabla \lambda|_k \, \Delta x \, \Delta y}
#'
#' @param grid    grid list from \code{GenGrid()}.
#' @param TT      transmissivity field \eqn{T} \code{[L^2/T]}, length-\code{n} vector.
#' @param iobs    scalar integer, the flat element index of the observation point.
#' @param lrw     real work array length for \code{steady.2D} (default 160000).
#' @return Numeric vector of adjoint variable \eqn{\lambda} (length \code{n}).
#' @keywords internal
adjoint2D <- function(grid, TT, iobs, lrw = 160000) {
  require('rootSolve')

  n   <- grid$n
  nx  <- grid$nx
  ny  <- grid$ny
  dx  <- grid$dx
  dy  <- grid$dy

  Tp <- if (length(TT) == 1) rep(TT, n) else TT
  Sp <- rep(1, n)
  Tpm <- matrix(Tp, nx, ny)
  Spm <- matrix(Sp, nx, ny)

  # Build the adjoint source: unit extraction at the observation cell
  # The adjoint PDE has source +delta(x - x_obs) on the RHS.
  # In the diffusion2D residual, Qp/(dx*dy*S) is subtracted from dY at the pump cell.
  # For adjoint, we want a unit source, so we set Qp = dx*dy (S=1 for steady-state)
  # and place it at the observation cell.

  # Convert flat index iobs to (Nxp, Nyp)
  Nxp <- ((iobs - 1) %% nx) + 1
  Nyp <- ((iobs - 1) %/% nx) + 1

  # Adjoint source: Qp = dx*dy so that Qp/(dx*dy*Spm) = 1 at the obs cell
  ### note that Qp_adj is treated as pumping out, so we need to negative it first.

  Qp_adj <- (-dx * dy) ### late we treat this as pumping...

  para <- list(dx = dx, dy = dy, nx = nx, ny = ny,
               Tpm = Tpm, Spm = Spm,
               Qp = Qp_adj, Nxp = Nxp, Nyp = Nyp)

  y_init <- rep(0, n)
  s_adj <- steady.2D(y = y_init, parms = para,
                     func = diffusion2D, dimens = c(nx, ny), lrw = lrw)$y

  return(s_adj)
}


#' Compute the Jacobian (sensitivity) matrix via the adjoint-state method
#'
#' For 2D steady-state hydraulic tomography, computes the Jacobian matrix
#' \eqn{\mathbf{J} = \partial \mathbf{h} / \partial (\ln \mathbf{T})}
#' (dimension \code{nobs × nelem}) using the adjoint-state method.
#'
#' For each pumping test \eqn{j} and each observation \eqn{i}:
#' \enumerate{
#'   \item Solve the forward problem to obtain head field \eqn{h}
#'   \item Solve the adjoint problem to obtain \eqn{\lambda_i}
#'   \item Compute sensitivity:
#'         \eqn{J_{ik} = -T_k \, (\nabla h|_k \cdot \nabla \lambda_i|_k) \, \Delta x \, \Delta y}
#' }
#'
#' @param grid    grid list from \code{GenGrid()}.
#' @param TT      transmissivity field \eqn{T} \code{[L^2/T]}, length-\code{n} vector.
#' @param qHT     list of pumping test data frames (columns \code{Qp, x, y}).
#' @param oHT     list of observation data frames (columns \code{data, x, y}).
#' @param lrw     real work array length for \code{steady.2D} (default 160000).
#' @return A matrix of dimension \code{nobs × nelem}.
#' @examples
#' # Setup grid and transmissivity field
#' domain <- c(20, 20, 0, 20, 0, 20)
#' grid   <- GenGrid(domain)
#' set.seed(123)
#' TT <- random2d(grid = grid, nsim = 1)$Tp
#'
#' # Define pumping tests and observations
#' qHT <- list(
#'   data.frame(Qp = 10, x = 10.5, y = 10.5),
#'   data.frame(Qp = 10, x = 15.5, y = 15.5)
#' )
#' oHT <- list(
#'   data.frame(data = -1, x = 5.5, y = 5.5),
#'   data.frame(data = -1, x = 12.5, y = 12.5)
#' )
#'
#' # Compute Jacobian
#' J <- jacobian2D(grid = grid, TT = TT, qHT = qHT, oHT = oHT)
#'
#' # Inspect results
#' # Jacobian dimensions: nobs × nelem
#' dim(J)
#' # Sensitivity map for the first observation (row 1)
#' image(matrix(J[1, ], grid$nx, grid$ny),
#'       main = "Sensitivity: obs 1")
#' @export
#' 
jacobian2D <- function(grid,
                       TT,
                       qHT = list(data.frame(Qp = 10, x = 20.5, y = 20.5)),
                       oHT = list(data.frame(data = -1, x = 11, y = 11)),
                       lrw = 160000) 
{

  require('rootSolve')

  nHT  <- length(qHT)
  n    <- grid$n
  nx   <- grid$nx
  ny   <- grid$ny
  dx   <- grid$dx
  dy   <- grid$dy
  dV   <- dx * dy   # cell volume (area in 2D)

  # ---- Collect observation element indices and counts ----
  loc_obsHT <- list()
  nobs_per_test <- integer(nHT)
  for (j in seq_len(nHT)) {
    oinf <- oHT[[j]]
    oelemdf <- getOelem(grid = grid, Oinf = oinf)
    loc_obsHT[[j]] <- oelemdf$nelem
    nobs_per_test[j] <- nrow(oinf)
  }
  nobs <- sum(nobs_per_test)

  Tp <- if (length(TT) == 1) rep(TT, n) else TT

  # ---- Pre-allocate Jacobian ----
  J <- matrix(0, nrow = nobs, ncol = n)

  # ---- Row offset tracker ----
  row_offset <- 0

  for (j in seq_len(nHT)) {
    qinf <- qHT[[j]]
    obs_elem <- loc_obsHT[[j]]
    nlocal <- length(obs_elem)

    # --- Forward solve for pumping test j ---
    h_forward <- Fsteady2dsim(grid = grid, TT = Tp, Qinf = qinf, lrw = lrw)$solution
    h_mat <- matrix(h_forward, nx, ny)

    # --- Compute forward head gradients (cell-centred) ---
    dhdx <- matrix(0, nx, ny)
    dhdy <- matrix(0, nx, ny)

    # Interior x-derivative (central difference)
    if (nx > 1) {
      dhdx[2:(nx - 1), ] <- (h_mat[3:nx, ] - h_mat[1:(nx - 2), ]) / (2 * dx)
      dhdx[1, ]  <- (h_mat[2, ] - h_mat[1, ]) / dx
      dhdx[nx, ] <- (h_mat[nx, ] - h_mat[nx - 1, ]) / dx
    }

    # Interior y-derivative (central difference)
    if (ny > 1) {
      dhdy[, 2:(ny - 1)] <- (h_mat[, 3:ny] - h_mat[, 1:(ny - 2)]) / (2 * dy)
      dhdy[, 1]  <- (h_mat[, 2] - h_mat[, 1]) / dy
      dhdy[, ny] <- (h_mat[, ny] - h_mat[, ny - 1]) / dy
    }

    # --- For each observation in this test, solve adjoint and compute sensitivity ---
    for (i_local in seq_len(nlocal)) {
      iobs <- obs_elem[i_local]
      global_row <- row_offset + i_local

      # Solve adjoint equation
      lambda <- adjoint2D(grid = grid, TT = Tp, iobs = iobs, lrw = lrw)
      lambda_mat <- matrix(lambda, nx, ny)

      # Compute adjoint gradients
      dldx <- matrix(0, nx, ny)
      dldy <- matrix(0, nx, ny)

      if (nx > 1) {
        dldx[2:(nx - 1), ] <- (lambda_mat[3:nx, ] - lambda_mat[1:(nx - 2), ]) / (2 * dx)
        dldx[1, ]  <- (lambda_mat[2, ] - lambda_mat[1, ]) / dx
        dldx[nx, ] <- (lambda_mat[nx, ] - lambda_mat[nx - 1, ]) / dx
      }
      if (ny > 1) {
        dldy[, 2:(ny - 1)] <- (lambda_mat[, 3:ny] - lambda_mat[, 1:(ny - 2)]) / (2 * dy)
        dldy[, 1]  <- (lambda_mat[, 2] - lambda_mat[, 1]) / dy
        dldy[, ny] <- (lambda_mat[, ny] - lambda_mat[, ny - 1]) / dy
      }

      # Sensitivity: J_{ik} = -T_k * (grad_h_k · grad_lambda_i_k) * dV
      grad_dot <- as.vector(dhdx * dldx + dhdy * dldy)
      J[global_row, ] <- -Tp * grad_dot * dV
    }

    row_offset <- row_offset + nlocal
  }

  return(J)
}


#' 2D Adjoint-based deterministic inversion for hydraulic tomography
#'
#' Performs deterministic Bayesian inversion to estimate the 2D
#' log-transmissivity field \eqn{\ln T(\mathbf{x})} from hydraulic tomography
#' (HT) data using the **adjoint-state method** to compute the Jacobian
#' (sensitivity) matrix analytically, instead of using ensemble statistics
#' as in the ESMDA-based \code{\link{Finverse3}}.
#'
#' The update follows the quasi-linear Bayesian formulation:
#' \deqn{\Delta \mathbf{m} = (\mathbf{J}^\top \mathbf{R}^{-1} \mathbf{J} +
#'       \mathbf{C}_{kk}^{-1})^{-1} \,
#'       \mathbf{J}^\top \mathbf{R}^{-1} \,
#'       (\mathbf{d}_{\text{obs}} - \mathbf{g}(\mathbf{m}))}
#' where \eqn{\mathbf{J}} is the Jacobian from \code{\link{jacobian2D}},
#' \eqn{\mathbf{R}} is the observation error covariance (assumed diagonal with
#' variance \code{sigma2obs}), and \eqn{\mathbf{C}_{kk}} is the prior parameter
#' covariance derived from the geostatistical model.
#'
#' A Levenberg-Marquardt damping parameter \code{lambda_lm} is added to the
#' diagonal of the Hessian approximation for stability.
#'
#' @param domain   a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)} \code{[L]}.
#' @param grid     grid list from \code{GenGrid()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT      list of pumping test data frames (columns \code{Qp, x, y}).
#' @param itermax  maximum number of iterations (default 5).
#' @param rmsemin  minimum RMSE stopping criterion (default 0).
#' @param oHT      list of observation data frames (columns \code{data, x, y}).
#' @param lrw      real work array length for \code{steady.2D} (default 160000).
#' @param sigma2obs observation error variance (default 1e-4).
#' @param lambda_lm Levenberg-Marquardt damping parameter (default 1.0); decays
#'   by \code{decay_lm} each iteration.
#' @param decay_lm  decay factor for \code{lambda_lm} (default 2.0).
#' @param geo      prior geostatistical parameters:
#'   \code{list(me, var, geomod, anis, range, nugget)}.
#' @param m0       initial log-T field, length \code{n} vector.  If \code{NULL}
#'   (default), uses the prior mean \code{geo$me} everywhere.
#' @param ifcor    logical; if \code{TRUE}, returns the Jacobian matrix and
#'   forward response after the first iteration without performing the update.
#' @param ifvarTh  logical; if \code{TRUE}, computes and stores the predicted
#'   head variance \code{varobsh} and the conditional (posterior) variance of
#'   log-T \code{varT} at each iteration.  \code{varobsh} is
#'   \code{diag(J \%*\% C_kk \%*\% t(J))}; \code{varT} follows the SLE
#'   conditional-covariance update (Yeh & Liu, 2000) and is analogous to a
#'   kriging variance.  Default \code{FALSE} because the covariance matrix
#'   operation can be expensive for large grids.
#' @return A list of per-iteration results, each a list with
#'   \code{meanT} (log-T estimate, same as \code{m}), \code{varT}
#'   (conditional variance of log-T, only when \code{ifvarTh = TRUE}),
#'   \code{meanobsh} (simulated heads), \code{varobsh} (predicted head variance,
#'   only when \code{ifvarTh = TRUE}),
#'   \code{m}, \code{h_sim}, \code{rmse}, \code{l2}, \code{l1}, and
#'   \code{J} (Jacobian, first iteration only).  The variable names
#'   \code{meanT}, \code{varT}, \code{meanobsh}, and \code{varobsh} are kept
#'   compatible with \code{\link{Finverse3}} so that \code{\link{inversePlot}}
#'   can be used directly.  When \code{ifcor = TRUE}, returns
#'   \code{list(J = J, h_sim = h_sim)}.
#' @export
#' @examples
#' \dontrun{
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT <- trueK[,-c(1,2)]
#' domain <- c(40,40,0,40,0,40)
#' grid <- GenGrid(domain)
#' Qinf1 <- data.frame(Qp=10, x=20.5, y=20.5)
#' qHT <- list(test1 = Qinf1)
#' trueh <- Fsteady2dsim(TT=TT, Qinf=Qinf1, domain=domain)
#' locx <- c(15,18,22,25,30)
#' locy <- c(15,18,22,25,30)
#' loc <- expand.grid(x=locx, y=locy)
#' Oinf1 <- data.frame(data=NA, x=loc$x, y=loc$y)
#' Oinf1 <- samData(grid=grid, Oinf=Oinf1, h=trueh$solution)
#' oHT <- list(test1 = Oinf1)
#' result <- FinverseAdj(grid=grid, qHT=qHT, oHT=oHT, itermax=3)
#' }
FinverseAdj <- function(
    domain      = c(40, 40, 0, 40, 0, 40),
    grid        = NULL,
    qHT         = list(data.frame(Qp = 10, x = 20.5, y = 20.5)),
    itermax     = 5,
    rmsemin     = 0,
    oHT         = list(data.frame(data = -1, x = 11, y = 11)),
    lrw         = 160000,
    sigma2obs   = 1e-4,
    lambda_lm   = 1.0,
    decay_lm    = 2.0,
    geo         = list(me = 0, var = 1, geomod = "Exp",
                       anis = c(90, 1), range = 30, nugget = 0),
    m0          = NULL,
    ifcor       = FALSE,
    ifvarTh     = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid(domain)

  n <- grid$n
  nHT <- length(qHT)

  # ---- Extract observation data ----
  trueobshHT <- list()
  loc_obsHT  <- list()
  for (i in seq_len(nHT)) {
    oinf <- oHT[[i]]
    oelemdf <- getOelem(grid = grid, Oinf = oinf)
    loc_obsHT[[i]] <- oelemdf$nelem
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs <- length(trueobsh)

  # ---- Build prior covariance C_kk from geostatistical model ----
  require('gstat')
  xy <- grid$grid

  dist_mat <- as.matrix(dist(xy))
  if (geo$geomod == "Exp") {
    C_kk <- geo$var * exp(-dist_mat / geo$range)
    diag(C_kk) <- geo$var + geo$nugget
  } else if (geo$geomod == "Gau") {
    C_kk <- geo$var * exp(-(dist_mat / geo$range)^2)
    diag(C_kk) <- geo$var + geo$nugget
  } else if (geo$geomod == "Sph") {
    hr <- dist_mat / geo$range
    C_kk <- geo$var * (1 - 1.5 * hr + 0.5 * hr^3)
    C_kk[hr > 1] <- 0
    diag(C_kk) <- geo$var + geo$nugget
  } else {
    C_kk <- geo$var * exp(-dist_mat / geo$range)
    diag(C_kk) <- geo$var + geo$nugget
  }

  # ---- Initial parameter estimate ----
  if (is.null(m0)) {
    m <- rep(geo$me, n)  # log-T
  } else {
    m <- m0
  }
  T_current <- exp(m)

  # ---- Iteration loop ----
  niter <- 1
  rmse <- 1e10
  iterdf <- list()

  while (niter <= itermax && rmse > rmsemin) {

    # --- Forward simulation at current T ---
    h_sim_vec <- samHT(grid = grid, TT = T_current,
                       qHT = qHT, oHT = oHT, lrw = lrw)

    # --- Misfit ---
    residual <- trueobsh - h_sim_vec
    rmse <- sqrt(mean(residual^2))
    l2 <- rmse
    l1 <- mean(abs(residual))

    msg <- paste('niter =', niter,
                 'rmse =', round(rmse, 6),
                 'l2 =', round(l2, 6),
                 'l1 =', round(l1, 6),
                 'lambda_lm =', round(lambda_lm, 4))
    print(msg)

    # --- Compute Jacobian ---
    J <- jacobian2D(grid = grid, TT = T_current,
                    qHT = qHT, oHT = oHT, lrw = lrw)

    if (ifcor) {
      covhk <- J %*% C_kk
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }

    # --- Compute update (dual / observation-space formulation) ---
    Jt <- t(J)                                     # nelem x nobs
    JCJt <- J %*% C_kk %*% Jt                      # nobs x nobs
    S <- JCJt + (sigma2obs + lambda_lm) * diag(nobs)

    beta <- solve(S, residual)                     # nobs
    dm <- as.vector(C_kk %*% Jt %*% beta)          # nelem

    # --- Estimate variances ---
    meanobsh <- h_sim_vec

    if (ifvarTh) {
      # Conditional (posterior) parameter covariance: C_kk - C_kk J^T S^{-1} J C_kk
      # (similar to kriging conditional variance; Yeh & Liu, 2000, eq. 4)
      KS <- C_kk %*% Jt %*% solve(S)               # nelem x nobs (Kalman-like gain)
      C_post <- C_kk - KS %*% J %*% C_kk           # nelem x nelem
      varT <- as.vector(diag(C_post))
      # Predicted head variance = diag(J C_kk J^T)
      JCJt_post <- J %*% C_post %*% Jt
      varobsh <- as.vector(diag(JCJt_post))
    } else {
      varT <- NULL
      varobsh <- NULL
    }

    # --- Update parameters ---
    m <- m + dm
    T_current <- exp(m)

    # --- Store iteration results ---
    # Use the same variable names as Finverse3 for compatibility with inversePlot
    iterdf[[niter]] <- list(
      meanT    = as.vector(m),
      varT     = varT,
      meanobsh = meanobsh,
      varobsh  = varobsh,
      m        = m,
      h_sim    = h_sim_vec,
      rmse     = rmse,
      l2       = l2,
      l1       = l1
    )

    if (niter == 1) {
      iterdf[[niter]]$J <- J
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    lambda_lm <- lambda_lm / decay_lm
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}


# ==============================================================================
# 2D Transient Adjoint-based Inversion
# ==============================================================================

# ---- Helper: build a "good" time grid (fine at start, coarser later) ----------
# For groundwater pumping problems, head changes rapidly at early times and
# slowly at late times.  A geometric progression in dt captures this efficiently.
#
# @param t_max   end time of the simulation (positive).
# @param npoints desired number of output time points (default 101).
# @param p       ratio of successive dt's (>1 ⇒ dt grows; <1 ⇒ dt shrinks).
#                Default 1.05 means each step is 5% longer than the previous.
# @return A strictly increasing numeric vector starting at 0 and ending
#         approximately at t_max.
.makeTimeGrid <- function(t_max, npoints = 101, p = 1.05) {
  if (t_max <= 0) t_max <- 1
  # Solve: dt1 * (p^0 + p^1 + ... + p^(npoints-2)) = t_max
  if (abs(p - 1) < 1e-8) {
    times <- seq(0, t_max, length.out = npoints)
  } else {
    # geometric sum: dt1 * (p^(npoints-1) - 1) / (p - 1) = t_max
    dt1 <- t_max * (p - 1) / (p^(npoints - 1) - 1)
    dts  <- dt1 * p^(seq(0, npoints - 2))
    times <- cumsum(c(0, dts))
    # ensure exact end
    times <- times / max(times) * t_max
  }
  return(times)
}


#' 2D Transient adjoint (Green's function) solver
#'
#' Computes the adjoint variable \eqn{\varphi(\mathbf{x}, t)} for a 2D
#' transient groundwater flow problem by solving the **forward** Green's
#' function equation with a unit injection at the observation cell
#' \eqn{\mathbf{x}_o}:
#' \deqn{S \frac{\partial \varphi}{\partial t} =
#'       \nabla \cdot (T \nabla \varphi) + \delta(\mathbf{x}-\mathbf{x}_o)\delta(t)}
#' subject to zero initial and boundary conditions.
#'
#' In practice, the impulse source is replaced by a continuous unit injection
#' (Q = -1 in the code convention), giving the step response \eqn{\Phi}.
#' The impulse response (adjoint) is then obtained by finite-differencing
#' \eqn{\Phi} in time:
#' \deqn{\varphi(\mathbf{x}, t) \approx \frac{\partial \Phi}{\partial t}}
#'
#' The time grid for the adjoint is generated independently of observation
#' times — it only needs to span from 0 to the end of pumping (\code{t_max}).
#' The convolution in \code{jacobian2DTr} handles time alignment via
#' interpolation.
#'
#' For a fixed observation location, the same \eqn{\varphi} time series can be
#' reused for every observation time at that location (Zha et al., 2020).
#'
#' @param grid  grid list from \code{GenGrid()}.
#' @param TT    transmissivity field \eqn{T} \code{[L^2/T]}, length-\code{n} vector.
#' @param SS    storage coefficient \eqn{S} [-], length-\code{n} vector or scalar.
#' @param iobs  flat index of the observation cell.
#' @param times numeric vector of output times for the ODE solver.
#'   If \code{NULL} (default), a geometric time grid is generated
#'   using \code{t_max} and \code{npoints}.
#' @param t_max end time of the adjoint simulation.  Only used when
#'   \code{times = NULL}.  Default: 100.
#' @param npoints number of time points for auto-generated grid (default 101).
#' @param p_grid growth factor for auto-generated grid (default 1.05).
#' @param lrw   real work array length for \code{ode.2D} (default 1600000).
#' @return A list with \code{phi} (matrix \code{nt \times n} of adjoint values) and
#'   \code{times} (the output time vector).
#' @keywords internal
adjoint2DTr <- function(grid, TT, SS, iobs, times = NULL,
                         t_max = 100, npoints = 101, p_grid = 1.05,
                         lrw = 1600000) {

  require('deSolve')

  # ---- Build time grid if not supplied ----------------------------------------
  if (is.null(times)) {
    times <- .makeTimeGrid(t_max, npoints, p_grid)
  } else {
    # Merge with a geometric grid to ensure early-time resolution for the
    # impulse response (adjoint), which varies rapidly near t = 0.
    auto_times <- .makeTimeGrid(max(times), npoints, p_grid)
    times <- sort(unique(c(times, auto_times)))
    if (length(times) < 2) {
      times <- seq(0, max(times, 1), length.out = 3)
    }
  }

  # Unit injection at the observation cell.
  # Note that Qp should be dx * dy, zyy20260718
dx <- grid$dx
dy <- grid$dy
qinf <- data.frame(Qp = -dx * dy, x = grid$grid$x[iobs], y = grid$grid$y[iobs]) 
#  qinf <- data.frame(Qp = -1, x = grid$grid$x[iobs], y = grid$grid$y[iobs])

  res <- Ftransient2dsim(grid = grid, TT = TT, SS = SS,
                          Qinf = qinf, times = times, lrw = lrw)

  out <- res$out            # nt x (n+1)
  sim_times <- out[, 1]
  Phi <- out[, -1, drop = FALSE]  # step response, nt x n
  nt <- nrow(out)
  n <- ncol(Phi)

  # Impulse response = time derivative of step response
  phi <- matrix(0, nt, n)

  if (nt == 1) {
    phi[1, ] <- 0
  } else {
    # Forward difference for all intervals
    phi[-nt, ] <- (Phi[2:nt, ] - Phi[1:(nt - 1), ]) / diff(sim_times)
    # Extend the last value
    phi[nt, ] <- phi[nt - 1, ]
  }

  list(phi = phi, times = sim_times)
}


#' Compute the Jacobian matrix for 2D transient hydraulic tomography
#' via the adjoint-state method
#'
#' For each pumping test \eqn{j} and each observation \eqn{i} (at time
#' \eqn{t_{\text{obs}}}), computes the sensitivity of the observed head
#' to log-transmissivity \eqn{\ln T} using the adjoint-state method of
#' Zha et al. (2020):
#' \deqn{\frac{\partial H(\mathbf{x}_o, t)}{\partial \ln T(\mathbf{x}_Y)} =
#'       T(\mathbf{x}_Y) \int_0^{t}
#'       \left[\nabla \varphi(\mathbf{x}_Y, t-\tau) \cdot
#'             \nabla H(\mathbf{x}_Y, \tau)\right] \, d\tau \, \Delta x \Delta y}
#'
#' **Time discretisation**: Forward and adjoint simulations use independent
#' geometric time grids (fine at start, coarser later).  The convolution
#' interpolates both onto a common uniform grid via \code{stats::approx()}
#' before integration.  Observation times only determine the upper limit
#' of the convolution integral — they do not constrain the simulation grids.
#'
#' The key efficiency is that the adjoint \eqn{\varphi} for a given observation
#' location is computed once and shared across all observation times at that
#' location.
#'
#' @param grid  grid list from \code{GenGrid()}.
#' @param TT    transmissivity field \eqn{T} \code{[L^2/T]}, length-\code{n} vector.
#' @param SS    storage coefficient \eqn{S} [-] (default \code{1e-4}).
#' @param qHT   list of pumping test data frames (columns \code{Qp, x, y}).
#' @param oHT   list of observation data frames (columns \code{data, x, y, time}).
#' @param times numeric vector of output times for the ODE solver.
#'   If \code{NULL} (default), a geometric grid is auto-generated.
#' @param t_max end time of the simulation.  If \code{NULL} (default),
#'   set to the maximum observation time across all tests.
#' @param npoints number of time points in auto-generated grids (default 151).
#' @param p_grid growth factor for auto-generated time grids (default 1.03).
#' @param nconv  number of points on the common convolution grid (default 201).
#' @param lrw   real work array length (default 1600000).
#' @return A matrix of dimension \code{nobs x nelem}.
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' grid <- GenGrid(c(20,20,0,20,0,20))
#' TT <- random2d(nsim=1, grid=grid)$Tp
#' qHT <- list(data.frame(Qp=10, x=10.5, y=10.5))
#' oHT <- list(data.frame(data=-1, x=5.5, y=5.5, time=50))
#' J <- jacobian2DTr(grid=grid, TT=TT, SS=1e-4, qHT=qHT, oHT=oHT)
#' }
jacobian2DTr <- function(grid,
                         TT,
                         SS      = 1e-4,
                         qHT     = list(data.frame(Qp = 10, x = 20.5, y = 20.5)),
                         oHT     = list(data.frame(data = -1, x = 11, y = 11,
                                                  time = 50)),
                         times   = NULL,
                         t_max   = NULL,
                         npoints = 151,
                         p_grid  = 1.03,
                         nconv   = 201,
                         lrw     = 1600000) {

  require('deSolve')

  nHT <- length(qHT)
  n   <- grid$n
  nx  <- grid$nx
  ny  <- grid$ny
  dx  <- grid$dx
  dy  <- grid$dy
  dV  <- dx * dy

  Tp <- if (length(TT) == 1) rep(TT, n) else TT
  Sp <- if (length(SS) == 1) rep(SS, n) else SS

  if (length(Tp) != n) {
    stop("length(TT) = ", length(TT), " does not match grid$n = ", n,
         ". Use TT of length ", n, " (or a scalar).")
  }
  if (length(Sp) != n) {
    stop("length(SS) = ", length(SS), " does not match grid$n = ", n,
         ". Use SS of length ", n, " (or a scalar).")
  }

  # ---- Collect observation element indices and times ----
  loc_obsHT <- list()
  nobs_per_test <- integer(nHT)
  obs_times_HT <- list()
  for (j in seq_len(nHT)) {
    oinf <- oHT[[j]]
    oelemdf <- getOelem(grid = grid, Oinf = oinf)
    loc_obsHT[[j]] <- oelemdf$nelem
    nobs_per_test[j] <- nrow(oinf)
    obs_times_HT[[j]] <- oinf$time
  }
  nobs <- sum(nobs_per_test)

  # ---- Determine simulation end time ----
  all_obs_times <- unique(unlist(obs_times_HT))
  if (length(all_obs_times) == 0) all_obs_times <- 0
  if (is.null(t_max)) {
    t_max <- max(all_obs_times)
  }
  if (t_max <= 0) t_max <- 1

  # ---- Build simulation time grids (independent of obs times) ----
  if (is.null(times)) {
    sim_times <- .makeTimeGrid(t_max, npoints, p_grid)
  } else {
    # Merge user-supplied times with a geometric grid to ensure sufficient
    # resolution at early times where the adjoint impulse response varies
    # rapidly.  Without this, sparse uniform grids (e.g. seq(0,50,1)) can
    # severely degrade the accuracy of the convolution integral.
    auto_times <- .makeTimeGrid(t_max, npoints, p_grid)
    sim_times <- sort(unique(c(times, auto_times)))
    if (length(sim_times) < 2) {
      sim_times <- seq(0, max(sim_times, t_max, 1), length.out = npoints)
    }
  }

  # ---- Helper: compute cell-centred 2D gradients ----
  grad2D <- function(v) {
    mat <- matrix(v, nx, ny)
    gx <- matrix(0, nx, ny)
    gy <- matrix(0, nx, ny)
    if (nx > 1) {
      gx[2:(nx - 1), ] <- (mat[3:nx, ] - mat[1:(nx - 2), ]) / (2 * dx)
      gx[1, ]  <- (mat[2, ] - mat[1, ]) / dx
      gx[nx, ] <- (mat[nx, ] - mat[nx - 1, ]) / dx
    }
    if (ny > 1) {
      gy[, 2:(ny - 1)] <- (mat[, 3:ny] - mat[, 1:(ny - 2)]) / (2 * dy)
      gy[, 1]  <- (mat[, 2] - mat[, 1]) / dy
      gy[, ny] <- (mat[, ny] - mat[, ny - 1]) / dy
    }
    list(x = as.vector(gx), y = as.vector(gy))
  }

  # ---- Pre-allocate Jacobian ----
  J <- matrix(0, nrow = nobs, ncol = n)

  # ---- Forward simulations for each pumping test ----
  fwd_cache <- list()
  for (j in seq_len(nHT)) {
    res_fwd <- Ftransient2dsim(grid = grid, TT = Tp, SS = Sp,
                                Qinf = qHT[[j]], times = sim_times, lrw = lrw)
    out_fwd <- res_fwd$out
    fwd_t <- out_fwd[, 1]
    h_mat  <- out_fwd[, -1, drop = FALSE]

    nt <- nrow(h_mat)
    grad_h_x <- matrix(0, nt, n)
    grad_h_y <- matrix(0, nt, n)
    for (it in seq_len(nt)) {
      g <- grad2D(h_mat[it, ])
      grad_h_x[it, ] <- g$x
      grad_h_y[it, ] <- g$y
    }

    fwd_cache[[j]] <- list(
      times = fwd_t, h = h_mat,
      grad_x = grad_h_x, grad_y = grad_h_y
    )
  }

  # ---- Adjoint (Green's function) for each unique observation location ----
  all_obs_locs <- unique(unlist(loc_obsHT))
  adj_cache <- list()
  for (iobs in all_obs_locs) {
    adj_res <- adjoint2DTr(grid = grid, TT = Tp, SS = Sp,
                            iobs = iobs, times = sim_times, lrw = lrw)
    phi       <- adj_res$phi
    adj_t     <- adj_res$times

    nt <- nrow(phi)
    grad_phi_x <- matrix(0, nt, n)
    grad_phi_y <- matrix(0, nt, n)
    for (it in seq_len(nt)) {
      g <- grad2D(phi[it, ])
      grad_phi_x[it, ] <- g$x
      grad_phi_y[it, ] <- g$y
    }

    adj_cache[[as.character(iobs)]] <- list(
      times = adj_t,
      grad_x = grad_phi_x, grad_y = grad_phi_y
    )
  }

  # ---- Convolution helper: interpolate grad vectors to a common time grid ----
  # Returns a matrix nt_common x n of interpolated values.
  .interpGrad <- function(t_old, grad_mat, t_new) {
    nt_new <- length(t_new)
    ncell  <- ncol(grad_mat)
    out <- matrix(0, nt_new, ncell)
    for (ic in seq_len(ncell)) {
      out[, ic] <- approx(t_old, grad_mat[, ic], xout = t_new, rule = 2)$y
    }
    out
  }

  # ---- Compute sensitivity by time convolution ----
  row_offset <- 0
  for (j in seq_len(nHT)) {
    obs_elem <- loc_obsHT[[j]]
    obs_t    <- obs_times_HT[[j]]
    nlocal   <- length(obs_elem)

    fwd <- fwd_cache[[j]]
    fwd_times <- fwd$times

    for (i_local in seq_len(nlocal)) {
      iobs <- obs_elem[i_local]
      tobs_i <- obs_t[i_local]
      global_row <- row_offset + i_local

      if (tobs_i <= 0) next  # no sensitivity at t=0

      adj <- adj_cache[[as.character(iobs)]]

      # ---- Build common convolution time grid [0, tobs] ----
      t_conv <- seq(0, tobs_i, length.out = nconv)
      dt_conv <- t_conv[2] - t_conv[1]

      # Interpolate forward and adjoint gradients onto t_conv
      # Forward:  grad_h(tau)    for tau in [0, tobs]
      grad_hx <- .interpGrad(fwd_times, fwd$grad_x, t_conv)
      grad_hy <- .interpGrad(fwd_times, fwd$grad_y, t_conv)

      # Adjoint:  grad_phi(tobs - tau)  for tau in [0, tobs]
      t_adj_target <- tobs_i - t_conv  # tobs, tobs-dt, ..., 0
      grad_px <- .interpGrad(adj$times, adj$grad_x, t_adj_target)
      grad_py <- .interpGrad(adj$times, adj$grad_y, t_adj_target)


      # ---- Trapezoidal integration over t_conv ----
      integrand <- grad_hx * grad_px + grad_hy * grad_py
      # Trapezoidal rule: sum_{i} (integrand[i] + integrand[i+1])/2 * dt
      ntc <- nrow(integrand)
      if (ntc > 1) {
        S_int <- dt_conv * (colSums(integrand[-ntc, , drop = FALSE]) +
                            colSums(integrand[-1,   , drop = FALSE])) / 2
      } else {
        S_int <- rep(0, n)
      }

      J[global_row, ] <- -Tp * S_int * dV
    }

    row_offset <- row_offset + nlocal
  }

  return(J)
}



#' 2D Transient adjoint-based deterministic inversion for hydraulic tomography
#'
#' Transient counterpart of \code{\link{FinverseAdj}}.  Performs deterministic
#' Bayesian inversion to estimate the 2D log-transmissivity field
#' \eqn{\ln T(\mathbf{x})} from transient hydraulic tomography (HT) data using
#' the **adjoint-state method** to compute the Jacobian (sensitivity) matrix
#' analytically via \code{\link{jacobian2DTr}}.
#'
#' The update follows the same quasi-linear Bayesian formulation as
#' \code{\link{FinverseAdj}}, with forward simulations handled by
#' \code{\link{Ftransient2dsim}} and observations sampled via
#' \code{\link{samDataTr}}.
#'
#' @param domain   a 6-element vector \code{c(nx, ny, x1, x2, y1, y2)} \code{[L]}.
#' @param grid     grid list from \code{GenGrid()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT      list of pumping test data frames (columns \code{Qp, x, y}).
#' @param itermax  maximum number of iterations (default 5).
#' @param rmsemin  minimum RMSE stopping criterion (default 0).
#' @param oHT      list of observation data frames (columns \code{data, x, y, time}).
#' @param SS       storage coefficient \eqn{S} [-] — scalar or length-\code{n}
#'   vector (default \code{1e-4}).
#' @param times    numeric vector of output times for the ODE solver.
#'   If \code{NULL} (default), a geometric grid is auto-generated.
#' @param npoints  number of time points for auto-generated grids (default 151).
#' @param p_grid   growth factor for auto-generated time grids (default 1.03).
#' @param nconv    number of points on the common convolution grid (default 201).
#' @param lrw      real work array length for \code{ode.2D} (default 1600000).
#' @param sigma2obs observation error variance (default 1e-4).
#' @param lambda_lm Levenberg-Marquardt damping parameter (default 1.0); decays
#'   by \code{decay_lm} each iteration.
#' @param decay_lm  decay factor for \code{lambda_lm} (default 2.0).
#' @param geo      prior geostatistical parameters:
#'   \code{list(me, var, geomod, anis, range, nugget)}.
#' @param m0       initial log-T field, length \code{n} vector.  If \code{NULL}
#'   (default), uses the prior mean \code{geo$me} everywhere.
#' @param ifcor    logical; if \code{TRUE}, returns the Jacobian matrix and
#'   forward response after the first iteration without performing the update.
#' @param ifvarTh  logical; if \code{TRUE}, computes and stores \code{varobsh}
#'   and \code{varT}.  Default \code{FALSE}.
#' @return A list of per-iteration results (same structure as \code{\link{FinverseAdj}}).
#' @export
#' @examples
#' \dontrun{
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT <- trueK$Tp
#' grid <- GenGrid(c(40,40,0,40,0,40))
#' Qinf1 <- data.frame(Qp=10, x=20.5, y=20.5)
#' qHT <- list(test1 = Qinf1)
#' loc <- expand.grid(x=c(15,18,22,25,30), y=c(15,18,22,25,30))
#' Oinf1 <- data.frame(data=NA, x=loc$x, y=loc$y, time=50)
#' res_fwd <- Ftransient2dsim(grid=grid, TT=TT, Qinf=Qinf1, times=c(0,50))
#' Oinf1 <- samDataTr(Oinf=Oinf1, grid=grid, result_tr=res_fwd)
#' oHT <- list(test1 = Oinf1)
#' result <- FinverseAdjTr(grid=grid, qHT=qHT, oHT=oHT, itermax=3)
#' }
FinverseAdjTr <- function(
    domain      = c(40, 40, 0, 40, 0, 40),
    grid        = NULL,
    qHT         = list(data.frame(Qp = 10, x = 20.5, y = 20.5)),
    itermax     = 5,
    rmsemin     = 0,
    oHT         = list(data.frame(data = -1, x = 11, y = 11, time = 50)),
    SS          = 1e-4,
    times       = NULL,
    npoints     = 151,
    p_grid      = 1.03,
    nconv       = 201,
    lrw         = 1600000,
    sigma2obs   = 1e-4,
    lambda_lm   = 1.0,
    decay_lm    = 2.0,
    geo         = list(me = 0, var = 1, geomod = "Exp",
                       anis = c(90, 1), range = 30, nugget = 0),
    m0          = NULL,
    ifcor       = FALSE,
    ifvarTh     = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid(domain)

  n <- grid$n
  nHT <- length(qHT)

  # ---- Extract observation data ----
  trueobshHT <- list()
  loc_obsHT  <- list()
  for (i in seq_len(nHT)) {
    oinf <- oHT[[i]]
    oelemdf <- getOelem(grid = grid, Oinf = oinf)
    loc_obsHT[[i]] <- oelemdf$nelem
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs <- length(trueobsh)

  # ---- Build simulation time grid (independent of obs times) ----
  all_obs_times <- unique(unlist(lapply(oHT, function(x) x$time)))
  if (length(all_obs_times) == 0) all_obs_times <- 0
  t_max <- max(all_obs_times)

  if (is.null(times)) {
    sim_times <- .makeTimeGrid(t_max, npoints, p_grid)
  } else {
    sim_times <- times
  }

  # ---- Build prior covariance C_kk from geostatistical model ----
  require('gstat')
  xy <- grid$grid

  dist_mat <- as.matrix(dist(xy))
  if (geo$geomod == "Exp") {
    C_kk <- geo$var * exp(-dist_mat / geo$range)
    diag(C_kk) <- geo$var + geo$nugget
  } else if (geo$geomod == "Gau") {
    C_kk <- geo$var * exp(-(dist_mat / geo$range)^2)
    diag(C_kk) <- geo$var + geo$nugget
  } else if (geo$geomod == "Sph") {
    hr <- dist_mat / geo$range
    C_kk <- geo$var * (1 - 1.5 * hr + 0.5 * hr^3)
    C_kk[hr > 1] <- 0
    diag(C_kk) <- geo$var + geo$nugget
  } else {
    C_kk <- geo$var * exp(-dist_mat / geo$range)
    diag(C_kk) <- geo$var + geo$nugget
  }

  # ---- Initial parameter estimate ----
  if (is.null(m0)) {
    m <- rep(geo$me, n)  # log-T
  } else {
    m <- m0
  }
  T_current <- exp(m)

  # ---- Iteration loop ----
  niter <- 1
  rmse <- 1e10
  iterdf <- list()

  while (niter <= itermax && rmse > rmsemin) {

    # --- Forward simulation at current T ---
    h_sim_vec <- samHTtr(grid = grid, TT = T_current, SS = SS,
                         qHT = qHT, oHT = oHT, times = sim_times, lrw = lrw)

    # --- Misfit ---
    residual <- trueobsh - h_sim_vec
    rmse <- sqrt(mean(residual^2))
    l2 <- rmse
    l1 <- mean(abs(residual))

    msg <- paste('niter =', niter,
                 'rmse =', round(rmse, 6),
                 'l2 =', round(l2, 6),
                 'l1 =', round(l1, 6),
                 'lambda_lm =', round(lambda_lm, 4))
    print(msg)

    # --- Compute Jacobian ---
    J <- jacobian2DTr(grid = grid, TT = T_current, SS = SS,
                      qHT = qHT, oHT = oHT,
                      times = sim_times, nconv = nconv, lrw = lrw)

    if (ifcor) {
      covhk <- J %*% C_kk
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }

    # --- Compute update (dual / observation-space formulation) ---
    Jt <- t(J)
    JCJt <- J %*% C_kk %*% Jt
    S_mat <- JCJt + (sigma2obs + lambda_lm) * diag(nobs)

    beta <- solve(S_mat, residual)
    dm <- as.vector(C_kk %*% Jt %*% beta)

    # --- Estimate variances ---
    meanobsh <- h_sim_vec

    if (ifvarTh) {
      KS <- C_kk %*% Jt %*% solve(S_mat)
      C_post <- C_kk - KS %*% J %*% C_kk
      varT <- as.vector(diag(C_post))
      JCJt_post <- J %*% C_post %*% Jt
      varobsh <- as.vector(diag(JCJt_post))
    } else {
      varT <- NULL
      varobsh <- NULL
    }

    # --- Update parameters ---
    m <- m + dm
    T_current <- exp(m)

    # --- Store iteration results ---
    iterdf[[niter]] <- list(
      meanT    = as.vector(m),
      varT     = varT,
      meanobsh = meanobsh,
      varobsh  = varobsh,
      m        = m,
      h_sim    = h_sim_vec,
      rmse     = rmse,
      l2       = l2,
      l1       = l1
    )

    if (niter == 1) {
      iterdf[[niter]]$J <- J
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    lambda_lm <- lambda_lm / decay_lm
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}

# ==============================================================================
# 3D Adjoint-based Deterministic Inversion
# ==============================================================================

#' 3D Adjoint-based deterministic inversion for hydraulic tomography
#'
#' 3D counterpart of \code{\link{FinverseAdj}}.  Performs deterministic Bayesian
#' inversion to estimate the 3D log-hydraulic-conductivity field
#' \eqn{\ln K(\mathbf{x})} from hydraulic tomography (HT) data using the
#' **adjoint-state method** to compute the Jacobian (sensitivity) matrix
#' analytically via \code{\link{jacobian3D}}.
#'
#' The update follows the same quasi-linear Bayesian formulation as
#' \code{\link{FinverseAdj}}, replacing the 2D transmissivity with 3D hydraulic
#' conductivity.  A Levenberg-Marquardt damping parameter \code{lambda_lm}
#' is added to the diagonal of the Hessian approximation for stability.
#'
#' @param domain   9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2, z1, z2)}.
#' @param grid     grid list from \code{GenGrid3D()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT      list of pumping test data frames (columns \code{Qp, x, y, z}).
#' @param itermax  maximum number of iterations (default 5).
#' @param rmsemin  minimum RMSE stopping criterion (default 0).
#' @param oHT      list of observation data frames (columns \code{data, x, y, z}).
#' @param lrw      real work array length for \code{steady.3D} (default 20000000).
#' @param sigma2obs observation error variance (default 1e-4).
#' @param lambda_lm Levenberg-Marquardt damping parameter (default 1.0); decays
#'   by \code{decay_lm} each iteration.
#' @param decay_lm  decay factor for \code{lambda_lm} (default 2.0).
#' @param geo      prior geostatistical parameters:
#'   \code{list(me, var, geomod, range, nugget, anis)}.
#'   Note \code{anis} is a 5-element vector for 3D.
#' @param m0       initial log-K field, length \code{n} vector.  If \code{NULL}
#'   (default), uses the prior mean \code{geo$me} everywhere.
#' @param ifcor    logical; if \code{TRUE}, returns the Jacobian matrix and
#'   forward response after the first iteration without performing the update.
#' @param ifvarTh  logical; if \code{TRUE}, computes and stores the predicted
#'   head variance \code{varobsh} and the conditional (posterior) variance of
#'   log-K \code{varT} at each iteration.  Default \code{FALSE} because the
#'   covariance matrix operation can be expensive for large grids.
#' @return A list of per-iteration results, each a list with
#'   \code{meanT} (log-K estimate, same as \code{m}), \code{varT}
#'   (conditional variance of log-K, only when \code{ifvarTh = TRUE}),
#'   \code{meanobsh} (simulated heads), \code{varobsh} (predicted head variance,
#'   only when \code{ifvarTh = TRUE}),
#'   \code{m}, \code{h_sim}, \code{rmse}, \code{l2}, \code{l1}, and
#'   \code{J} (Jacobian, first iteration only).  When \code{ifcor = TRUE},
#'   returns \code{list(J = J, h_sim = h_sim)}.
#' @export
#' @examples
#' \dontrun{
#' domain3d <- c(15, 15, 5, 0, 15, 0, 15, 0, 5)
#' grid3d   <- GenGrid3D(domain3d)
#' set.seed(42)
#' trueK3d  <- random3d(nsim=1, grid=grid3d)
#' Qinf3d   <- data.frame(Qp=10, x=7.5, y=7.5, z=2.5)
#' qHT3d    <- list(test1 = Qinf3d)
#' trueh3d  <- Fsteady3dsim(grid=grid3d, KK=trueK3d$Kp, Qinf=Qinf3d)
#' loc      <- expand.grid(x=c(3,6,9,12), y=c(3,6,9,12))
#' Oinf3d   <- data.frame(data=NA, x=loc$x, y=loc$y, z=2.5)
#' Oinf3d   <- samData3D(Oinf=Oinf3d, grid=grid3d, h=trueh3d$solution)
#' oHT3d    <- list(test1 = Oinf3d)
#' result   <- Finverse3DAdj(grid=grid3d, qHT=qHT3d, oHT=oHT3d, itermax=3)
#' }
Finverse3DAdj <- function(
    domain      = c(15, 15, 5, 0, 15, 0, 15, 0, 5),
    grid        = NULL,
    qHT         = list(data.frame(Qp = 10, x = 7.5, y = 7.5, z = 2.5)),
    itermax     = 5,
    rmsemin     = 0,
    oHT         = list(data.frame(data = -1, x = 5, y = 5, z = 2.5)),
    lrw         = 20000000,
    sigma2obs   = 1e-4,
    lambda_lm   = 1.0,
    decay_lm    = 2.0,
    geo         = list(me = 0, var = 1, geomod = "Exp",
                       range = 10, nugget = 0,
                       anis = c(0, 0, 0, 1, 1)),
    m0          = NULL,
    ifcor       = FALSE,
    ifvarTh     = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid3D(domain)

  n <- grid$n
  nHT <- length(qHT)

  # ---- Extract observation data ----
  trueobshHT <- list()
  loc_obsHT  <- list()
  for (i in seq_len(nHT)) {
    oinf <- oHT[[i]]
    oelemdf <- getOelem3D(grid = grid, Oinf = oinf)
    loc_obsHT[[i]] <- oelemdf$nelem
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs <- length(trueobsh)

  # ---- Build prior covariance C_kk from geostatistical model (3D) ----
  require('gstat')
  xyz <- grid$grid

  dist_mat <- as.matrix(dist(xyz))
  if (geo$geomod == "Exp") {
    C_kk <- geo$var * exp(-dist_mat / geo$range)
    diag(C_kk) <- geo$var + geo$nugget
  } else if (geo$geomod == "Gau") {
    C_kk <- geo$var * exp(-(dist_mat / geo$range)^2)
    diag(C_kk) <- geo$var + geo$nugget
  } else if (geo$geomod == "Sph") {
    hr <- dist_mat / geo$range
    C_kk <- geo$var * (1 - 1.5 * hr + 0.5 * hr^3)
    C_kk[hr > 1] <- 0
    diag(C_kk) <- geo$var + geo$nugget
  } else {
    C_kk <- geo$var * exp(-dist_mat / geo$range)
    diag(C_kk) <- geo$var + geo$nugget
  }

  # ---- Initial parameter estimate ----
  if (is.null(m0)) {
    m <- rep(geo$me, n)  # log-K
  } else {
    m <- m0
  }
  K_current <- exp(m)

  # ---- Iteration loop ----
  niter <- 1
  rmse <- 1e10
  iterdf <- list()

  while (niter <= itermax && rmse > rmsemin) {

    # --- Forward simulation at current K ---
    h_sim_vec <- samHT3D(grid = grid, TT = K_current,
                         qHT = qHT, oHT = oHT, lrw = lrw)

    # --- Misfit ---
    residual <- trueobsh - h_sim_vec
    rmse <- sqrt(mean(residual^2))
    l2 <- rmse
    l1 <- mean(abs(residual))

    msg <- paste('niter =', niter,
                 'rmse =', round(rmse, 6),
                 'l2 =', round(l2, 6),
                 'l1 =', round(l1, 6),
                 'lambda_lm =', round(lambda_lm, 4))
    print(msg)

    # --- Compute Jacobian ---
    J <- jacobian3D(grid = grid, KK = K_current,
                    qHT = qHT, oHT = oHT, lrw = lrw)

    if (ifcor) {
      covhk <- J %*% C_kk
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }

    # --- Compute update (dual / observation-space formulation) ---
    Jt <- t(J)                                     # nelem x nobs
    JCJt <- J %*% C_kk %*% Jt                      # nobs x nobs
    S <- JCJt + (sigma2obs + lambda_lm) * diag(nobs)

    beta <- solve(S, residual)                     # nobs
    dm <- as.vector(C_kk %*% Jt %*% beta)          # nelem

    # --- Estimate variances ---
    meanobsh <- h_sim_vec

    if (ifvarTh) {
      # Predicted head variance = diag(J C_kk J^T)
      varobsh <- as.vector(diag(JCJt))

      # Conditional (posterior) parameter covariance
      KS <- C_kk %*% Jt %*% solve(S)               # nelem x nobs
      C_post <- C_kk - KS %*% J %*% C_kk           # nelem x nelem
      varT <- as.vector(diag(C_post))
    } else {
      varT <- NULL
      varobsh <- NULL
    }

    # --- Update parameters ---
    m <- m + dm
    K_current <- exp(m)

    # --- Store iteration results ---
    iterdf[[niter]] <- list(
      meanT    = as.vector(m),
      varT     = varT,
      meanobsh = meanobsh,
      varobsh  = varobsh,
      m        = m,
      h_sim    = h_sim_vec,
      rmse     = rmse,
      l2       = l2,
      l1       = l1
    )

    if (niter == 1) {
      iterdf[[niter]]$J <- J
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    lambda_lm <- lambda_lm / decay_lm
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}


#' 3D Adjoint-based deterministic inversion — well-screen version
#'
#' 3D counterpart of \code{\link{FinverseAdj}} for well-screen (vertical
#' interval) observations.  Uses \code{\link{jacobian3DScreen}} to compute
#' the Jacobian for interval-averaged head observations, where each
#' observation well has a screen defined by \code{z_top} and \code{z_bottom}.
#'
#' The inversion algorithm is otherwise identical to \code{\link{Finverse3DAdj}}.
#'
#' @param domain   9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2, z1, z2)}.
#' @param grid     grid list from \code{GenGrid3D()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT      list of pumping test data frames with columns
#'   \code{Qp, x, y, z_top, z_bottom} (well-screen pumping) or
#'   \code{Qp, x, y, z} (point pumping).
#' @param itermax  maximum number of iterations (default 5).
#' @param rmsemin  minimum RMSE stopping criterion (default 0).
#' @param oHT      list of observation data frames with columns
#'   \code{data, x, y, z_top, z_bottom}.
#' @param lrw      real work array length for \code{steady.3D} (default 20000000).
#' @param sigma2obs observation error variance (default 1e-4).
#' @param lambda_lm Levenberg-Marquardt damping parameter (default 1.0); decays
#'   by \code{decay_lm} each iteration.
#' @param decay_lm  decay factor for \code{lambda_lm} (default 2.0).
#' @param geo      prior geostatistical parameters:
#'   \code{list(me, var, geomod, range, nugget, anis)}.
#' @param m0       initial log-K field, length \code{n} vector.  If \code{NULL}
#'   (default), uses the prior mean \code{geo$me} everywhere.
#' @param ifcor    logical; if \code{TRUE}, returns the Jacobian matrix and
#'   forward response after the first iteration without performing the update.
#' @param ifvarTh  logical; if \code{TRUE}, computes and stores the predicted
#'   head variance \code{varobsh} and the conditional (posterior) variance of
#'   log-K \code{varT} at each iteration.  Default \code{FALSE}.
#' @return A list of per-iteration results, each a list with
#'   \code{meanT}, \code{varT}, \code{meanobsh}, \code{varobsh},
#'   \code{m}, \code{h_sim}, \code{rmse}, \code{l2}, \code{l1}, and
#'   \code{J} (Jacobian, first iteration only).  When \code{ifcor = TRUE},
#'   returns \code{list(J = J, h_sim = h_sim)}.
#' @export
#' @examples
#' \dontrun{
#' domain3d <- c(15, 15, 5, 0, 15, 0, 15, 0, 5)
#' grid3d   <- GenGrid3D(domain3d)
#' set.seed(42)
#' trueK3d  <- random3d(nsim=1, grid=grid3d)
#' Qinf3d   <- data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)
#' qHT3d    <- list(test1 = Qinf3d)
#' res3d    <- Fsteady3dsim(grid=grid3d, KK=trueK3d$Kp, Qinf=Qinf3d)
#' loc      <- expand.grid(x=c(3,6,9,12), y=c(3,6,9,12))
#' Oinf3d   <- data.frame(data=NA, x=loc$x, y=loc$y,
#'                        z_top=1, z_bottom=4)
#' Oinf3d   <- samData3DScreen(Oinf=Oinf3d, grid=grid3d, h=res3d$solution)
#' oHT3d    <- list(test1 = Oinf3d)
#' result   <- Finverse3DScreenAdj(grid=grid3d, qHT=qHT3d, oHT=oHT3d, itermax=3)
#' }
Finverse3DScreenAdj <- function(
    domain      = c(15, 15, 5, 0, 15, 0, 15, 0, 5),
    grid        = NULL,
    qHT         = list(data.frame(Qp = 10, x = 7.5, y = 7.5,
                                  z_top = 2, z_bottom = 3)),
    itermax     = 5,
    rmsemin     = 0,
    oHT         = list(data.frame(data = -1, x = 5, y = 5,
                                  z_top = 1, z_bottom = 4)),
    lrw         = 20000000,
    sigma2obs   = 1e-4,
    lambda_lm   = 1.0,
    decay_lm    = 2.0,
    geo         = list(me = 0, var = 1, geomod = "Exp",
                       range = 10, nugget = 0,
                       anis = c(0, 0, 0, 1, 1)),
    m0          = NULL,
    ifcor       = FALSE,
    ifvarTh     = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid3D(domain)

  n <- grid$n
  nHT <- length(qHT)

  # ---- Extract observation data ----
  trueobshHT <- list()
  for (i in seq_len(nHT)) {
    oinf <- oHT[[i]]
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs <- length(trueobsh)

  # ---- Build prior covariance C_kk from geostatistical model (3D) ----
  require('gstat')
  xyz <- grid$grid

  dist_mat <- as.matrix(dist(xyz))
  if (geo$geomod == "Exp") {
    C_kk <- geo$var * exp(-dist_mat / geo$range)
    diag(C_kk) <- geo$var + geo$nugget
  } else if (geo$geomod == "Gau") {
    C_kk <- geo$var * exp(-(dist_mat / geo$range)^2)
    diag(C_kk) <- geo$var + geo$nugget
  } else if (geo$geomod == "Sph") {
    hr <- dist_mat / geo$range
    C_kk <- geo$var * (1 - 1.5 * hr + 0.5 * hr^3)
    C_kk[hr > 1] <- 0
    diag(C_kk) <- geo$var + geo$nugget
  } else {
    C_kk <- geo$var * exp(-dist_mat / geo$range)
    diag(C_kk) <- geo$var + geo$nugget
  }

  # ---- Initial parameter estimate ----
  if (is.null(m0)) {
    m <- rep(geo$me, n)  # log-K
  } else {
    m <- m0
  }
  K_current <- exp(m)

  # ---- Iteration loop ----
  niter <- 1
  rmse <- 1e10
  iterdf <- list()

  while (niter <= itermax && rmse > rmsemin) {

    # --- Forward simulation at current K (screen-averaged) ---
    h_sim_vec <- samHT3DScreen(grid = grid, TT = K_current,
                               qHT = qHT, oHT = oHT, lrw = lrw)

    # --- Misfit ---
    residual <- trueobsh - h_sim_vec
    rmse <- sqrt(mean(residual^2))
    l2 <- rmse
    l1 <- mean(abs(residual))

    msg <- paste('niter =', niter,
                 'rmse =', round(rmse, 6),
                 'l2 =', round(l2, 6),
                 'l1 =', round(l1, 6),
                 'lambda_lm =', round(lambda_lm, 4))
    print(msg)

    # --- Compute Jacobian (screen-averaged) ---
    J <- jacobian3DScreen(grid = grid, KK = K_current,
                          qHT = qHT, oHT = oHT, lrw = lrw)

    if (ifcor) {
      covhk <- J %*% C_kk
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }

    # --- Compute update (dual / observation-space formulation) ---
    Jt <- t(J)                                     # nelem x nobs
    JCJt <- J %*% C_kk %*% Jt                      # nobs x nobs
    S <- JCJt + (sigma2obs + lambda_lm) * diag(nobs)

    beta <- solve(S, residual)                     # nobs
    dm <- as.vector(C_kk %*% Jt %*% beta)          # nelem

    # --- Estimate variances ---
    meanobsh <- h_sim_vec

    if (ifvarTh) {
      varobsh <- as.vector(diag(JCJt))

      KS <- C_kk %*% Jt %*% solve(S)
      C_post <- C_kk - KS %*% J %*% C_kk
      varT <- as.vector(diag(C_post))
    } else {
      varT <- NULL
      varobsh <- NULL
    }

    # --- Update parameters ---
    m <- m + dm
    K_current <- exp(m)

    # --- Store iteration results ---
    iterdf[[niter]] <- list(
      meanT    = as.vector(m),
      varT     = varT,
      meanobsh = meanobsh,
      varobsh  = varobsh,
      m        = m,
      h_sim    = h_sim_vec,
      rmse     = rmse,
      l2       = l2,
      l1       = l1
    )

    if (niter == 1) {
      iterdf[[niter]]$J <- J
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    lambda_lm <- lambda_lm / decay_lm
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}


# ==============================================================================
# 2D Transient Hydraulic Tomography — Inverse Functions
# ==============================================================================

#' Run forward simulations for one T field across all 2D transient pumping tests
#'
#' Transient counterpart of \code{samHT}. Calls \code{Ftransient2dsim} for each
#' pumping test and samples simulated heads at observation well locations and
#' times using \code{samDataTr}.
#'
#' @param grid   grid from \code{GenGrid()}.
#' @param TT     length-\code{n} vector of transmissivity T \code{[L^2/T]}.
#' @param SS     storage coefficient S [-] — scalar or length-\code{n} vector.
#' @param qHT    list of pumping test data frames, each with columns
#'   \code{Qp, x, y}.
#' @param oHT    list of observation data frames, each with columns
#'   \code{data, x, y, time}.
#' @param times  numeric vector of output times for the ODE solver.  Must
#'   include all observation times in \code{oHT}.
#' @param lrw    real work array length for \code{ode.2D}; default 1600000.
#' @param simplify if \code{TRUE} return a plain vector of sampled heads
#'   (all tests concatenated); if \code{FALSE} return updated \code{oHT}.
#' @return vector of sampled drawdown values (when \code{simplify = TRUE}).
#' @export
#' @examples
#' grid  <- GenGrid()
#' Qinf  <- data.frame(Qp=10, x=20.5, y=20.5)
#' Oinf  <- data.frame(data=NA, x=c(15,25), y=c(20,20), time=c(0.5, 0.5))
#' times <- seq(0, 1, by=0.05)
#' da    <- samHTtr(grid=grid, TT=0.1, SS=1e-4,
#'                  qHT=list(Qinf), oHT=list(Oinf), times=times)
samHTtr <- function(grid,
                    TT       = 0.1,
                    SS       = 1e-4,
                    qHT      = list(data.frame(Qp=10, x=20.5, y=20.5)),
                    oHT      = list(data.frame(data=NA, x=11, y=11, time=50)),
                    times    = seq(0, 100, by=10),
                    lrw      = 1600000,
                    simplify = TRUE) {

  nHT <- length(qHT)
  for (i in seq_len(nHT)) {
    qinf <- qHT[[i]]
    oinf <- oHT[[i]]
    res  <- Ftransient2dsim(grid = grid, TT = TT, SS = SS,
                            Qinf = qinf, times = times, lrw = lrw)
    oinf <- samDataTr(Oinf = oinf, grid = grid, result_tr = res)
    oHT[[i]] <- oinf[, c('data', 'x', 'y', 'time')]
  }
  if (simplify) {
    oHTdf <- dplyr::bind_rows(oHT, .id = 'id')
    return(oHTdf$data)
  } else {
    return(oHT)
  }
}


#' Parallel Monte Carlo forward runs for 2D transient HT (ensemble)
#'
#' Transient counterpart of \code{samHTmcPar}. Runs \code{samHTtr} in parallel
#' across all ensemble members using \pkg{foreach} + \pkg{doParallel}.
#'
#' @param grid   grid from \code{GenGrid()}.
#' @param TT     \code{n × nsim} matrix of T realisations.
#' @param SS     storage coefficient S [-] — scalar or length-\code{n} vector.
#' @param qHT    list of pumping test data frames (columns \code{Qp, x, y}).
#' @param oHT    list of observation data frames (columns \code{data, x, y, time}).
#' @param times  numeric vector of output times for the ODE solver.
#' @param lrw    real work array length; default 1600000.
#' @param ncore  number of parallel cores.
#' @return \code{nsim × nobs} matrix of simulated heads at observation wells.
#' @export
#' @examples
#' grid  <- GenGrid()
#' TTmat <- random2d(nsim=5, grid=grid)
#' TTmat <- as.matrix(TTmat[,-c(1,2)])
#' Qinf  <- data.frame(Qp=10, x=20.5, y=20.5)
#' Oinf  <- data.frame(data=NA, x=c(15,25), y=c(20,20), time=c(0.5, 0.5))
#' times <- seq(0, 1, by=0.05)
#' da    <- samHTmcParTr(grid=grid, TT=TTmat, SS=1e-4,
#'                       qHT=list(Qinf), oHT=list(Oinf), times=times, ncore=2)
samHTmcParTr <- function(grid,
                         TT    = 0.1,
                         SS    = 1e-4,
                         qHT   = list(data.frame(Qp=10, x=20.5, y=20.5)),
                         oHT   = list(data.frame(data=NA, x=11, y=11, time=50)),
                         times = seq(0, 100, by=10),
                         lrw   = 1600000,
                         ncore = 4) {

  library('foreach')
  library('doParallel')
  registerDoParallel(cores = ncore)

  nsim <- ncol(TT)
  x <- foreach(i = 1:nsim, .combine = rbind) %dopar% {
    HydroTomo::samHTtr(grid = grid, TT = TT[, i], SS = SS,
                       qHT = qHT, oHT = oHT, times = times, lrw = lrw)
  }

  stopImplicitCluster()
  return(x)
}


#' 2D Transient Ensemble-based Bayesian inverse for hydraulic tomography (parallel)
#'
#' Transient extension of \code{Finverse3}.  Replaces the steady-state forward
#' solver with \code{Ftransient2dsim} and uses \code{samDataTr} for time-aware
#' observation sampling.  The Bayesian updating equations (ensemble smoother)
#' are identical to the steady-state version.
#'
#' @param domain  6-element vector \code{c(nx, ny, x1, x2, y1, y2)}.
#' @param grid    grid from \code{GenGrid()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT     list of pumping test data frames with columns
#'   \code{Qp, x, y}.
#' @param nsim    ensemble size (default 50).
#' @param itermax maximum iterations.
#' @param varmeanTmax  convergence threshold on variance of mean ln(T).
#' @param rmsemin minimum RMSE stopping criterion.
#' @param mul     stabiliser (default 1.0).
#' @param decay   stabiliser decay per iteration (default 1.05).
#' @param oHT     list of observation data frames with columns
#'   \code{data, x, y, time}.
#' @param SS      storage coefficient S [-] — scalar or length-\code{n}
#'   vector.  Default \code{1e-4}.
#' @param times   numeric vector of output times for the ODE solver.  Must
#'   cover all observation times in \code{oHT}.
#' @param lrw     real work array length for \code{ode.2D}
#'   (default 1600000; increase for larger grids).
#' @param ncore   number of parallel cores (default 10).
#' @param geo     prior geostatistical parameters — same list structure as
#'   \code{random2d()}: \code{list(me, var, geomod, anis, range, nugget)}.
#' @return list of per-iteration results, each a list with
#'   \code{meanT} (mean ln T vector), \code{varT} (variance of ln T),
#'   \code{meanobsh}, \code{varobsh}.
#' @export
#' @examples
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT    <- trueK[,-c(1,2)]
#' domain <- c(40,40,0,40,0,40)
#' grid   <- GenGrid(domain)
#' Qinf1  <- data.frame(Qp=10, x=20.5, y=20.5)
#' qHT    <- list(test1 = Qinf1)
#' times  <- seq(0, 1, by=0.05)
#' res    <- Ftransient2dsim(grid=grid, TT=TT, SS=1e-4,
#'                           Qinf=Qinf1, times=times)
#' locx   <- c(15,18,22,25,30)
#' locy   <- c(15,18,22,25,30)
#' loc    <- expand.grid(x=locx, y=locy)
#' Oinf1  <- data.frame(data=NA, x=loc$x, y=loc$y, time=0.5)
#' Oinf1  <- samDataTr(Oinf=Oinf1, grid=grid, result_tr=res)
#' oHT    <- list(test1 = Oinf1)
#' result <- Finverse3Tr(grid=grid, qHT=qHT, oHT=oHT,
#'                       SS=1e-4, times=times, nsim=20, itermax=1, ncore=2)
Finverse3Tr <- function(
    domain      = c(40, 40, 0, 40, 0, 40),
    grid        = NULL,
    qHT         = list(data.frame(Qp=10, x=20.5, y=20.5)),
    nsim        = 50,
    itermax     = 5,
    varmeanTmax = 5,
    rmsemin     = 0,
    mul         = 1.0,
    decay       = 1.05,
    oHT         = list(data.frame(data=-1, x=11, y=11, time=50)),
    SS          = 1e-4,
    times       = seq(0, 100, by=10),
    lrw         = 1600000,
    ncore       = 10,
    geo         = list(me=0, var=1, geomod="Exp",
                       anis=c(90,1), range=30, nugget=0),
    ifcor       = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid(domain)

  nHT <- length(qHT)

  # ---- extract observation data ----------------------------------------------
  trueobshHT <- list()
  for (i in seq_len(nHT)) {
    oinf            <- oHT[[i]]
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs     <- length(trueobsh)

  # ---- initial ensemble ------------------------------------------------------
  yy   <- random2d(nsim = nsim, grid = grid, geo = geo)
  Tnew <- as.matrix(yy[, -c(1, 2)])   # n × nsim matrix of T values

  # ---- iteration loop --------------------------------------------------------
  niter    <- 1
  varmeanT <- 0
  rmse     <- 1e10
  msgdf    <- data.frame(niter = niter, varmeanT = varmeanT, rmse = rmse)
  iterdf   <- list()

  while (niter <= itermax & varmeanT < varmeanTmax & rmse > rmsemin) {

    varT     <- apply(log(Tnew), 1, var)
    meanT    <- apply(log(Tnew), 1, mean)
    varmeanT <- var(meanT)

    # parallel transient forward runs for all ensemble members
    obsh <- samHTmcParTr(grid  = grid,
                         TT    = Tnew,
                         SS    = SS,
                         qHT   = qHT,
                         oHT   = oHT,
                         times = times,
                         lrw   = lrw,
                         ncore = ncore)

    # ---- statistics ----------------------------------------------------------
    varobsh  <- apply(obsh, 2, var)
    meanobsh <- apply(obsh, 2, mean)
    weigs    <- 1 / varobsh / sum(1 / varobsh)
    rmse     <- mean((trueobsh - meanobsh)^2 * weigs)^0.5
    l2       <- mean((trueobsh - meanobsh)^2)^0.5
    l1       <- mean(abs(trueobsh - meanobsh))

    # ---- Bayesian update (same as steady-state) ------------------------------
    covh  <- cov(obsh)
    covhk <- cov(obsh, t(log(Tnew)))
    if (ifcor) {
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }
    covh1 <- covh
    diag(covh1) <- rep((1 + mul) * max(diag(covh)), nobs)
    a <- solve(covh1, covhk)   # nobs × n

    for (i in seq_len(nsim)) {
      Tnew[, i] <- Tnew[, i] * exp(t(a) %*% (trueobsh - obsh[i, ]))
    }

    msg <- paste('niter =', niter,
                 'varmeanT =', round(varmeanT, 4),
                 'rmse =', round(rmse, 4),
                 'l2 =', round(l2, 4),
                 'l1 =', round(l1, 4))
    print(msg)

    iterdf[[niter]] <- list(meanT    = as.vector(meanT),
                            varT     = as.vector(varT),
                            meanobsh = meanobsh,
                            varobsh  = varobsh)

    if (niter == 1) {
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    mul   <- mul / decay
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}


#' Run forward simulations for one K field across all 3D transient pumping tests
#'
#' Transient counterpart of \code{samHT3D}. Calls \code{Ftransient3dsim} for each
#' pumping test and samples simulated heads at observation well locations and
#' times using \code{samData3DTr}.
#'
#' @param grid   grid from \code{GenGrid3D()}.
#' @param KK     length-\code{n} vector of hydraulic conductivity K \code{[L/T]}.
#' @param Ss     specific storage Ss \code{[1/L]} — scalar or length-\code{n}
#'   vector.
#' @param qHT    list of pumping test data frames, each with columns
#'   \code{Qp, x, y, z}.
#' @param oHT    list of observation data frames, each with columns
#'   \code{data, x, y, z, time}.
#' @param times  numeric vector of output times for the ODE solver.  Must
#'   include all observation times in \code{oHT}.
#' @param lrw    real work array length for \code{ode.3D}; default 20000000.
#' @param simplify if \code{TRUE} return a plain vector of sampled heads
#'   (all tests concatenated); if \code{FALSE} return updated \code{oHT}.
#' @return vector of sampled drawdown values (when \code{simplify = TRUE}).
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15, 15, 5, 0, 15, 0, 15, 0, 5))
#' Qinf   <- data.frame(Qp=10, x=7.5, y=7.5, z=2.5)
#' Oinf   <- data.frame(data=NA, x=c(5,10), y=c(7.5,7.5),
#'                      z=c(2.5,2.5), time=c(0.5, 0.5))
#' times  <- seq(0, 1, by=0.05)
#' da     <- samHT3Dtr(grid=grid3d, KK=0.1, Ss=1e-4,
#'                     qHT=list(Qinf), oHT=list(Oinf), times=times)
samHT3Dtr <- function(grid,
                      KK       = 0.1,
                      Ss       = 1e-4,
                      qHT      = list(data.frame(Qp=10, x=7.5, y=7.5, z=2.5)),
                      oHT      = list(data.frame(data=NA, x=5, y=5, z=2.5, time=50)),
                      times    = seq(0, 100, by=10),
                      lrw      = 20000000,
                      simplify = TRUE) {

  nHT <- length(qHT)
  for (i in seq_len(nHT)) {
    qinf <- qHT[[i]]
    oinf <- oHT[[i]]
    res  <- Ftransient3dsim(grid = grid, KK = KK, Ss = Ss,
                            Qinf = qinf, times = times, lrw = lrw)
    oinf <- samData3DTr(Oinf = oinf, grid = grid, result_tr = res)
    oHT[[i]] <- oinf[, c('data', 'x', 'y', 'z', 'time')]
  }
  if (simplify) {
    oHTdf <- dplyr::bind_rows(oHT, .id = 'id')
    return(oHTdf$data)
  } else {
    return(oHT)
  }
}


#' Parallel Monte Carlo forward runs for 3D transient HT (ensemble)
#'
#' Transient counterpart of \code{samHTmcPar3D}. Runs \code{samHT3Dtr} in
#' parallel across all ensemble members using \pkg{foreach} + \pkg{doParallel}.
#'
#' @param grid   grid from \code{GenGrid3D()}.
#' @param KK     \code{n × nsim} matrix of K realisations.
#' @param Ss     specific storage Ss \code{[1/L]} — scalar or length-\code{n}
#'   vector.
#' @param qHT    list of pumping test data frames (columns \code{Qp, x, y, z}).
#' @param oHT    list of observation data frames (columns
#'   \code{data, x, y, z, time}).
#' @param times  numeric vector of output times for the ODE solver.
#' @param lrw    real work array length; default 20000000.
#' @param ncore  number of parallel cores.
#' @return \code{nsim × nobs} matrix of simulated heads at observation wells.
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15,15,5,0,15,0,15,0,5))
#' KKmat  <- random3d(nsim=5, grid=grid3d)
#' KKmat  <- as.matrix(KKmat[,-c(1,2,3)])
#' Qinf   <- data.frame(Qp=10, x=7.5, y=7.5, z=2.5)
#' Oinf   <- data.frame(data=NA, x=c(5,10), y=c(7.5,7.5),
#'                      z=c(2.5,2.5), time=c(0.5, 0.5))
#' times  <- seq(0, 1, by=0.05)
#' da     <- samHTmcPar3DTr(grid=grid3d, KK=KKmat, Ss=1e-4,
#'                          qHT=list(Qinf), oHT=list(Oinf), times=times, ncore=2)
samHTmcPar3DTr <- function(grid,
                           KK    = 0.1,
                           Ss    = 1e-4,
                           qHT   = list(data.frame(Qp=10, x=7.5, y=7.5, z=2.5)),
                           oHT   = list(data.frame(data=NA, x=5, y=5, z=2.5, time=50)),
                           times = seq(0, 100, by=10),
                           lrw   = 20000000,
                           ncore = 4) {

  library('foreach')
  library('doParallel')
  registerDoParallel(cores = ncore)

  nsim <- ncol(KK)
  x <- foreach(i = 1:nsim, .combine = rbind) %dopar% {
    HydroTomo::samHT3Dtr(grid = grid, KK = KK[, i], Ss = Ss,
                         qHT = qHT, oHT = oHT, times = times, lrw = lrw)
  }

  stopImplicitCluster()
  return(x)
}


#' 3D Transient Ensemble-based Bayesian inverse for hydraulic tomography (parallel)
#'
#' Transient extension of \code{Finverse3D}.  Replaces the steady-state forward
#' solver with \code{Ftransient3dsim} and uses \code{samData3DTr} for time-aware
#' observation sampling.  The Bayesian updating equations (ensemble smoother)
#' are identical to the steady-state version.
#'
#' @param domain  9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2, z1, z2)}.
#'   **Keep small for testing** (e.g. \code{c(15,15,5,0,15,0,15,0,5)}).
#' @param grid    grid from \code{GenGrid3D()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT     list of pumping test data frames with columns
#'   \code{Qp, x, y, z_top, z_bottom} (well-screen pumping) or
#'   \code{Qp, x, y, z} (point pumping).
#' @param nsim    ensemble size (default 50; use ~20 for quick tests).
#' @param itermax maximum iterations (set to 1 for initial testing).
#' @param varmeanTmax  convergence threshold on variance of mean ln(K).
#' @param rmsemin minimum RMSE stopping criterion.
#' @param mul     stabiliser (default 1.0).
#' @param decay   stabiliser decay per iteration (default 1.05).
#' @param oHT     list of observation data frames with columns
#'   \code{data, x, y, z, time}.
#' @param Ss      specific storage Ss \code{[1/L]} — scalar or length-\code{n}
#'   vector.  Default \code{1e-4}.
#' @param times   numeric vector of output times for the ODE solver.  Must
#'   cover all observation times in \code{oHT}.
#' @param lrw     real work array length for \code{ode.3D}
#'   (default 20000000; increase for larger grids).
#' @param ncore   number of parallel cores (default 4).
#' @param geo     prior geostatistical parameters — same list structure as
#'   \code{random3d()}: \code{list(me, var, geomod, range, nugget, anis)}.
#'   Note \code{anis} is a 5-element vector for 3D.
#' @param ifcor    logical; if \code{TRUE}, returns the ensemble results after
#'   the first iteration without performing the update. Useful for checking
#'   forward simulations before running the full inversion.
#' @return list of per-iteration results, each a list with
#'   \code{meanT} (mean ln K vector), \code{varT} (variance of ln K),
#'   \code{meanobsh}, \code{varobsh}.  When \code{ifcor = TRUE}, returns
#'   list with \code{obsh} (simulated heads), \code{Tnew} (ensemble), and
#'   \code{meanobsh}.
#' @export
#' @examples
#' # --- small synthetic test (15x15x5, 1 iteration) ---
#' domain3d <- c(15, 15, 5, 0, 15, 0, 15, 0, 5)
#' grid3d   <- GenGrid3D(domain3d)
#' set.seed(42)
#' trueK3d  <- random3d(nsim=1, grid=grid3d)
#' Qinf3d   <- data.frame(Qp=10, x=7.5, y=7.5, z=2.5)
#' qHT3d    <- list(test1 = Qinf3d)
#' times    <- seq(0, 1, by=0.05)
#' res3d    <- Ftransient3dsim(grid=grid3d, KK=trueK3d$Kp, Ss=1e-4,
#'                             Qinf=Qinf3d, times=times)
#' loc      <- expand.grid(x=c(3,6,9,12), y=c(3,6,9,12))
#' Oinf3d   <- data.frame(data=NA, x=loc$x, y=loc$y, z=2.5, time=0.5)
#' Oinf3d   <- samData3DTr(Oinf=Oinf3d, grid=grid3d, result_tr=res3d)
#' oHT3d    <- list(test1 = Oinf3d)
#' result3d <- Finverse3DTr(grid=grid3d, qHT=qHT3d, oHT=oHT3d,
#'                          Ss=1e-4, times=times, nsim=20, itermax=1, ncore=2)
Finverse3DTr <- function(
    domain      = c(15, 15, 5, 0, 15, 0, 15, 0, 5),
    grid        = NULL,
    qHT         = list(data.frame(Qp=10, x=7.5, y=7.5, z=2.5)),
    nsim        = 50,
    itermax     = 5,
    varmeanTmax = 5,
    rmsemin     = 0,
    mul         = 1.0,
    decay       = 1.05,
    oHT         = list(data.frame(data=-1, x=5, y=5, z=2.5, time=50)),
    Ss          = 1e-4,
    times       = seq(0, 100, by=10),
    lrw         = 20000000,
    ncore       = 4,
    geo         = list(me=0, var=1, geomod="Exp",
                       range=10, nugget=0,
                       anis=c(0, 0, 0, 1, 1)),
    ifcor       = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid3D(domain)

  nHT <- length(qHT)

  # ---- extract observation data ----------------------------------------------
  trueobshHT <- list()
  for (i in seq_len(nHT)) {
    oinf            <- oHT[[i]]
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs     <- length(trueobsh)

  # ---- initial ensemble ------------------------------------------------------
  yy   <- random3d(nsim = nsim, grid = grid, geo = geo)
  Knew <- as.matrix(yy[, -c(1, 2, 3)])   # n × nsim matrix of K values

  # ---- iteration loop --------------------------------------------------------
  niter    <- 1
  varmeanT <- 0
  rmse     <- 1e10
  msgdf    <- data.frame(niter = niter, varmeanT = varmeanT, rmse = rmse)
  iterdf   <- list()

  while (niter <= itermax & varmeanT < varmeanTmax & rmse > rmsemin) {

    varT     <- apply(log(Knew), 1, var)
    meanT    <- apply(log(Knew), 1, mean)
    varmeanT <- var(meanT)

    # parallel transient forward runs for all ensemble members
    obsh <- samHTmcPar3DTr(grid  = grid,
                           KK    = Knew,
                           Ss    = Ss,
                           qHT   = qHT,
                           oHT   = oHT,
                           times = times,
                           lrw   = lrw,
                           ncore = ncore)

    # ---- statistics ----------------------------------------------------------
    varobsh  <- apply(obsh, 2, var)
    meanobsh <- apply(obsh, 2, mean)
    weigs    <- 1 / varobsh / sum(1 / varobsh)
    rmse     <- mean((trueobsh - meanobsh)^2 * weigs)^0.5
    l2       <- mean((trueobsh - meanobsh)^2)^0.5
    l1       <- mean(abs(trueobsh - meanobsh))

    # ---- Bayesian update (same as steady-state) ------------------------------
    covh  <- cov(obsh)
    covhk <- cov(obsh, t(log(Knew)))
    if (ifcor) {
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }
    covh1 <- covh
    diag(covh1) <- rep((1 + mul) * max(diag(covh)), nobs)
    a <- solve(covh1, covhk)   # nobs × n

    for (i in seq_len(nsim)) {
      Knew[, i] <- Knew[, i] * exp(t(a) %*% (trueobsh - obsh[i, ]))
    }

    msg <- paste('niter =', niter,
                 'varmeanT =', round(varmeanT, 4),
                 'rmse =', round(rmse, 4),
                 'l2 =', round(l2, 4),
                 'l1 =', round(l1, 4))
    print(msg)

    iterdf[[niter]] <- list(meanT    = as.vector(meanT),
                            varT     = as.vector(varT),
                            meanobsh = meanobsh,
                            varobsh  = varobsh)

    if (niter == 1) {
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    mul   <- mul / decay
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}


# ==============================================================================




# ==============================================================================
# 3D Transient HT — Well-Screen (vertical interval) Inverse Functions
# ==============================================================================

#' Run forward simulations for one K field — 3D transient well-screen version
#'
#' Like \code{\link{samHT3Dtr}} but for wells with vertical screens.
#' Calls \code{\link{Ftransient3dsim}} and samples interval-averaged heads
#' via \code{\link{samData3DTrScreen}}.
#'
#' @param grid   grid from \code{GenGrid3D()}.
#' @param KK     length-\code{n} vector of hydraulic conductivity K \code{[L/T]}.
#' @param Ss     specific storage Ss \code{[1/L]} — scalar or length-\code{n}
#'   vector.
#' @param qHT    list of pumping test data frames, each with columns
#'   \code{Qp, x, y, z_top, z_bottom} (well-screen pumping) or
#'   \code{Qp, x, y, z} (point pumping).
#' @param oHT    list of observation data frames, each with columns
#'   \code{data, x, y, z_top, z_bottom, time}.
#' @param times  numeric vector of output times for the ODE solver.
#' @param lrw    real work array length; default 20000000.
#' @param simplify if \code{TRUE} return a plain vector; \code{FALSE} returns
#'   updated \code{oHT}.
#' @return vector of interval-averaged drawdown values (when \code{simplify = TRUE}).
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15, 15, 5, 0, 15, 0, 15, 0, 5))
#' Qinf   <- data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)
#' Oinf   <- data.frame(data=NA, x=c(5,10), y=c(7.5,7.5),
#'                      z_top=c(1,1), z_bottom=c(4,4), time=c(0.5, 0.5))
#' times  <- seq(0, 1, by=0.05)
#' da     <- samHT3DtrScreen(grid=grid3d, KK=0.1, Ss=1e-4,
#'                           qHT=list(Qinf), oHT=list(Oinf), times=times)
samHT3DtrScreen <- function(grid,
                            KK       = 0.1,
                            Ss       = 1e-4,
                            qHT      = list(data.frame(Qp=10, x=7.5, y=7.5,
                                                      z_top=2, z_bottom=3)),
                            oHT      = list(data.frame(data=NA, x=5, y=5,
                                                      z_top=1, z_bottom=4, time=50)),
                            times    = seq(0, 100, by=10),
                            lrw      = 20000000,
                            simplify = TRUE) {

  nHT <- length(qHT)
  for (i in seq_len(nHT)) {
    qinf <- qHT[[i]]
    oinf <- oHT[[i]]
    res  <- Ftransient3dsim(grid = grid, KK = KK, Ss = Ss,
                            Qinf = qinf, times = times, lrw = lrw)
    oinf <- samData3DTrScreen(Oinf = oinf, grid = grid, result_tr = res)
    oHT[[i]] <- oinf[, c('data', 'x', 'y', 'z_top', 'z_bottom', 'time')]
  }
  if (simplify) {
    oHTdf <- dplyr::bind_rows(oHT, .id = 'id')
    return(oHTdf$data)
  } else {
    return(oHT)
  }
}


#' Parallel MC forward runs for 3D transient HT — well-screen version (ensemble)
#'
#' Like \code{\link{samHTmcPar3DTr}} but uses \code{\link{samHT3DtrScreen}}
#' for interval-averaged well-screen sampling.
#'
#' @param grid   grid from \code{GenGrid3D()}.
#' @param KK     \code{n × nsim} matrix of K realisations.
#' @param Ss     specific storage Ss \code{[1/L]} — scalar or length-\code{n}
#'   vector.
#' @param qHT    list of pumping test data frames (columns
#'   \code{Qp, x, y, z_top, z_bottom} for well-screen pumping, or
#'   \code{Qp, x, y, z} for point pumping).
#' @param oHT    list of observation data frames (columns
#'   \code{data, x, y, z_top, z_bottom, time}).
#' @param times  numeric vector of output times for the ODE solver.
#' @param lrw    real work array length; default 20000000.
#' @param ncore  number of parallel cores.
#' @return \code{nsim × nobs} matrix of interval-averaged heads.
#' @export
#' @examples
#' grid3d <- GenGrid3D(c(15,15,5,0,15,0,15,0,5))
#' KKmat  <- random3d(nsim=5, grid=grid3d)
#' KKmat  <- as.matrix(KKmat[,-c(1,2,3)])
#' Qinf   <- data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)
#' Oinf   <- data.frame(data=NA, x=c(5,10), y=c(7.5,7.5),
#'                      z_top=c(1,1), z_bottom=c(4,4), time=c(0.5, 0.5))
#' times  <- seq(0, 1, by=0.05)
#' da     <- samHTmcPar3DTrScreen(grid=grid3d, KK=KKmat, Ss=1e-4,
#'                                qHT=list(Qinf), oHT=list(Oinf),
#'                                times=times, ncore=2)
samHTmcPar3DTrScreen <- function(grid,
                                 KK    = 0.1,
                                 Ss    = 1e-4,
                                 qHT   = list(data.frame(Qp=10, x=7.5, y=7.5,
                                                         z_top=2, z_bottom=3)),
                                 oHT   = list(data.frame(data=NA, x=5, y=5,
                                                        z_top=1, z_bottom=4, time=50)),
                                 times = seq(0, 100, by=10),
                                 lrw   = 20000000,
                                 ncore = 4) {

  library('foreach')
  library('doParallel')
  registerDoParallel(cores = ncore)

  nsim <- ncol(KK)
  x <- foreach(i = 1:nsim, .combine = rbind) %dopar% {
    HydroTomo::samHT3DtrScreen(grid = grid, KK = KK[, i], Ss = Ss,
                               qHT = qHT, oHT = oHT, times = times, lrw = lrw)
  }

  stopImplicitCluster()
  return(x)
}


#' 3D Transient HT ensemble inverse — well-screen (vertical interval) version
#'
#' Like \code{\link{Finverse3DTr}} but for wells with vertical screens
#' (filter sections).  Observation wells are defined by \code{(x, y)} and a
#' screen interval \code{[z_top, z_bottom]}.  The forward solver
#' \code{\link{Ftransient3dsim}} computes full 3D transient heads, and
#' \code{\link{samData3DTrScreen}} averages heads over all grid cells within
#' the screen interval.
#'
#' @param domain  9-element vector \code{c(nx, ny, nz, x1, x2, y1, y2, z1, z2)}.
#' @param grid    grid from \code{GenGrid3D()}; generated from \code{domain}
#'   if \code{NULL}.
#' @param qHT     list of pumping test data frames with columns
#'   \code{Qp, x, y, z_top, z_bottom} (well-screen pumping) or
#'   \code{Qp, x, y, z} (point pumping).
#' @param nsim    ensemble size (default 50).
#' @param itermax maximum iterations.
#' @param varmeanTmax  convergence threshold on variance of mean ln(K).
#' @param rmsemin minimum RMSE stopping criterion.
#' @param mul     stabiliser (default 1.0).
#' @param decay   stabiliser decay per iteration (default 1.05).
#' @param oHT     list of observation data frames with columns
#'   \code{data, x, y, z_top, z_bottom, time}.
#' @param Ss      specific storage Ss \code{[1/L]} — scalar or length-\code{n}
#'   vector.
#' @param times   numeric vector of output times for the ODE solver.
#' @param lrw     real work array length for \code{ode.3D} (default 20000000).
#' @param ncore   number of parallel cores (default 4).
#' @param geo     prior geostatistical parameters (list).
#' @param ifcor    logical; if \code{TRUE}, returns the ensemble results after
#'   the first iteration without performing the update.
#' @return list of per-iteration results.  When \code{ifcor = TRUE}, returns
#'   list with \code{obsh}, \code{Tnew}, and \code{meanobsh}.
#' @export
#' @examples
#' # --- small synthetic test with well screens ---
#' domain3d <- c(15, 15, 5, 0, 15, 0, 15, 0, 5)
#' grid3d   <- GenGrid3D(domain3d)
#' set.seed(42)
#' trueK3d  <- random3d(nsim=1, grid=grid3d)
#' Qinf3d   <- data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)
#' qHT3d    <- list(test1 = Qinf3d)
#' times    <- seq(0, 1, by=0.05)
#' res3d    <- Ftransient3dsim(grid=grid3d, KK=trueK3d$Kp, Ss=1e-4,
#'                             Qinf=Qinf3d, times=times)
#' loc      <- expand.grid(x=c(3,6,9,12), y=c(3,6,9,12))
#' Oinf3d   <- data.frame(data=NA, x=loc$x, y=loc$y,
#'                        z_top=1, z_bottom=4, time=0.5)
#' Oinf3d   <- samData3DTrScreen(Oinf=Oinf3d, grid=grid3d, result_tr=res3d)
#' oHT3d    <- list(test1 = Oinf3d)
#' result3d <- Finverse3DTrScreen(grid=grid3d, qHT=qHT3d, oHT=oHT3d,
#'                                Ss=1e-4, times=times, nsim=20, itermax=1, ncore=2)
Finverse3DTrScreen <- function(
    domain      = c(15, 15, 5, 0, 15, 0, 15, 0, 5),
    grid        = NULL,
    qHT         = list(data.frame(Qp=10, x=7.5, y=7.5, z_top=2, z_bottom=3)),
    nsim        = 50,
    itermax     = 5,
    varmeanTmax = 5,
    rmsemin     = 0,
    mul         = 1.0,
    decay       = 1.05,
    oHT         = list(data.frame(data=-1, x=5, y=5,
                                 z_top=1, z_bottom=4, time=50)),
    Ss          = 1e-4,
    times       = seq(0, 100, by=10),
    lrw         = 20000000,
    ncore       = 4,
    geo         = list(me=0, var=1, geomod="Exp",
                       range=10, nugget=0,
                       anis=c(0, 0, 0, 1, 1)),
    ifcor       = FALSE) {

  set.seed(200)
  startTime <- Sys.time()

  if (is.null(grid)) grid <- GenGrid3D(domain)

  nHT <- length(qHT)

  # ---- extract observation data ----------------------------------------------
  trueobshHT <- list()
  for (i in seq_len(nHT)) {
    oinf            <- oHT[[i]]
    trueobshHT[[i]] <- oinf$data
  }
  trueobsh <- unlist(trueobshHT)
  nobs     <- length(trueobsh)

  # ---- initial ensemble ------------------------------------------------------
  yy   <- random3d(nsim = nsim, grid = grid, geo = geo)
  Knew <- as.matrix(yy[, -c(1, 2, 3)])

  # ---- iteration loop --------------------------------------------------------
  niter    <- 1
  varmeanT <- 0
  rmse     <- 1e10
  msgdf    <- data.frame(niter = niter, varmeanT = varmeanT, rmse = rmse)
  iterdf   <- list()

  while (niter <= itermax & varmeanT < varmeanTmax & rmse > rmsemin) {

    varT     <- apply(log(Knew), 1, var)
    meanT    <- apply(log(Knew), 1, mean)
    varmeanT <- var(meanT)

    # parallel transient forward runs — screen-averaged
    obsh <- samHTmcPar3DTrScreen(grid  = grid,
                                 KK    = Knew,
                                 Ss    = Ss,
                                 qHT   = qHT,
                                 oHT   = oHT,
                                 times = times,
                                 lrw   = lrw,
                                 ncore = ncore)

    # ---- statistics ----------------------------------------------------------
    varobsh  <- apply(obsh, 2, var)
    meanobsh <- apply(obsh, 2, mean)
    weigs    <- 1 / varobsh / sum(1 / varobsh)
    rmse     <- mean((trueobsh - meanobsh)^2 * weigs)^0.5
    l2       <- mean((trueobsh - meanobsh)^2)^0.5
    l1       <- mean(abs(trueobsh - meanobsh))

    # ---- Bayesian update -----------------------------------------------------
    covh  <- cov(obsh)
    covhk <- cov(obsh, t(log(Knew)))
    if (ifcor) {
      print("ifcor = TRUE: returning cross-covariance covhk.")
      return(covhk)
    }
    covh1 <- covh
    diag(covh1) <- rep((1 + mul) * max(diag(covh)), nobs)
    a <- solve(covh1, covhk)

    for (i in seq_len(nsim)) {
      Knew[, i] <- Knew[, i] * exp(t(a) %*% (trueobsh - obsh[i, ]))
    }

    msg <- paste('niter =', niter,
                 'varmeanT =', round(varmeanT, 4),
                 'rmse =', round(rmse, 4),
                 'l2 =', round(l2, 4),
                 'l1 =', round(l1, 4))
    print(msg)

    iterdf[[niter]] <- list(meanT    = as.vector(meanT),
                            varT     = as.vector(varT),
                            meanobsh = meanobsh,
                            varobsh  = varobsh)

    if (niter == 1) {
      print("--- time for one iteration ---")
      print(difftime(Sys.time(), startTime))
    }

    niter <- niter + 1
    mul   <- mul / decay
  }

  print("--- total time ---")
  print(difftime(Sys.time(), startTime))
  return(iterdf)
}