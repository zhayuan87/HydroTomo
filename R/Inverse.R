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
#' @return A list of per-iteration results, each element a list with
#'   \code{meanT} (mean ln(T) vector), \code{varT} (variance of ln(T)),
#'   \code{meanobsh} (mean simulated head), \code{varobsh} (variance of
#'   simulated head).
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
    geo=list(me=0,var=1,geomod="Exp",anis=c(90,1),range=30,nugget=0)) # should be list since multiple pumping test.
{
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
    geo=list(me=0,var=1,geomod="Exp",anis=c(90,1),range=30,nugget=0)) # should be list since multiple pumping test.
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
    geo=list(me=0,var=1,geomod="Exp",anis=c(90,1),range=30,nugget=0)) # should be list since multiple pumping test.
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
#' @return list of per-iteration results, each a list with
#'   \code{meanT} (mean ln K vector), \code{varT} (variance of ln K),
#'   \code{meanobsh}, \code{varobsh}.
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
                       anis=c(0, 0, 0, 1, 1))) {

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
#' @return list of per-iteration results.
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
                       anis=c(0, 0, 0, 1, 1))) {

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
  Qp_adj <- dx * dy

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
#' # ---- 1. Setup grid and transmissivity field ----
#' domain <- c(20, 20, 0, 20, 0, 20)
#' grid   <- GenGrid(domain)
#' set.seed(123)
#' TT <- random2d(grid = grid, nsim = 1)$Tp
#'
#' # ---- 2. Define pumping tests and observations ----
#' qHT <- list(
#'   data.frame(Qp = 10, x = 10.5, y = 10.5),
#'   data.frame(Qp = 10, x = 15.5, y = 15.5)
#' )
#' oHT <- list(
#'   data.frame(data = -1, x = 5.5, y = 5.5),
#'   data.frame(data = -1, x = 12.5, y = 12.5)
#' )
#'
#' # ---- 3. Compute Jacobian ----
#' J <- jacobian2D(grid = grid, TT = TT, qHT = qHT, oHT = oHT)
#'
#' # ---- 4. Inspect results ----
#' # Jacobian dimensions: nobs × nelem
#' dim(J)
#' # Sensitivity map for the first observation (row 1)
#' image(matrix(J[1, ], grid$nx, grid$ny),
#'       main = "Sensitivity: obs 1")
#' @export
jacobian2D <- function(grid,
                       TT,
                       qHT = list(data.frame(Qp = 10, x = 20.5, y = 20.5)),
                       oHT = list(data.frame(data = -1, x = 11, y = 11)),
                       lrw = 160000) {

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
#' @return A list of per-iteration results, each a list with
#'   \code{m} (log-T estimate), \code{h_sim} (simulated heads),
#'   \code{rmse}, \code{l2}, \code{l1}, and \code{J} (Jacobian, first iteration
#'   only).  When \code{ifcor = TRUE}, returns list(J = J, h_sim = h_sim).
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
    ifcor       = FALSE) {

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
  # Use the same variogram model as random2d to compute covariance
  require('gstat')
  xy <- grid$grid

  # Compute prior covariance matrix analytically from variogram model
  # C_kk(i,j) = sill - variogram(|x_i - x_j|)
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
      print("ifcor = TRUE: returning Jacobian and simulated heads without update.")
      return(list(J = J, h_sim = h_sim_vec))
    }

    # --- Compute update (dual / observation-space formulation) ---
    # dm = C_kk * J^T * (J * C_kk * J^T + R + lambda_lm * I)^{-1} * residual
    Jt <- t(J)                                     # nelem x nobs
    JCJt <- J %*% C_kk %*% Jt                      # nobs x nobs
    S <- JCJt + (sigma2obs + lambda_lm) * diag(nobs)

    beta <- solve(S, residual)                     # nobs
    dm <- as.vector(C_kk %*% Jt %*% beta)          # nelem

    # --- Update parameters ---
    m <- m + dm
    T_current <- exp(m)

    # --- Store iteration results ---
    iterdf[[niter]] <- list(
      m       = m,
      h_sim   = h_sim_vec,
      rmse    = rmse,
      l2      = l2,
      l1      = l1
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
                       anis=c(90,1), range=30, nugget=0)) {

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
  Qp_adj <- dx * dy

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
#' # ---- 1. Setup grid and transmissivity field ----
#' domain <- c(20, 20, 0, 20, 0, 20)
#' grid   <- GenGrid(domain)
#' set.seed(123)
#' TT <- random2d(grid = grid, nsim = 1)$Tp
#'
#' # ---- 2. Define pumping tests and observations ----
#' qHT <- list(
#'   data.frame(Qp = 10, x = 10.5, y = 10.5),
#'   data.frame(Qp = 10, x = 15.5, y = 15.5)
#' )
#' oHT <- list(
#'   data.frame(data = -1, x = 5.5, y = 5.5),
#'   data.frame(data = -1, x = 12.5, y = 12.5)
#' )
#'
#' # ---- 3. Compute Jacobian ----
#' J <- jacobian2D(grid = grid, TT = TT, qHT = qHT, oHT = oHT)
#'
#' # ---- 4. Inspect results ----
#' # Jacobian dimensions: nobs × nelem
#' dim(J)
#' # Sensitivity map for the first observation (row 1)
#' image(matrix(J[1, ], grid$nx, grid$ny),
#'       main = "Sensitivity: obs 1")
#' @export
jacobian2D <- function(grid,
                       TT,
                       qHT = list(data.frame(Qp = 10, x = 20.5, y = 20.5)),
                       oHT = list(data.frame(data = -1, x = 11, y = 11)),
                       lrw = 160000) {

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
#' @return A list of per-iteration results, each a list with
#'   \code{m} (log-T estimate), \code{h_sim} (simulated heads),
#'   \code{rmse}, \code{l2}, \code{l1}, and \code{J} (Jacobian, first iteration
#'   only).  When \code{ifcor = TRUE}, returns list(J = J, h_sim = h_sim).
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
    ifcor       = FALSE) {

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
  # Use the same variogram model as random2d to compute covariance
  require('gstat')
  xy <- grid$grid

  # Compute prior covariance matrix analytically from variogram model
  # C_kk(i,j) = sill - variogram(|x_i - x_j|)
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
      print("ifcor = TRUE: returning Jacobian and simulated heads without update.")
      return(list(J = J, h_sim = h_sim_vec))
    }

    # --- Compute update (dual / observation-space formulation) ---
    # dm = C_kk * J^T * (J * C_kk * J^T + R + lambda_lm * I)^{-1} * residual
    Jt <- t(J)                                     # nelem x nobs
    JCJt <- J %*% C_kk %*% Jt                      # nobs x nobs
    S <- JCJt + (sigma2obs + lambda_lm) * diag(nobs)

    beta <- solve(S, residual)                     # nobs
    dm <- as.vector(C_kk %*% Jt %*% beta)          # nelem

    # --- Update parameters ---
    m <- m + dm
    T_current <- exp(m)

    # --- Store iteration results ---
    iterdf[[niter]] <- list(
      m       = m,
      h_sim   = h_sim_vec,
      rmse    = rmse,
      l2      = l2,
      l1      = l1
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
#' @return list of per-iteration results, each a list with
#'   \code{meanT} (mean ln K vector), \code{varT} (variance of ln K),
#'   \code{meanobsh}, \code{varobsh}.
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
                       anis=c(0, 0, 0, 1, 1))) {

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
  Qp_adj <- dx * dy

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
#' # ---- 1. Setup grid and transmissivity field ----
#' domain <- c(20, 20, 0, 20, 0, 20)
#' grid   <- GenGrid(domain)
#' set.seed(123)
#' TT <- random2d(grid = grid, nsim = 1)$Tp
#'
#' # ---- 2. Define pumping tests and observations ----
#' qHT <- list(
#'   data.frame(Qp = 10, x = 10.5, y = 10.5),
#'   data.frame(Qp = 10, x = 15.5, y = 15.5)
#' )
#' oHT <- list(
#'   data.frame(data = -1, x = 5.5, y = 5.5),
#'   data.frame(data = -1, x = 12.5, y = 12.5)
#' )
#'
#' # ---- 3. Compute Jacobian ----
#' J <- jacobian2D(grid = grid, TT = TT, qHT = qHT, oHT = oHT)
#'
#' # ---- 4. Inspect results ----
#' # Jacobian dimensions: nobs × nelem
#' dim(J)
#' # Sensitivity map for the first observation (row 1)
#' image(matrix(J[1, ], grid$nx, grid$ny),
#'       main = "Sensitivity: obs 1")
#' @export
jacobian2D <- function(grid,
                       TT,
                       qHT = list(data.frame(Qp = 10, x = 20.5, y = 20.5)),
                       oHT = list(data.frame(data = -1, x = 11, y = 11)),
                       lrw = 160000) {

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
#' @return A list of per-iteration results, each a list with
#'   \code{m} (log-T estimate), \code{h_sim} (simulated heads),
#'   \code{rmse}, \code{l2}, \code{l1}, and \code{J} (Jacobian, first iteration
#'   only).  When \code{ifcor = TRUE}, returns list(J = J, h_sim = h_sim).
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
    ifcor       = FALSE) {

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
  # Use the same variogram model as random2d to compute covariance
  require('gstat')
  xy <- grid$grid

  # Compute prior covariance matrix analytically from variogram model
  # C_kk(i,j) = sill - variogram(|x_i - x_j|)
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
      print("ifcor = TRUE: returning Jacobian and simulated heads without update.")
      return(list(J = J, h_sim = h_sim_vec))
    }

    # --- Compute update (dual / observation-space formulation) ---
    # dm = C_kk * J^T * (J * C_kk * J^T + R + lambda_lm * I)^{-1} * residual
    Jt <- t(J)                                     # nelem x nobs
    JCJt <- J %*% C_kk %*% Jt                      # nobs x nobs
    S <- JCJt + (sigma2obs + lambda_lm) * diag(nobs)

    beta <- solve(S, residual)                     # nobs
    dm <- as.vector(C_kk %*% Jt %*% beta)          # nelem

    # --- Update parameters ---
    m <- m + dm
    T_current <- exp(m)

    # --- Store iteration results ---
    iterdf[[niter]] <- list(
      m       = m,
      h_sim   = h_sim_vec,
      rmse    = rmse,
      l2      = l2,
      l1      = l1
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
#' @return list of per-iteration results.
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
                       anis=c(0, 0, 0, 1, 1))) {

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
  Qp_adj <- dx * dy

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
#' # ---- 1. Setup grid and transmissivity field ----
#' domain <- c(20, 20, 0, 20, 0, 20)
#' grid   <- GenGrid(domain)
#' set.seed(123)
#' TT <- random2d(grid = grid, nsim = 1)$Tp
#'
#' # ---- 2. Define pumping tests and observations ----
#' qHT <- list(
#'   data.frame(Qp = 10, x = 10.5, y = 10.5),
#'   data.frame(Qp = 10, x = 15.5, y = 15.5)
#' )
#' oHT <- list(
#'   data.frame(data = -1, x = 5.5, y = 5.5),
#'   data.frame(data = -1, x = 12.5, y = 12.5)
#' )
#'
#' # ---- 3. Compute Jacobian ----
#' J <- jacobian2D(grid = grid, TT = TT, qHT = qHT, oHT = oHT)
#'
#' # ---- 4. Inspect results ----
#' # Jacobian dimensions: nobs × nelem
#' dim(J)
#' # Sensitivity map for the first observation (row 1)
#' image(matrix(J[1, ], grid$nx, grid$ny),
#'       main = "Sensitivity: obs 1")
#' @export
jacobian2D <- function(grid,
                       TT,
                       qHT = list(data.frame(Qp = 10, x = 20.5, y = 20.5)),
                       oHT = list(data.frame(data = -1, x = 11, y = 11)),
                       lrw = 160000) {

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
#' @return A list of per-iteration results, each a list with
#'   \code{m} (log-T estimate), \code{h_sim} (simulated heads),
#'   \code{rmse}, \code{l2}, \code{l1}, and \code{J} (Jacobian, first iteration
#'   only).  When \code{ifcor = TRUE}, returns list(J = J, h_sim = h_sim).
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
    ifcor       = FALSE) {

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
  # Use the same variogram model as random2d to compute covariance
  require('gstat')
  xy <- grid$grid

  # Compute prior covariance matrix analytically from variogram model
  # C_kk(i,j) = sill - variogram(|x_i - x_j|)
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
      print("ifcor = TRUE: returning Jacobian and simulated heads without update.")
      return(list(J = J, h_sim = h_sim_vec))
    }

    # --- Compute update (dual / observation-space formulation) ---
    # dm = C_kk * J^T * (J * C_kk * J^T + R + lambda_lm * I)^{-1} * residual
    Jt <- t(J)                                     # nelem x nobs
    JCJt <- J %*% C_kk %*% Jt                      # nobs x nobs
    S <- JCJt + (sigma2obs + lambda_lm) * diag(nobs)

    beta <- solve(S, residual)                     # nobs
    dm <- as.vector(C_kk %*% Jt %*% beta)          # nelem

    # --- Update parameters ---
    m <- m + dm
    T_current <- exp(m)

    # --- Store iteration results ---
    iterdf[[niter]] <- list(
      m       = m,
      h_sim   = h_sim_vec,
      rmse    = rmse,
      l2      = l2,
      l1      = l1
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
