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
#'   \code{Qp, x, y, z}.
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
#'   \code{Qp, x, y, z}.
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
