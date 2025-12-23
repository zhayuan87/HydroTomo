#' make a Monte Carlo inverse simulation.
#' HT simulations.
#' forward information. domain, bc, pump loc,
#' observation location and value.
#' geostatistical prior information.
#' hyper parameters --- nsim, itermax, etc.
#' @param domain grid c(nx,ny,x1,x2,y1,y2)
#' @param qHT list of pumping test information.
#' @param nsim number of ensemble, default 50.
#' @param itermax maximum number of iterations.
#' @param varmeanTmax maximum variance of mean T.
#' @param rmsemin minimum rmse.
#' @param mul stablizer.
#' @param loc_obsHT list of observation location.
#' @param trueobshHT list of observation value.
#' @return list of iteration data.
#' @export
#' @examples
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT <- trueK[,-c(1,2)]
#' Qinf1=list(Qp=10,xp=20.5,yp=20.5)
#' trueh <- Fsteady2dsim(TT=TT,Qinf=Qinf1)
#' trueh <- trueh$solution
#' n <- 1600
#' nobs <- 25
#' locx = c(15,18,22,25,30)
#' locy = c(15,18,22,25,30)
#' loc= expand.grid(x=locx,y=locy)
#' loc_obs <- (loc$y-1)*40 + loc$x
#' trueobsh <- trueh[loc_obs]
#' loc_obsHT <- list(loc_obs)
#' trueobshHT <- list(trueobsh)
#' result <- Finverse(loc_obs=loc_obsHT,trueobsh=trueobshHT)
#' Qinf2=list(Qp=10,xp=10.5,yp=10.5)
#' trueh2 <- Fsteady2dsim(TT=TT,Qinf=Qinf2)
#' trueh2 <- trueh2$solution
#' trueobsh2 <- trueh2[loc_obs]
#' loc_obsHT <- list(loc_obs,loc_obs)
#' trueobshHT <- list(trueobsh,trueobsh2)
#' result <- Finverse(loc_obs=loc_obsHT,trueobsh=trueobshHT,qHT = list(Qinf1,Qinf2))

Finverse <- function(
    domain=c(40,40,0,40,0,40),
    qHT=list(list(Qp=10,xp=20.5,yp=20.5)), # should be double list.
    nsim=50,
    itermax=5,
    varmeanTmax =5,
    rmsemin = 0,
    mul=1.0,
    loc_obsHT, # should be list since multiple pumping test.
    trueobshHT) # should be list since multiple pumping test.
{
  set.seed(200)
  ### record the time so to see how long it takes.
  startTime <- Sys.time()

  #nsim = 50
  nHT <- length(qHT) # number of pumping test.
  yy <- random2d(nsim=nsim,domain = domain)
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
  trueobsh <- unlist(trueobshHT) # the data format in HT is a list.
  nobs <- length(trueobsh)
  # start of the itertion loop.
  while(niter<=itermax & varmeanT<varmeanTmax & rmse>rmsemin){
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
    hHT <- lapply(1:nHT, function(j) {
      Qinf <- qHT[[j]]
      # 使用 lapply 生成所有模拟结果
      h_list <- lapply(1:nsim, function(i)
        Fsteady2dsim(TT = Tnew[, i], Qinf = Qinf)$solution
      )
      # 一次性合并所有结果
      h_matrix <- do.call(rbind, h_list)
      # 提取观测点数据
      h_matrix[, loc_obsHT[[j]]]
    })


    # obsh: nsim*nobs
    if(nHT>1) obsh <- do.call("cbind",hHT) else obsh <- hHT[[1]]

    # get the head variance for each observation.
    varobsh <- apply(obsh,2,var)
    meanobsh <- apply(obsh,2,mean)
    # get the misfit.
    weigs <- 1/varobsh/sum(1/varobsh)
    rmse <- mean((trueobsh - meanobsh)^2*weigs)^0.5
    l2 <- mean((trueobsh - meanobsh)^2)^0.5
    l1 <- mean(abs(trueobsh - meanobsh))

    # get the covariance of h and h-T
    covh <- cov(obsh)
    data = list(obsh,t(log(Tnew)))
    # notice it is lnK rather than K.
    cc <- function(obsh,TT){
      covhk = cov(obsh,TT)
    }
    covhk = suppressWarnings(multiApply::Apply(data,target_dims = list(1,1),cc))

    # add stablizer term.
    covh1 <- covh
    #diag(covh1) <-  (1+mul)*diag(covh)
    diag(covh1) <-  rep((1+mul)*max(diag(covh)),nobs)
    #### now we need to invert the covariance function and get the covh-1* %*% covhf
    a = solve(covh1,covhk[[1]])
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
    print(msg)
    ### we need to store the iteration data.
    iterdf[[niter]] <- list(meanT = as.vector(meanT), varT = as.vector(varT), meanobsh = meanobsh, varobsh = varobsh)
    niter <- niter + 1
    mul <- mul/1.05

    ### iteration over time.
  }
  endTime <- Sys.time()
  # get the time in seconds or mins.
  print(paste("The computational time takes ", round(endTime-startTime,digits=2),"seconds."))
  return(iterdf)
}
