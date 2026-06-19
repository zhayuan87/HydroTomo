library(HydroTomo)
domain = c(41,41,0,41,0,41)
grid = GenGrid(domain=domain)
for(iseed in 44:60){
  set.seed(iseed)
  if(!dir.exists(as.character(iseed)))dir.create(as.character(iseed))
  geo = list(me = 3.2, var = 3, geomod = "Sph", anis = c(90, 1), range = 10, nugget = 0.02)
  trueK <- random2d(nsim=1,grid=grid,geo = geo)
  # plot the generated log -T field.
  #Plotparameter2d(trueK,iflog=T)
  x = c(11,21,31) - 0.5
  y = x
  xy = expand.grid(x=x,y=y)

  qHT1 = list()
  for(i in 1:nrow(xy)){
    Qp = 200
    x =xy$x[i]
    y =xy$y[i]
    qHT1[[i]] = data.frame(x=x,y=y,Qp=Qp)
  }

  ### in the obs. pumping at one, and obs. at the rest...

  oHT1 = list()

  for(i in 1:nrow(xy)){
    data = NA
    x =xy$x[-i]
    y =xy$y[-i]
    oHT1[[i]] = data.frame(x=x,y=y,data=data)
  }

  ##### now sample the data.
  for (i in 1:length(oHT1)){
    Qinf <- qHT1[[i]]
    trueh <- Fsteady2dsim(TT=trueK$Tp,Qinf=Qinf,grid=grid,lrw=160000)$solution
    oinf <- oHT1[[i]]
    oinf = samData(grid = grid,Oinf = oinf,h=trueh)
    oHT1[[i]] <- oinf[,c('data','x','y')]
  }

  dfoHT1 = do.call(rbind,oHT1)


  # figure 1.......
  library(ggplot2)
  library(showtext)
  showtext_auto()
  df = trueK
  df$Tp = log(trueK$Tp)
  qHTdf = data.frame(id= paste0("W",1:9),xy)
  p1 = ggplot() +
    geom_tile(data = df, aes(x, y, fill = Tp)) +
    scale_fill_viridis_c(option = "D")+
    geom_point(data = qHTdf, aes(x = x, y = y), color = "white", size = 3, shape=17)+
    geom_text(data = qHTdf, aes(x = x, y = y, label = id), vjust = -1, color = "black")+
    labs(title = "(a) 参考lnT及井布置", x = "X/m", y = "Y/m", fill = "log(T)")+
    theme_minimal()
  p1

  s <- Fsteady2dsim(TT=trueK$Tp,Qinf=qHT1[[1]],grid=grid,lrw=160000)
  s$`降深/m` = -s$solution
  p2 = ggplot() +
    geom_tile(data = s, aes(x, y, fill = `降深/m`)) +
    scale_fill_viridis_c(option = "D",trans="log10")+
    labs(title = "(b) 模拟降深及观测井", x = "X/m", y = "Y/m")+
    theme_minimal()
  p2

  library(ggpubr)

  p = ggarrange(p1, p2, ncol = 2, nrow = 1)

  ggsave(filename = paste0(iseed,"/figure1.pdf"),plot = p,
         width = 10, height = 5)
  #### now do the inverse.

  lire = list()
  ### W1, 3,5,7
  se = c(1,3,5,7)
  # we select qHT1 and oHT1 from blh_synthetic, first test1, then 1,2, then 1,2,3
  for(i in 1:4){
    ii = se[1:i]
    #ii =1:i
    qHT11 = qHT1[ii]
    oHT11 = oHT1[ii]
    result <- Finverse3(grid =grid, qHT = qHT11, oHT = oHT11, nsim=500, itermax=20, lrw=160000, decay=1.4,ncore=10,
                        geo = geo)
    lire[[i]] = result
  }
  re = list(trueK=trueK,qHT1=qHT1,oHT1=oHT1,lire=lire,grid=grid)
  ### save the result
  saveRDS(re, file = paste0(iseed,"/lire_syn.rds"))
  re = readRDS(file = paste0(iseed,"/lire_syn.rds"))
  lire = re$lire
  trueK = re$trueK
  oHT1 =re$oHT1
  qHT1 =re$qHT1
  p1k = list()
  p2k = list()
  p3k = list()
  p4k = list()
  i1 = paste0("(",letters[1:4],") " ,1:4,"次抽水试验")
  inn = c(25,25,25,25)
  for(i in 1:4){
    niterm = length(lire[[i]])
    #niterm = inn[i]
    ii = se[1:i]
    #ii =1:i
    oHT11 = oHT1[ii]
    tmp = inversePlot(niterm=niterm,grid=grid,iterdf = lire[[i]],oHT = oHT11,trueK=trueK,
                      p1title = i1[i],
                      p1z = "ln(T) [L2/T]",
                      p1xlab = "X/m",
                      p1ylab = "Y/m",
                      p2title = i1[i],
                      p2z = "Var[ln(T)] [L2/T]",
                      p2xlab = "X/m",
                      p2ylab = "Y/m",
                      p3title = i1[i],
                      p3xlab = "模拟降深/m",
                      p3ylab = "观测降深/m",
                      p4title = i1[i],
                      p4xlab = "参考lnT",
                      p4ylab = "估计lnT")
    p1k[[i]] = tmp$lnk
    p2k[[i]] = tmp$varlnk
    p3k[[i]] = tmp$headscatter
    p4k[[i]] = tmp$lnkscatter
  }
  ncol =2
  nrow = 2
  p1 = ggpubr::ggarrange(plotlist = p1k, ncol = ncol,nrow =nrow)
  ggsave(filename = paste0(iseed,"/figure2_lnT.pdf"),plot = p1,
         width = 10, height = 8)
  p2 = ggpubr::ggarrange(plotlist = p2k, ncol = ncol,nrow =nrow)
  ggsave(filename = paste0(iseed,"/figure3_varlnT.pdf"),plot = p2,
         width = 10, height = 8)
  p3 = ggpubr::ggarrange(plotlist = p3k, ncol = ncol,nrow =nrow)
  ggsave(filename = paste0(iseed,"/figure4_headscatter.pdf"),plot = p3,
         width = 10, height = 8)
  p4 = ggpubr::ggarrange(plotlist = p4k, ncol = ncol,nrow =nrow)
  ggsave(filename = paste0(iseed,"/figure5_lnkscatter.pdf"),plot = p4,
         width = 10, height = 8)

  # forward prediction.

  ### now we use that to do the prediction on test 5-6
  plist = list()
  # test 5-9 -se
  qHT_pred = qHT1[-se]
  oHT_pred = oHT1[-se]
  for(j in 1:4){
    pred_results = list()
    niterm = length(lire[[j]])
    #niterm = inn[j]
    invertedK = lire[[j]][[niterm]]$meanT
    TT_inv = exp(invertedK)
    for(i in 1:length(qHT_pred)){
      Qinf = qHT_pred[[i]]
      simh <- Fsteady2dsim(TT=TT_inv,Qinf=Qinf,grid=grid,lrw=160000)$solution
      oinf <- oHT_pred[[i]]
      pred_df = data.frame(observed = -oinf$data,
                           simulated = -samData(grid=grid,Oinf=oinf,h=simh)$data)
      pred_results[[i]] = pred_df
    }
    # make pred_results a long data frame
    library(dplyr)
    pred_results_df = dplyr::bind_rows(pred_results,.id='test_id')

    plist[[j]] = pred_results_df
  }

  p_pred = list()
  tmp = paste0("(",letters[1:4],") ")
  for(i in 1:4){
    pred_results_df = plist[[i]]
    p_pred[[i]] = predictPlot(df = pred_results_df,
                              title = paste0(tmp[i],"1-",i,"次试验结果预测"),
                              xlab = "观测降深/m",
                              ylab = "预测降深/m")
  }
  p_all = ggpubr::ggarrange(plotlist = p_pred, ncol=2, nrow=2)
  p_all
  ggsave(filename = paste0(iseed,"/figure6_prediction.pdf"),plot = p_all,
         width = 10, height = 8)

}

