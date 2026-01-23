library(HydroTomo)
domain = c(41,41,0,41,0,41)
grid = GenGrid(domain=domain)
set.seed(1999)
geo = list(me = 3.2, var = 3, geomod = "Sph", anis = c(90, 1), range = 10, nugget = 0)
trueK <- random2d(nsim=1,grid=grid,geo = geo)
# plot the generated log -T field.
Plotparameter2d(trueK,iflog=T)
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
mean_oHT1 = mean(dfoHT1$data)
print(paste("The mean value of oHT1 is ", round(mean_oHT1,4)))



#### now do the inverse.

lire = list()
# we select qHT1 and oHT1 from blh_synthetic, first test1, then 1,2, then 1,2,3
for(i in 1:9){
  qHT11 = qHT1[1:i]
  oHT11 = oHT1[1:i]
  result <- Finverse3(grid =grid, qHT = qHT11, oHT = oHT11, nsim=50, itermax=20, lrw=160000, decay=2,ncore=5,
                      geo = geo)
  lire[[i]] = result
}
re = list(trueK=trueK,qHT1=qHT1,oHT1=oHT1,lire=lire,grid=grid)
### save the result
saveRDS(re, file = "lire_syn_20260119-gau.rds")


re = readRDS(file = "lire_syn_20260119.rds")
lire = re$lire
trueK = re$trueK
oHT1 =re$oHT1
qHT1 =re$qHT1
p1k = list()
p2k = list()
p3k = list()
p4k = list()
i1 = paste0("(",letters[1:9],") " ,1:9,"次抽水试验")
for(i in 1:9){
  niterm = length(lire[[i]])
  oHT11 = oHT1[1:i]
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
p1 = ggpubr::ggarrange(plotlist = p1k, ncol = 3, nrow = 3)
p1
p2 = ggpubr::ggarrange(plotlist = p2k, ncol = 3, nrow = 3)
p2
p3 = ggpubr::ggarrange(plotlist = p3k, ncol = 3, nrow = 3)
p3
p4 = ggpubr::ggarrange(plotlist = p4k, ncol = 3, nrow = 3)
p4

### we also want to
set.seed(401)
oHT_noised2 <- oHT1
for (i in 1:length(oHT_noised2)){
  noise = rnorm(n=nrow(oHT_noised2[[i]]), mean=0, sd=0.05)
  oHT_noised2[[i]]$data <- oHT_noised2[[i]]$data *(1+noise) # relative error.
}

# we select qHT1 and oHT1 from blh_synthetic, first test1, then 1,2, then 1,2,3
for(i in 1:9){
  qHT11 = qHT1[1:i]
  oHT11 = oHT_noised2[1:i]
  result <- Finverse3(grid =grid, qHT = qHT11, oHT = oHT11, nsim=50, itermax=20, lrw=160000, decay=2,ncore=5,
                      geo = geo)
  lire[[i]] = result
}
re = list(trueK=trueK,qHT1=qHT1,oHT1=oHT_noised2,lire=lire,grid=grid)

saveRDS(re, file = "lire_syn_20260119-noise0.2.rds")
re = readRDS(file = "lire_syn_20260119-noise0.2.rds")
lire = re$lire
trueK = re$trueK
oHT_noised2 =re$oHT1
qHT1 =re$qHT1
p1k = list()
p2k = list()
p3k = list()
p4k = list()
i1 = paste0("(",letters[1:9],") " ,1:9,"次抽水试验")
for(i in 1:9){
  niterm = length(lire[[i]])
  oHT11 = oHT_noised2[1:i]
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
p1 = ggpubr::ggarrange(plotlist = p1k, ncol = 3, nrow = 3)
p1
p2 = ggpubr::ggarrange(plotlist = p2k, ncol = 3, nrow = 3)
p2
p3 = ggpubr::ggarrange(plotlist = p3k, ncol = 3, nrow = 3)
p3
p4 = ggpubr::ggarrange(plotlist = p4k, ncol = 3, nrow = 3)
p4
