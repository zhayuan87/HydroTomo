#' Plot the generated head distribution from Fsteady2dsim.
#' @return a plot showing contour.
#' @param steady2d a data.frame, the result of Fsteady2dsim (x,y,headchange).
#' @param plotshow a logical value, if TRUE, the plot will be shown.
#' @param plotfile a string, the path to save the plot.
#' @param plotwidth a numeric value, the width of the plot.
#' @param plotheight a numeric value, the height of the plot.
#' @param ifdrawdown a logical value, if true,  use drawdown, else using head change.
#' @param palette a string, the palette to use. could be magma, inferno, plasma, viridis, cividis, rocket, marko, turbo, or simply A-G.
#' @export

Plotsteady2d <- function(s, palette="viridis",plotshow = TRUE,
                         plotfile = NULL, plotwidth = 10, plotheight = 10, ifdrawdown = TRUE) {

  if (ifdrawdown) {
    s$solution = - s$solution
    title = "2D Steady-State Drawdown Distribution"
    z = "Drawdown"
  } else {
    title = "2D Steady-State Head Distribution"
    z = "Head"
  }

  require('ggplot2')
  p = ggplot(s, aes(x, y, fill = solution)) +
    geom_tile() +
    scale_fill_viridis_c(option = palette)
    labs(title = title, x = "X/m", y = "Y/m", fill = z)



  if(!is.null(plotfile)){
    ggsave(plotfile, width = plotwidth, height = plotheight)
  }
  return(p)
}



#' now we can choose the best iteration to plot.
#' estimated K plot.
#' @export
#' @param niter the number of iterations.
#' @param iterdf the data.frame of iterations.
#' @param trueK the true K field (used in synthetic case).
#' @examples
#' set.seed(100)
#' trueK <- random2d(nsim=1)
#' TT <- trueK[,-c(1,2)]
#' trueh <- Fsteady2dsim(TT=TT)
#' trueh <- trueh$solution
#' n <- 1600
#' nobs <- 25
#' locx = c(15,18,22,25,30)
#' locy = c(15,18,22,25,30)
#' loc= expand.grid(x=locx,y=locy)
#' loc_obs <- (loc$y-1)*40 + loc$x
#' trueobsh <- trueh[loc_obs]
#' loc_obs = list(loc_obs)
#' trueobsh = list(trueobsh)
#' result <- Finverse(loc_obs=loc_obs,trueobsh=trueobsh)
#' inversePlot(niterm=5,iterdf = result,trueh = trueobsh,trueK=trueK)
inversePlot <- function(niterm=1,iterdf,trueh,trueK=NULL){
  # niterm = 5
  # note that trueh is actually a list contains many pumping tests.
  np <- sapply(trueh,length) # get the number of observations for each test.
  ntest <- rep("test1",np[1])
  if(length(np)>=2){
    for (i in 2:length(np)){
      ntest <- append(ntest,rep(paste0("test",i),np[i]))
    }
  }

  trueh <- unlist(trueh) # unlist
  df <- data.frame(x=trueK$x,y=trueK$y,v=iterdf[[niterm]]$meanT,var=iterdf[[niterm]]$varT)
  library(ggplot2)
  p1<- ggplot(df) +
    aes(x = x, y = y, fill = v) +
    geom_tile()+
    scale_fill_viridis_c(option = "viridis", direction = 1)+
    ggtitle(label="estimated lnK value")
  p2 <- ggplot(df) +
    aes(x = x, y = y, fill = var) +
    geom_tile()+
    scale_fill_viridis_c(option = "viridis", direction = 1)+
    ggtitle(label="estimated lnK variance")

  dfh <- data.frame(trueh=trueh,simh=iterdf[[niterm]]$meanobsh,varh=iterdf[[niterm]]$varobsh,ntest=ntest)
  #https://www.roelpeters.be/how-to-add-a-regression-equation-and-r-squared-in-ggplot2/
  #https://rpkgs.datanovia.com/ggpubr/reference/stat_regline_equation.html
  require(ggpubr)
  p4 <- ggplot(dfh) +
    geom_errorbar(aes(x = trueh,ymin=simh-varh^0.5,ymax=simh+varh^0.5),width=0.02,color="gray")+
    geom_point(size=4,aes(x = trueh, y = simh,colour=ntest))+
    geom_smooth(method = "lm", aes(x = trueh, y = simh),se=FALSE) +
    stat_regline_equation(aes(x = trueh, y = simh),label.y = -0.5) +
    stat_cor( aes(x = trueh, y = simh),label.y = -1) +
    ggtitle(label="observation versus simulated head with variance")

  p1
  p2
  p4

  if(!is.null(trueK)){
    scatterks <- data.frame(x=log(trueK$Tp),y=log(iterdf[[niterm]]$meanT))
    p3 <- ggplot(scatterks) +
      aes(x,y) +
      geom_point(size=1)+
      geom_smooth(method = "lm", se=T) +
      stat_regline_equation(label.y = -0.5)+
      ggtitle(label="true versus inverted lnT")
  }
  if(!is.null(trueK))
    return(list(lnk=p1,varlnk=p2,lnkscatter=p3,headscatter=p4))
  else
    return(list(lnk=p1,varlnk=p2,headscatter=p4))
}
