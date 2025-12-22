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
