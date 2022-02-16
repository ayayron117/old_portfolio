#' Creat a CH-CDF plot for multiple proteins
#'
#' @param delta_CDF list of delta CDF values
#' @param delta_CH List of delta CH values
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' plot_CH_CDF(delta_CDF,delta_CH)
#' }
plot_CH_CDF <- function(delta_CDF,delta_CH){

  CH_CDF <- data.frame(delta_CDF=delta_CDF,delta_CH=delta_CH)
  plot <- ggplot(CH_CDF, aes(x=delta_CDF, y= delta_CH)) +
    geom_point(shape=1, size =2) +
    lims(x=c(-1,1),y=c(-1,1)) +
    theme_minimal() +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

  print(plot)

}
