#' Create CDF plot for multiple proteins
#'
#' @param CDF List of CDF values
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' plot_CDF(CDF)
#' }
plot_CDF <- function(CDF){

  names <- names(CDF)

  for (i in 1:length(CDF)){
    plot <-  plot(seq(0,1,by=0.05),CDF[[i]]$Frequency,main=names[i], xlab="PONDR Score",ylab="Cumulative Frequency")

    print(plot)

  }

}
