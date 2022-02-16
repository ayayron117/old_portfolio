#' Delta Cumulative Distribution Frequency
#'
#' @param CDF List of CDF values for each protein
#'
#' @return
#' @export
#'
#' @examples
#'\dontrun{
#' d_CDF <- delta_CDF(CDF)
#' }
delta_CDF <- function(CDF){
  delta_CDF <- data.frame(matrix(nrow=length(PROTEINS$ID),ncol=1))
  colnames(delta_CDF) <- c("delta_CDF")
  rownames(delta_CDF) <- PROTEINS$ID
  for (i in 1:length(PROTEINS$ID)){
    breaks <- seq(0,1,by=0.05)
    diff <- breaks[13:19]-CDF[[i]]$Frequency[13:19]
    sum <- sum(diff)
    delta_CDF$delta_CDF[i] <- sum/7
  }
  return(delta_CDF)
}
