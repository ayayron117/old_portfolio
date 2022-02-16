#' Calculate Cumulative Distribution Frequency
#'
#' @param PONDR_VLXT List of disorder scores for each protein
#' @param ID Chracter vector of protein names
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' CDF <- get_CDF(PONDR_VLXT, PROTEINS$ID)
#' }
get_CDF <- function(PONDR_VLXT,ID){
  CDF_ <- {}
  CDF <- vector(mode='list', length=length(ID))
  names(CDF) <- ID
  for (i in 1:length(PONDR_VLXT)){
    P <-as.numeric(unlist(PONDR_VLXT[[i]][3]))
    breaks <- seq(0,1,by=0.05)
    cut <- cut(P,breaks,right=F)
    freq <- table(cut)
    cumfreq <- cumsum(freq)/length(as.numeric(unlist(PONDR_VLXT[[i]][3])))
    Frequency = c(0, cumfreq)
    cumfreq<- as.data.frame(Frequency)
    CDF_ <- cumfreq
    # CDF_$Intervals <- rownames(CDF_)
    row.names(CDF_) <- seq(from=1,to=21,by=1)
    CDF_ -> CDF[[i]]
  }
  return(CDF)
}
