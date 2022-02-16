#' Run PONDR VLXT predictor for each protein
#'
#' @param FASTA List of FASTA formatted sequences
#' @param ID Character vector of protein names
#'
#' @return
#' @import RSelenium
#' @export
#'
#' @examples
#' \dontrun{
#' PONDR_VLXT <- VLXT(FASTA,ID)
#' }
#'
VLXT <- function(FASTA, ID){

  PONDR_VLXT <<- vector(mode='list', length=length(ID))

  system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)
  rD <- rsDriver(browser="firefox", port=4545L, verbose=T)
  remDr <- rD[["client"]]

  for (i in 1:length(ID)) {

    remDr$navigate("http://www.pondr.com/")

    # Check Raw Output
    remDr$findElements("name", "wcwraw")[[1]]$clickElement()

    # Uncheck Graphic, Statistics, Sequence Report
    remDr$findElements("name", "graphic")[[1]]$clickElement()
    remDr$findElements("name", "seq")[[1]]$clickElement()
    remDr$findElements("name", "stats")[[1]]$clickElement()

    # Paste FASTA
    remDr$findElement(using = 'name', value = 'Sequence')$sendKeysToElement(as.list(FASTA[[i]][[1]][1]))
    remDr$findElement(using = 'name', value = 'Sequence')$sendKeysToElement(as.list("\n"))
    remDr$findElement(using = 'name', value = 'Sequence')$sendKeysToElement(as.list(FASTA[[i]][[1]][-1]))

    remDr$findElements("name", "submit_result")[[1]]$clickElement()

    Sys.sleep(2)

    rm(webElem)
    rm(remDr)

    remDr <- rD[["client"]]

    webElem <- remDr$findElement(value ='/html/body/pre[6]')

    data <- webElem$getElementText()[[1]]

    write.table(data,file="data.txt",sep="\n",col.names=F,row.names=F,quote=FALSE)

    PONDR <- data.frame(read.table("data.txt",header = T))

    PONDR_VLXT[i] <<- list(PONDR)

  }

  names(PONDR_VLXT) <<- ID
  return(PONDR_VLXT)

}
