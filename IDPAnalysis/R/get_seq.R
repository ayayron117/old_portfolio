#'Get sequences from UniProt, also stores them in FASTA format
#'
#' @param ID Character vector of UniProt IDs
#'
#' @return
#' @import UniprotR
#' @export
#'
#' @examples
#' \dontrun{
#' PROTEINS <- get_seq(ID)
#' }
get_seq <- function(ID) {

  ProtNames <- GetNamesTaxa(ID[[1]])
  ProtSeqs <- GetSequences(ID[[1]])

  ProtNames$ID <- row.names(ProtNames)
  ProtSeqs$ID <- row.names(ProtSeqs)

  old.column.names <- c("ID", "Entry.name")
  ProtNames_cleaned <- ProtNames[old.column.names]

  column.names <- c("ID", "Sequence")
  ProtSeqs.cleaned <- ProtSeqs[column.names]

  PROTEINS <- merge(ProtNames_cleaned,ProtSeqs.cleaned,sort=F)

  colnames(PROTEINS) <- c("ID","NAME","SEQUENCE")

  fasta = {}

  imported.fastas <- {}

  get_fastas <- function(i) {
    pb <- progress::progress_bar$new(total = length(i))
    imported.fastas[i] <- read.csv(paste0("https://www.uniprot.org/uniprot/",i,".Fasta"),header= F, sep= "\t")
    return(imported.fastas)
    pb$tick()
  }

  skip_errors <- function(i) return(tryCatch(get_fastas(i), error= function(e) e))

  test <- lapply(ID[[1]], skip_errors)

  PROTEINS$FASTA <- as.list(test)

  PROTEINS <<- PROTEINS

}
