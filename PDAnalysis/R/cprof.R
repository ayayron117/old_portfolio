#' @title cprof
#' @description Plots the residue profile based on reference data sets
#' @param sequences Character vector of sequences
#' @param names Character vector of protein names
#' @examples
#'
#' s <- c('MGNHAGKRELNAEKASTNSETNRGESEKKRNLGELSRTTSEDNEVFGEADANQNNGTSSQ
#' DTAVTDSKRTADPKNAWQDAHPADPGSRPHLIRLFSRDAPGREDNTFKDRPSESDELQTI
#' QEDSAATSESLDVMASQKRPSQRHGSKYLATASTMDHARHGFLPRHRDTGILDSIGRFFG
#' GDRGAPKRGSGKDSHHPARTAHYGSLPQKSHGRTQDENPVVHFFKNIVTPRTPPPSQGKG
#' RGLSLSRFSWGAEGQRPGFGYGGRASDYKSAHKGFKGVDAQGTLSKIFKLGGRDSRSGSP
#' MARR')
#'
#' n <- c('MBP P02686')
#'
#' cprof(s,n)

cprof <- function(sequences,names) {

  # Saves the input sequences into a data frame
  sequence <- data.frame(sequences)

  # Adds the protein names to the sequence data frame
  row.names(sequence) <- names

  # Data frame which contains the reference data sets. Each data set contains
  # the average frequency of each residue (amino acid) in a set proteins.
  # SwissProt is a set which attempts to represent the frequency of residues
  # in nature, PDB_S25 contains crystalizable proteins, and DisProt contains
  # intrinically disordered proteins.
  profile_ref <- data.frame(SwissProt = c(1.5, 1.13, 5.9, 3.03, 3.96, 9.65, 2.29,
  6.73, 4.13, 2.38, 5.4, 5.41, 5.35, 6.96, 7.89, 5.92, 3.95, 6.83, 6.67, 4.83),
  PDB_S25= c(1.74, 1.44, 5.61, 3.5, 3.98, 8.68, 2.41, 6.72, 4.58, 2.22, 4.93,
  5.63, 5.83, 7.16, 7.7, 6.37, 3.95, 6.19, 6.65, 4.57), DisProt = c(0.8, 0.67,
  3.24, 2.13, 2.44, 6.22, 1.93, 5.41, 3.82, 1.87, 4.82, 5.56, 5.8, 7.41, 8.1,
  7.85, 5.27, 8.65, 9.89, 8.11))

  # Names of each residue
  row.names(profile_ref) <- c('C','W','I','Y','F','L','H','V','N','M','R','T',
                              'D','G', 'A','K','Q','S','E','P')

  # Counts the number of each residue in each protein
  count <- data.frame(C = str_count(sequence$s,'C'),W = str_count(sequence$s,'W'),
                      I = str_count(sequence$s,'I'),Y = str_count(sequence$s,'Y'),
                      F = str_count(sequence$s,'F'),L = str_count(sequence$s,'L'),
                      H = str_count(sequence$s,'H'),V = str_count(sequence$s,'V'),
                      N = str_count(sequence$s,'N'),M = str_count(sequence$s,'M'),
                      R = str_count(sequence$s,'R'),T = str_count(sequence$s,'T'),
                      D = str_count(sequence$s,'D'),G = str_count(sequence$s,'G'),
                      A = str_count(sequence$s,'A'),K = str_count(sequence$s,'K'),
                      Q = str_count(sequence$s,'Q'),S = str_count(sequence$s,'S'),
                      E = str_count(sequence$s,'K'),P = str_count(sequence$s,'P'))

  # Counts the total number of residues in each protein
  n_s <- nchar(sequences)

  # Calculates the frequency of each residue in each protein
  freq <- count/n_s

  # Calculates the percent frequency
  freq <- freq*100

  Swiss_prof <- list()
  PDB_S25_prof <- list()
  DisProt_prof <- list()
  prof <- vector(mode = 'list', length = length(s))

  # Calculates the relative frequency of each residue based on each data set and
  # stores them in lists
  for (i in 1:length(s)) {

    Swiss_prof[[i]] <- (freq[i,] - profile_ref[,1])/profile_ref[,1]

    PDB_S25_prof[[i]] <- (freq[i,] - profile_ref[,2])/profile_ref[,2]

    DisProt_prof[[i]] <- (freq[i,] - profile_ref[,3])/profile_ref[,3]

    prof[[i]][[1]] <- Swiss_prof[[i]]
    prof[[i]][[2]] <- PDB_S25_prof[[i]]
    prof[[i]][[3]] <- PDB_S25_prof[[i]]

  }

  for (i in 1:length(s)) {

    # Extracts data from the lists and stores them in arrays
    array1 <- array(as.numeric(unlist(prof[[i]][[1]])), dim = c(20,1))
    array2 <- array(as.numeric(unlist(prof[[i]][[2]])), dim = c(20,1))
    array3 <- array(as.numeric(unlist(prof[[i]][[3]])), dim = c(20,1))

    prof_array <- cbind(array1,array2,array3)

    row.names(prof_array) <- c('C','W','I','Y','F','L','H','V','N','M','R','T',
                               'D','G', 'A','K','Q','S','E','P')

    colnames(prof_array) <- c('SwissProt','PDB_S25','DisProt')

    t <- aperm(prof_array)

    title <- paste(names[i],'Compostion Profile',sep=" ")

    # Plot
    par(mar = c(5, 5, 2, 2),xpd = TRUE)

    barplot(t,beside = TRUE,ylim = c(-1,2),col=c("red","blue","green"),xlab='Amino Acid',
            ylab='Percent Frequency',main=title,cex.lab =1.2)

    legend('bottomleft',legend=c('SwissProt','PDB_S25','DisProt'),fill=c("red","blue","green"),
           cex=0.8, ncol=3,inset=c(-0.13,-0.265))

  }

}
