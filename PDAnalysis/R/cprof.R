#' Composition Profile
#'
#' @param sequence Character vector of sequences, do not use FASTA format
#' @param ID Character vector containing the names of the proteins
#'
#' @return
#' @import purrr ggplot2 stringr
#' @export
#'
#' @examples
#' \dontrun{
#' cprof(PROTEINS$SEQUENCE, PROTEINS$ID)
#' }
cprof <- function(sequence,ID) {

  # Saves the input sequences into a data frame
  seq <- data.frame(sequence)

  # Adds the protein names to the sequence data frame
  row.names(seq) <- ID

  # Data frame which contains the reference data sets. Each data set contains
  # the average frequency of each residue (amino acid) in a set proteins.
  # SwissProt is a set which attempts to represent the frequency of residues
  # in nature, PDB_S25 contains crystalizable proteins, and DisProt contains
  # intrinically disordered proteins.
  profile_ref <- data.frame(SwissProt = c(1.5, 1.13, 5.9, 3.03, 3.96, 9.65, 2.29,
                                          6.73, 4.13, 2.38, 5.4, 5.41, 5.35, 6.96, 7.89, 5.92, 3.95, 6.83, 6.67, 4.83),
                            PDB_S25= c(1.74, 1.44, 5.61, 3.5, 3.98, 8.68, 2.41, 6.72, 4.58, 2.22, 4.93, 5.63,
                                       5.83, 7.16, 7.7, 6.37, 3.95, 6.19, 6.65, 4.57),
                            DisProt = c(0.8, 0.67, 3.24, 2.13, 2.44, 6.22, 1.93, 5.41, 3.82, 1.87, 4.82, 5.56,
                                        5.8, 7.41, 8.1,7.85, 5.27, 8.65, 9.89, 8.11))

  # Names of each residue
  row.names(profile_ref) <- c('C','W','I','Y','F','L','H','V','N','M','R','T',
                              'D','G', 'A','K','Q','S','E','P')

  # Counts the number of each residue in each protein
  count <- data.frame(W = str_count(seq,'W'),Y = str_count(seq,'Y'),
                      I = str_count(seq,'I'),F = str_count(seq,'F'),
                      C = str_count(seq,'C'),L = str_count(seq,'L'),
                      V = str_count(seq,'V'),M = str_count(seq,'M'),
                      N = str_count(seq,'N'),T = str_count(seq,'T'),
                      A = str_count(seq,'A'),R = str_count(seq,'R'),
                      G = str_count(seq,'G'),D = str_count(seq,'D'),
                      Q = str_count(seq,'Q'),S = str_count(seq,'S'),
                      H = str_count(seq,'H'),E = str_count(seq,'E'),
                      K = str_count(seq,'K'),P = str_count(seq,'P'))

  # Counts the total number of residues in each protein
  n_s <- nchar(seq)

  # Calculates the frequency of each residue in each protein
  freq <- count/n_s

  # Calculates the percent frequency
  freq <- freq*100

  Swiss_prof <- list()
  PDB_S25_prof <- list()
  DisProt_prof <- list()
  prof <- vector(mode = 'list', length = length(seq))

  # Calculates the relative frequency of each residue based on each data set and
  # stores them in lists
  for (i in 1:length(seq)) {

    Swiss_prof[[i]] <- (freq[i,] - profile_ref[,1])/profile_ref[,1]

    PDB_S25_prof[[i]] <- (freq[i,] - profile_ref[,2])/profile_ref[,2]

    DisProt_prof[[i]] <- (freq[i,] - profile_ref[,3])/profile_ref[,3]

    prof[[i]][[1]] <- Swiss_prof[[i]]
    prof[[i]][[2]] <- PDB_S25_prof[[i]]
    prof[[i]][[3]] <- PDB_S25_prof[[i]]

  }

  for (i in 1:length(seq)) {

    # Extracts data from the lists and stores them in arrays
    array1 <- array(as.numeric(unlist(prof[[i]][[1]])), dim = c(20,1))
    array2 <- array(as.numeric(unlist(prof[[i]][[2]])), dim = c(20,1))
    array3 <- array(as.numeric(unlist(prof[[i]][[3]])), dim = c(20,1))

    # Combines them into one array
    prof_array <- rbind(array1,array2,array3)

    # Stores the names of the residues in an object
    residues <- row.names(profile_ref)

    # Stores the names of the reference data sets in an object, each is repeated 20 times
    sets <- rep(c('SwissProt','PDB_S25','DisProt'),each=20)

    # Creates a data frame that can be used to generate a grouped barplot
    profile <- data.frame(freq=prof_array,residues=residues,Set=sets)

    # Classifies the names of the residues as factors so that to prevent ggplot from alphabetizing them
    profile$residues <- factor(profile$residues, levels = row.names(profile_ref))

    # Stores the name of the protein for the title during each iteration
    title <- paste(ID[i],'Compostion Profile',sep=" ")

    # Plot
    plot<-ggplot(profile, aes(x=residues, y=freq, fill=Set))+
      geom_bar(stat="identity", position=position_dodge(),color='black')+
      theme(panel.background = element_rect(fill = "aliceblue",color = "aliceblue"),
            panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "gray89"),
            panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "gray89"))+
      xlab('Amino Acid')+
      ylab('Relative Frequency')+
      ylim(-1.5,2)+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5),title =element_text(size=12,face='bold'))

    # # Saves the plot as an image file, it can be found in the working directory
    # n <- as.character(i)
    # dir <- getwd()
    # dir <- paste0(dir,'/cprof_',n,'.png')
    # ggsave(dir, width= 9, height= 5)
    #
     print(plot)

  }

}
