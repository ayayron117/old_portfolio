#' @title chplot
#' @description Creates a charge-hydropathy plot for multiple proteins
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
#' chplot(s,n)

chplot <- function(sequences,names) {

  # Saves the input sequences into a data frame
  sequence <- data.frame(sequences)

  # Adds the protein names to the sequence data frame
  row.names(sequence) <- names

  # Data frame which contains the charges and hydropathies of each residue
  # (amino acid)
  charge_hydropathy <- data.frame(charge = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                             0, -1, 0, 0, 0, -1, 1, 0), IDP_hydropathy = c(0.4, 0.355, 1, 0.811, 0.777,
                                                                                           0.922, 0.966, 0.711, 0.111, 0.422, 0.7, 0, 0.455, 0.111, 0.111, 0.411, 0.144,
                                                                                           0.111, 0.066, 0.322))

  # Names of each residue
  row.names(charge_hydropathy) <- c('W','Y','I','F','C','L','V','M','N','T',
                                    'A','R','G','D','Q','S','H','E','K','P')

  # Counts the number of each residue in each protein
  count <- data.frame(W = str_count(sequence$s,'W'),Y = str_count(sequence$s,'Y'),
                      I = str_count(sequence$s,'I'),F = str_count(sequence$s,'F'),
                      C = str_count(sequence$s,'C'),L = str_count(sequence$s,'L'),
                      V = str_count(sequence$s,'V'),M = str_count(sequence$s,'M'),
                      N = str_count(sequence$s,'N'),T = str_count(sequence$s,'T'),
                      A = str_count(sequence$s,'A'),R = str_count(sequence$s,'R'),
                      G = str_count(sequence$s,'G'),D = str_count(sequence$s,'D'),
                      Q = str_count(sequence$s,'Q'),S = str_count(sequence$s,'S'),
                      H = str_count(sequence$s,'H'),E = str_count(sequence$s,'E'),
                      K = str_count(sequence$s,'K'),P = str_count(sequence$s,'P'))

  # Calculates the total hydropathy contributed by each residue
  H <- data.frame(W_H = count['W']*charge_hydropathy['W',"IDP_hydropathy"],
                  Y_H = count['Y']*charge_hydropathy['Y',"IDP_hydropathy"],
                  I_H = count['I']*charge_hydropathy['I',"IDP_hydropathy"],
                  F_H = count['F']*charge_hydropathy['F',"IDP_hydropathy"],
                  C_H = count['C']*charge_hydropathy['C',"IDP_hydropathy"],
                  L_H = count['L']*charge_hydropathy['L',"IDP_hydropathy"],
                  V_H = count['V']*charge_hydropathy['V',"IDP_hydropathy"],
                  M_H = count['M']*charge_hydropathy['M',"IDP_hydropathy"],
                  N_H = count['N']*charge_hydropathy['N',"IDP_hydropathy"],
                  T_H = count['T']*charge_hydropathy['T',"IDP_hydropathy"],
                  A_H = count['A']*charge_hydropathy['A',"IDP_hydropathy"],
                  R_H = count['R']*charge_hydropathy['R',"IDP_hydropathy"],
                  G_H = count['G']*charge_hydropathy['G',"IDP_hydropathy"],
                  D_H = count['D']*charge_hydropathy['D',"IDP_hydropathy"],
                  Q_H = count['Q']*charge_hydropathy['Q',"IDP_hydropathy"],
                  S_H = count['S']*charge_hydropathy['S',"IDP_hydropathy"],
                  H_H = count['H']*charge_hydropathy['H',"IDP_hydropathy"],
                  E_H = count['E']*charge_hydropathy['E',"IDP_hydropathy"],
                  K_H = count['K']*charge_hydropathy['K',"IDP_hydropathy"],
                  P_H = count['P']*charge_hydropathy['P',"IDP_hydropathy"])

  # Counts the total number of residues in each protein
  n_s <- nchar(sequences)

  # Calculates the net hydropathy of each protein
  s_H <- rowSums(H)

  # Calculates the mean net hydropathy of each protein
  mean_H <- s_H/n_s

  # Calculates the total charge contributed by each residue
  C <- data.frame(W_C = count['W']*charge_hydropathy['W',"charge"],
                  Y_C = count['Y']*charge_hydropathy['Y',"charge"],
                  I_C = count['I']*charge_hydropathy['I',"charge"],
                  F_C= count['F']*charge_hydropathy['F',"charge"],
                  C_C= count['C']*charge_hydropathy['C',"charge"],
                  L_C= count['L']*charge_hydropathy['L',"charge"],
                  V_C= count['V']*charge_hydropathy['V',"charge"],
                  M_C= count['M']*charge_hydropathy['M',"charge"],
                  N_C= count['N']*charge_hydropathy['N',"charge"],
                  T_C= count['T']*charge_hydropathy['T',"charge"],
                  A_C= count['A']*charge_hydropathy['A',"charge"],
                  R_C= count['R']*charge_hydropathy['R',"charge"],
                  G_C= count['G']*charge_hydropathy['G',"charge"],
                  D_C= count['D']*charge_hydropathy['D',"charge"],
                  Q_C= count['Q']*charge_hydropathy['Q',"charge"],
                  S_C= count['S']*charge_hydropathy['S',"charge"],
                  H_C= count['H']*charge_hydropathy['H',"charge"],
                  E_C= count['E']*charge_hydropathy['E',"charge"],
                  K_C= count['K']*charge_hydropathy['K',"charge"],
                  P_C= count['P']*charge_hydropathy['P',"charge"])

  # Calculates the total charge of each protein
  s_C <- rowSums(C)

  # Calculates the mean net charge of each protein
  mean_C <- s_C/n_s

  # Takes the absolute value of the mean net charge
  mean_C <- abs(mean_C)

  # Stores the values of the mean hydropathy and mean charge in a data frame
  d <- data.frame(mean_H,mean_C,Proteins=names,line_color='black')

  plot <- ggplot(d,aes(mean_H,mean_C,colour=Proteins))+
    ggtitle("Charge VS Hydropathy")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Mean Scaled Hydropathy") +
    ylab("Absolute Net Charge")+
    coord_cartesian(xlim = c(0.2,0.65), ylim = c(0,0.65))+
    geom_point()+
    geom_segment(aes(x = c(0.4043), y = c(0), xend = c(0.6777), yend = c(0.75)),            color='black')+
    theme(panel.background = element_rect(fill = "aliceblue",color = "aliceblue"),
          title =element_text(size=12,face='bold'))+
    geom_point(size=3)

  dir <- getwd()
  dir <- paste0(dir,'/chplot.png')
  ggsave(dir, width= 9, height= 5)

  print(plot)

}
