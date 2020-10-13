#' Simulate Parents Genotypes
#'
#' This function simulate a genotype data for one chromosome in two parents P1 and P2
#'
#' @param nsnps Number of requested SNPs
#' @param chrname Name for the simulated chromosome
#' @param minpos Initial position for the simulated chromosome
#' @param maxpos Final position for the simulated chromosome
#' @param Hom1 Proportion of REF homozygotes
#' @param Hom2 Proportion of ALT homozygotes
#' @param seed Seed for the sampling steps
#' @return A data frame with the simulated genotypes and genomic information
#' @export
parentsim <- function(nsnps = nsnps, chrname = chrname, minpos = 1, maxpos = 500000, Hom1 = 0.25, Hom2 = 0.25, seed = 123456){
  set.seed(seed = seed)
  Het <- 1-Hom1-Hom2
  P1 <- sample(x = c(0,1,2), size = nsnps, replace = TRUE, prob = c(Hom1, Het, Hom2))
  P2 <- sample(x = c(0,1,2), size = nsnps, replace = TRUE, prob = c(Hom1, Het, Hom2))
  CHROM <- rep(chrname, nsnps)
  POS <- sort(x = sample(x = seq(from = minpos, to = maxpos), size = nsnps, replace = FALSE))
  SNPID <- paste0("SNP", seq(from = 1, to = nsnps))
  ALLELES <- c()
  for(i in seq(from = 1, to = nsnps)){
    ALLELES[i] <- paste(sample(x = c("A","T","G","C"), size = 2, replace = FALSE), collapse = "/")
  }
  parentmat <- data.frame(CHROM,POS,SNPID,ALLELES,P1,P2)
  return(parentmat)
}
