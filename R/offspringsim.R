#' Simulate Offspring Genotpyes
#'
#' This function simulate a genotype data for one chromosome for a given number of samples based on parents genotypes
#'
#' @param parentdf Data frame containing the parents genotype - generated from parentsim function
#' @param prop Proportion of original number os SNPs retained in the final simulated genotype - subset function
#' @param nsamples Number of offsprings
#' @return A data frame with the simulated genotypes and genomic information for the offsprings
#' @export
offspringsim <- function(parentdf = parentdf, prop =  prop, nsamples = nsamples){
  nsnps <- nrow(x = parentdf)
  nsubsnps <- round(x = nsnps * prop, digits = 0)
  subsnpsidx <- sort(x = sample(x = seq(from = 1, to = nsnps), size = nsubsnps, replace = FALSE))
  offspringmat <- data.frame(matrix(data = NA, nrow = nsubsnps, ncol = nsamples+6))
  colnames(offspringmat) <- c("CHROM", "POS", "SNPID", "ALLELES", paste0("OS", seq(from = 1, to = nsamples)))
  offspringmat[,1:4] <- parentdf[subsnpsidx,1:4]
  for(i in seq(from = 1, to = nsubsnps)){
    for(j in seq(from = 5, to = nrow(x = offspringmat))){
      GP1 <- parentdf[subsnpsidx[i],5]
      GP2 <- parentdf[subsnpsidx[i],6]
      if(GP1 == 0 & GP2 == 0){
        offspringmat[i,j] <- 0
      } else {
        if(GP1 == 0 & GP2 == 1){
          offspringmat[i,j] <- sample(x = c(0,1), size = 1, replace = FALSE, prob = c(0.5, 0.5))
        } else {
          if(GP1 == 0 & GP2 == 2){
            offspringmat[i,j] <- 1
          } else {
            if(GP1 == 1 & GP2 == 0){
              offspringmat[i,j] <- sample(x = c(1,0), size = 1, replace = FALSE, prob = c(0.5, 0.5))
            } else {
              if(GP1 == 1 & GP2 == 1){
                offspringmat[i,j] <- sample(x = c(0,1,2), size = 1, replace = FALSE, prob = c(0.25, 0.5, 0.25))
              } else {
                if(GP1 == 1 & GP2 == 2){
                  offspringmat[i,j] <- sample(x = c(1,2), size = 1, replace = FALSE, prob = c(0.5, 0.5))
                } else {
                  if(GP1 == 2 & GP2 == 0){
                    offspringmat[i,j] <- 1
                  } else {
                    if(GP1 == 2 & GP2 == 1){
                      offspringmat[i,j] <- sample(x = c(2,1), size = 1, replace = FALSE, prob = c(0.5, 0.5))
                    } else {
                      if(GP1 == 2 & GP2 == 2){
                        offspringmat[i,j] <- 2
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(offspringmat)
}
