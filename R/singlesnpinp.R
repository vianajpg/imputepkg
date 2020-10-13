#' Impute the Genotypes per SNP
#'
#' This function impute data in the offspring based on the parents genotypes
#'
#' @param OS Offspring genotypic data
#' @param P Parents genotypic data
#' @return Imputed Offspring Data
#' @export
singlesnpinp <- function(OS = OS, P = P) {
  nsnps_parents <- nrow(x = P)
  nvar_offspring <- ncol(x = OS)
  newOS <- data.frame(matrix(data = NA, nrow = nsnps_parents, ncol = nvar_offspring))
  newOS[,1:4] <- P[,1:4]
  P_OS_snpidx <- P[,2] %in% OS[,2]
  newOS[P_OS_snpidx,5:nvar_offspring] <- OS[,5:nvar_offspring]
  for(i in seq(from = 1, to = nsnps_parents)){
    for(j in seq(from = 5, to = nvar_offspring)){
      if(is.na(newOS[i,j]) & P[i,5] == 0 & P[i,6] == 0) {
        newOS[i,j] <- 0
      } else {
        if(is.na(newOS[i,j]) & P[i,5] == 2 & P[i,6] == 2) {
          newOS[i,j] <- 2
        } else {
          if(is.na(newOS[i,j]) & P[i,5] == 0 & P[i,6] == 2) {
            newOS[i,j] <- 1
          } else {
            if(is.na(newOS[i,j]) & P[i,5] == 2 & P[i,6] == 0) {
              newOS[i,j] <- 1
            } else {
              newOS[i,j] <- 3
            }
          }
        }
      }
    }
  }
  return(newOS)
}
