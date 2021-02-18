#' Impute Inbred VCF per defined window
#'
#' This function impute VCF data from the inbred populations based on the VCF data from the parents
#'
#' @param P_in_vcf Parents original VCF
#' @param IP_in_vcf Inbred Populations original VCF
#' @param IP_new_out_vcf Name for an output - Imputed Inbred VCF
#' @param win Size of the discrete window
#' @return Imputed inbred populations VCF
#' @export

imputeVCF_haplo <- function(P_in_vcf = "", IP_in_vcf = "", IP_new_out_vcf = "", win = 5){
  wd <- getwd()
  if(!require(BiocManager)){
    install.packages("BiocManager")
  }
  if(!require(SNPRelate)){
    require(BiocManager)
    BiocManager::install(pkgs = "SNPRelate")
  }
  if(!require(SeqArray)){
    require(BiocManager)
    BiocManager::install(pkgs = "SeqArray")
  }

  if(!file.exists(paste0(wd,"/","P_in_vcf.GDS"))){
    snpgdsVCF2GDS(vcf.fn = P_in_vcf, out.fn = paste0(wd,"/","P_in_vcf.GDS"), snpfirstdim = TRUE)
  }

  if(!file.exists(paste0(wd, "/","IP_in_vcf.GDS"))){
    snpgdsVCF2GDS(vcf.fn = IP_in_vcf, out.fn = paste0(wd,"/","IP_in_vcf.GDS"), snpfirstdim = TRUE)
  }

  Pgdsinp <- snpgdsOpen(filename = paste0(wd,"/","P_in_vcf.GDS"))
  IPgdsinp <- snpgdsOpen(filename = paste0(wd,"/","IP_in_vcf.GDS"))

  P_sample_id <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "sample.id"))
  P_snp_id <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "snp.id"))
  P_snp_pos <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "snp.position"))
  P_snp_chr <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "snp.chromosome"))
  P_snp_allele <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "snp.allele"))
  P_genotype <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "genotype"))

  IP_sample_id <- read.gdsn(node = index.gdsn(node = IPgdsinp, index = "sample.id"))
  IP_snp_id <- read.gdsn(node = index.gdsn(node = IPgdsinp, index = "snp.id"))
  IP_snp_pos <- read.gdsn(node = index.gdsn(node = IPgdsinp, index = "snp.position"))
  IP_snp_chr <- read.gdsn(node = index.gdsn(node = IPgdsinp, index = "snp.chromosome"))
  IP_snp_allele <- read.gdsn(node = index.gdsn(node = IPgdsinp, index = "snp.allele"))
  IP_genotype <- read.gdsn(node = index.gdsn(node = IPgdsinp, index = "genotype"))
  ##########

  # Extracting reference haplotypes

  refgen <- P_genotype
  targetgen <- IP_genotype
  refgen[refgen == 3] <- NA
  targetgen[targetgen == 3] <- NA
  newtargetgen <- targetgen

  # Creating the sliding window

  win <- win
  step <- win
  pos <- seq(from = 1, to = nrow(x = refgen))
  sp <- seq(from = 1, to = c(nrow(x = refgen)-step), by = step)
  ep <- sp + win -1

  if(tail(x = ep, n = 1) != nrow(x = refgen)){
    ep[length(x = ep) + 1] <- nrow(x = refgen)
    sp[length(x = sp) + 1] <- tail(x = ep, n = 1) - win
  }


  if(length(x = sp) != length(x = ep)){stop("Sliding-window with different number of coordinates")}

  for(i in seq(from = 1, to = length(x = sp))){

    # Setting the window

    rf1 <- refgen[sp[i]:ep[i],1]
    rf2 <- refgen[sp[i]:ep[i],2]
    lc = i

    for(j in seq(from = 1, to = ncol(x = targetgen))){

      # Calculating the distances

      tg <- targetgen[sp[i]:ep[i],j]
      wna <- is.na(x = tg)
      drf1 <- dist(rbind(tg, rf1))
      drf2 <- dist(rbind(tg, rf2))

      # Windows-based imputation

      ## Empty window

      if(is.na(x = drf1 == drf2)){

        if(lc == 1){

          lcdrf1 <- NA
          lcdrf2 <- NA
          nlc <- lc

          while(is.na(lcdrf1 == lcdrf2)){
            nlc <- nlc+1
            if(nlc > length(x = sp)){stop("Error in the loop control #1")}
            lcrf1 <- refgen[sp[nlc]:ep[nlc],1]
            lcrf2 <- refgen[sp[nlc]:ep[nlc],2]
            lctg <- targetgen[sp[nlc]:ep[nlc],j]
            lcwna <- is.na(x = lctg)
            lcdrf1 <- dist(rbind(lctg, lcrf1))
            lcdrf2 <- dist(rbind(lctg, lcrf2))
          }

          if(lcdrf1 < lcdrf2){
            lcntg <- lctg
            lcntg[lcwna] <- lcrf1[lcwna]
          }

          if(lcdrf1 > lcdrf2){
            lcntg <- lctg
            lcntg[lcwna] <- lcrf2[lcwna]
          }

          if(lcdrf1 == lcdrf2){
            lcntg <- lctg
            lcntg[lcwna] <- c(lcrf1[lcwna]+lcrf2[lcwna])/2
          }

          newtargetgen[sp[i]:ep[i],j] <- lcntg

        }

        #

        if(lc > 1){

          lcdrf1 <- NA
          lcdrf2 <- NA
          nlc <- lc

          while(is.na(lcdrf1 == lcdrf2)){
            nlc <- nlc-1
            if(nlc < 1){stop("Error in the loop control #2")}
            lcrf1 <- refgen[sp[nlc]:ep[nlc],1]
            lcrf2 <- refgen[sp[nlc]:ep[nlc],2]
            lctg <- targetgen[sp[nlc]:ep[nlc],j]
            lcwna <- is.na(x = lctg)
            lcdrf1 <- dist(rbind(lctg, lcrf1))
            lcdrf2 <- dist(rbind(lctg, lcrf2))
          }

          if(lcdrf1 < lcdrf2){
            lcntg <- lctg
            lcntg[lcwna] <- lcrf1[lcwna]
          }

          if(lcdrf1 > lcdrf2){
            lcntg <- lctg
            lcntg[lcwna] <- lcrf2[lcwna]
          }

          if(lcdrf1 == lcdrf2){
            lcntg <- lctg
            lcntg[lcwna] <- c(lcrf1[lcwna]+lcrf2[lcwna])/2
          }

          newtargetgen[sp[i]:ep[i],j] <- lcntg

        }

        #


      } else {

        if(drf1 < drf2){
          ntg <- tg
          ntg[wna] <- rf1[wna]
        }

        if(drf1 > drf2){
          ntg <- tg
          ntg[wna] <- rf2[wna]
        }

        if(drf1 == drf2){
          ntg <- tg
          ntg[wna] <- c(rf1[wna]+rf2[wna])/2
        }
      }

      # Imputing the data

      newtargetgen[sp[i]:ep[i],j] <- ntg
    }
  }

  newtargetgen[is.na(x = newtargetgen)] <- 3


  ##########
  snpgdsCreateGeno(gds.fn = paste0(wd,"/",IP_new_out_vcf,".GDS"), genmat = newtargetgen, sample.id = IP_sample_id, snp.id = IP_snp_id, snp.chromosome = IP_snp_chr, snp.position = IP_snp_pos, snp.allele = IP_snp_allele, snpfirstdim = TRUE)
  seqSNP2GDS(gds.fn = paste0(wd,"/",IP_new_out_vcf,".GDS"), out.fn = paste0(wd,"/",IP_new_out_vcf,".SEQ"))
  seqGDS2VCF(gdsfile = paste0(wd,"/",IP_new_out_vcf,".SEQ"), vcf.fn = paste0(wd,"/",IP_new_out_vcf,".VCF"))
}
