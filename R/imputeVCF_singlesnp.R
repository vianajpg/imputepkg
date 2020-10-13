#' Impute Offspring VCF per SNP
#'
#' This function impute VCF data from the offsprings based on the VCF data from the parents
#'
#' @param P_in_vcf Parents original VCF
#' @param OS_in_vcf Offspring original VCF
#' @param OS_new_out_vcf Name for an output - Imputed Offspring VCF
#' @return Imputed Offspring VCF
#' @export
imputeVCF_singlesnp <- function(P_in_vcf = "", OS_in_vcf = "", OS_new_out_vcf = ""){
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
  snpgdsVCF2GDS(vcf.fn = P_in_vcf, out.fn = paste0(wd,"P_in_vcf.GDS"), snpfirstdim = TRUE)
  snpgdsVCF2GDS(vcf.fn = OS_in_vcf, out.fn = paste0(wd,"OS_in_vcf.GDS"), snpfirstdim = TRUE)

  Pgdsinp <- snpgdsOpen(filename = paste0(wd,"P_in_vcf.GDS"))
  OSgdsinp <- snpgdsOpen(filename = paste0(wd,"OS_in_vcf.GDS"))

  P_sample_id <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "sample.id"))
  P_snp_id <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "snp.id"))
  P_snp_pos <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "snp.position"))
  P_snp_chr <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "snp.chromosome"))
  P_snp_allele <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "snp.allele"))
  P_genotype <- read.gdsn(node = index.gdsn(node = Pgdsinp, index = "genotype"))

  OS_sample_id <- read.gdsn(node = index.gdsn(node = OSgdsinp, index = "sample.id"))
  OS_snp_id <- read.gdsn(node = index.gdsn(node = OSgdsinp, index = "snp.id"))
  OS_snp_pos <- read.gdsn(node = index.gdsn(node = OSgdsinp, index = "snp.position"))
  OS_snp_chr <- read.gdsn(node = index.gdsn(node = OSgdsinp, index = "snp.chromosome"))
  OS_snp_allele <- read.gdsn(node = index.gdsn(node = OSgdsinp, index = "snp.allele"))
  OS_genotype <- read.gdsn(node = index.gdsn(node = OSgdsinp, index = "genotype"))

  newOS <- data.frame(matrix(data = NA, nrow = 0, ncol = length(OS_sample_id) + 4))
  for(i in unique(x = P_snp_chr)){
    P_subchr <- P_snp_chr[P_snp_chr == i]
    P_subpos <- P_snp_pos[P_snp_chr == i]
    P_subid <- P_snp_id[P_snp_chr == i]
    P_suball <- P_snp_allele[P_snp_chr == i]
    P_subgen <- P_genotype[P_snp_chr == i,]
    P_subdata <- data.frame(P_subchr, P_subpos, P_subid, P_suball, P_subgen)
    rm(c("P_subchr", "P_subpos", "P_subid", "P_suball", "P_subgen"))
    OS_subchr <- OS_snp_chr[OS_snp_chr == i]
    OS_subpos <- OS_snp_pos[OS_snp_chr == i]
    OS_subid <- OS_snp_id[OS_snp_chr == i]
    OS_suball <- OS_snp_allele[OS_snp_chr == i]
    OS_subgen <- OS_genotype[OS_snp_chr == i,]
    OS_subdata <- data.frame(OS_subchr, OS_subpos, OS_subid, OS_suball, OS_subgen)
    rm(c("OS_subchr", "OS_subpos", "OS_subid", "OS_suball", "OS_subgen"))
    rbind(newOS,singlesnpinp(OS = OS_subdata, P = P_subdata))
  }

  snpgdsCreateGeno(gds.fn = paste0(wd,OS_new_out_vcf,".GDS"), genmat = newOS[,5:ncol(newOS)], sample.id = OS_sample_id, snp.id = newOS[,3], snp.chromosome = newOS[,1], snp.position = newOS[,2], snp.allele = newOS[,4], snpfirstdim = TRUE)
  seqSNP2GDS(gds.fn = paste0(wd,OS_new_out_vcf,".GDS"), out.fn = paste0(wd,OS_new_out_vcf,".SEQ"))
  seqGDS2VCF(gdsfile = paste0(wd,OS_new_out_vcf,".SEQ"), vcf.fn = paste0(wd,OS_new_out_vcf,".VCF"))
}
