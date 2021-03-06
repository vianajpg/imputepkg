% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputeVCF_haploslider.R
\name{imputeVCF_haploslider}
\alias{imputeVCF_haploslider}
\title{Impute Inbred VCF using sliding-window}
\usage{
imputeVCF_haploslider(
  P_in_vcf = "",
  IP_in_vcf = "",
  IP_new_out_vcf = "",
  win = 5,
  step = 2
)
}
\arguments{
\item{P_in_vcf}{Parents original VCF}

\item{IP_in_vcf}{Inbred Populations original VCF}

\item{IP_new_out_vcf}{Name for an output - Imputed Inbred VCF}

\item{win}{Size of the sliding window}

\item{step}{Amount of SNPs jumped after each step}
}
\value{
Imputed inbred populations VCF
}
\description{
This function impute VCF data from the inbred populations based on the VCF data from the parents
}
