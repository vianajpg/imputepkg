% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parentsim.R
\name{parentsim}
\alias{parentsim}
\title{Simulate Parents Genotypes}
\usage{
parentsim(
  nsnps = nsnps,
  chrname = chrname,
  minpos = 1,
  maxpos = 5e+05,
  Hom1 = 0.25,
  Hom2 = 0.25,
  seed = 123456
)
}
\arguments{
\item{nsnps}{Number of requested SNPs}

\item{chrname}{Name for the simulated chromosome}

\item{minpos}{Initial position for the simulated chromosome}

\item{maxpos}{Final position for the simulated chromosome}

\item{Hom1}{Proportion of REF homozygotes}

\item{Hom2}{Proportion of ALT homozygotes}

\item{seed}{Seed for the sampling steps}
}
\value{
A data frame with the simulated genotypes and genomic information
}
\description{
This function simulate a genotype data for one chromosome in two parents P1 and P2
}
