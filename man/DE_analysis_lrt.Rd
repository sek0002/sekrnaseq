% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEscript_lrt.R
\name{DE_analysis_lrt}
\alias{DE_analysis_lrt}
\title{Limma/voom and EdgeR differential analysis}
\usage{
DE_analysis_lrt(
  countdata,
  sampleinfo,
  thresh = 0.5,
  keep,
  group_param,
  exp_name,
  design,
  DEgroup,
  cont.matrix,
  species,
  genekeytype,
  archived = "FALSE",
  plots = F
)
}
\arguments{
\item{countdata}{Raw counts table with entrez or ensembl IDs (output from merged bam files)}

\item{sampleinfo}{Sampleinfo file (see template)}

\item{thresh}{Minimum CPM threshold to exclude}

\item{keep}{Minimum number of samples below thresh before exclusion (value usually set to number of replciates per group)}

\item{group_param}{select group_param column to compare from sampleinfo}

\item{exp_name}{Name to call experiment}

\item{design}{Design matrix (use model.matrix to generate, see template)}

\item{DEgroup}{Specify comparison of interest for DE}

\item{cont.matrix}{Contrast matrix (use makeContrast to generate, see template)}

\item{species}{"human" or "mouse"}

\item{genekeytype}{"entrez" or "ensembl"}

\item{archived}{Should be set to FALSE, unless using archived Human Ensembl database}

\item{plots}{Set to TRUE to generate TIFF and PDF volcano and MAplots in Volcano/exp_name}

\item{log}{Set to true to output normalised LogCPM}

\item{qlftest}{Default=TRUE, set to FALSE for edgeR LRT test}

\item{ntop}{default top 100, for MDSplot function only}
}
\value{
Table of normalised and filtered DGElist objects for limma/voom and EdgeR
}
\description{
DEanalysis with Limma/voom and EdgeR, with outputs.
}
