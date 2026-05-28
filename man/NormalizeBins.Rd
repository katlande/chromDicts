\name{NormalizeBins}
\alias{NormalizeBins}
\title{Normalize bin data from bam files}
\usage{
NormalizeBins(input_file, output_file=NULL, constant_binsize=T, sample_names=NULL)
}
\description{
Takes a bed file output from multibamsummary and converts it into a BPM-normalized matrix ready to be turned into a chromDict.
}
\arguments{
\item{input_file}{Path to the output bed file from multibamsummary.}
\item{output_file}{Path to the output save location; if NULL, will only write to working environment. Default=NULL.}
\item{constant_binsize}{Whether or not a uniform binsize was used. If TRUE, centers bins. If FALSE, uses the start position of the bin. Default=TRUE.}
\item{sample_names}{Optional character vector to rename the sample columns of the output data.frame.}
}
\examples{
# basic:
# BPM_file <- NormalizeBins("multibam_out.bed")

# rename samples:
# BPM_file <- NormalizeBins("multibam_out.bed", sample_names=c("trt1", "trt2", "ctl1", "ctl2"))

# save output to file:
# BPM_file <- NormalizeBins("multibam_out.bed", output_file="/path/to/BPMs.tsv")
}
