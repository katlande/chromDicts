\name{chromDictAny}
\alias{chromDictAny}
\title{Turn any value into a chromDict.}
\usage{
chromDictAny(mat, IDs=NULL, multiple.samples=T, remove.extra=T, n.autosome=22)
}
\description{
Takes a data.frame containing genome-wide values and creates a chromDict. Can either take a single sample or take an average across samples.
}
\arguments{
\item{mat}{Input matrix, ideally prepped by NormalizeBins().}
\item{IDs}{Name(s) of sample(s) to use, as character or vector or characters. If NULL, takes the average of all samples.}
\item{multiple.samples}{Boolean, whether or not multiple samples are being averaged. Default=TRUE.}
\item{remove.extra}{Removes abnormal chromosome names (anything except autosomes & x/y). Suggested to turn this off if using an genome with different naming conventions from mouse and human. Default=TRUE.}
\item{n.autosome}{If remove.extra==TRUE, specifies how many autosomes the chromosome contains. For human: 22, for mouse: 19. Default=22.}
}
\examples{
# average all trt samples into a single value:
# chromDict <- chromDictAny(mat, IDs=c("trt1", "trt2", "trt3"))

# one chromDict per sample:
# chromDict_trt1 <- chromDictAny(mat, IDs="trt1", multiple.samples=F)
# chromDict_trt2 <- chromDictAny(mat, IDs="trt2", multiple.samples=F)
# chromDict_trt3 <- chromDictAny(mat, IDs="trt3", multiple.samples=F)

}
