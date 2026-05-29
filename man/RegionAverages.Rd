\name{RegionAverages}
\alias{RegionAverages}
\title{Plot the mean value across a set of regions.}
\usage{
RegionAverages(DictList, regions, cols=NULL)
}
\description{
Makes a boxplot of the mean value across a user specified set of regions for all samples in a DictList.
}
\arguments{
\item{DictList}{list of chromDicts to use. One dict becomes one track.}
\item{regions}{Bed-formatted data.frame containing loci to pull data from.}
\item{cols}{Optional vector of colours to use for boxes.}}
\examples{
# RegionAverages(dictList, data.frame(c=c("chr10", "chr10", "chr2", "chr1", "chr10"),
#                                     s=c(42200000,42203000,42204000, 42204000, 42290000),
#                                     e=c(42202000,42203500,42205000, 42205000, 42295000)))
}
