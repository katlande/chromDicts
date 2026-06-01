\name{MakeGeneData}
\alias{MakeGeneData}
\title{Pull data from a list of chromDicts}
\usage{
MakeGeneData(chromDictList, c, s, e)
}
\description{
Pull data from a list of chromDicts at a specified locus.
}
\arguments{
\item{chromDictList}{Named list of chromDicts.}
\item{c}{Chromosome of locus.}
\item{s}{Integer, start position of locus.}
\item{e}{Integer, end position of locus.}
}
\examples{
# average all trt samples into a single value:
# df <- chromDictAny(chromDictList, "chr1", 1500000, 1600000)

}
