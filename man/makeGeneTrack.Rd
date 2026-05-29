\name{makeGeneTrack}
\alias{makeGeneTrack}
\title{Make a peak track at a specified gene.}
\usage{
makeGeneTrack(DictList, gtf, gene, annotation_bed=NULL, rel=NULL, cols=NULL,
              meta=NULL, facet=T, include.all.annotations=F, buffer=1000)
}
\description{
Creates peak tracks as grobs in ggplot using chromDicts at a specific gene.
}
\arguments{
\item{DictList}{list of chromDicts to use. One dict becomes one track.}
\item{gtf}{GTF containing genomic annotations. Designed for gtfs from ensembl.}
\item{gene}{Name of gene to plot as character.}
\item{annotation_bed}{Optional: bed-formatted data.frame containing loci to highlight.}
\item{rel}{Optional: numeric, the relative window ratio of the track windows to the annotation window. Larger value means the peak tracks take up a larger portion of vertical space. Will estimate based on the number of tracks/transcripts if not set.}
\item{cols}{Vector of colours to use for peak tracks.}
\item{meta}{Two column data.frame containing track names (names(DictList)) and Condition; will colour tracks by the 'Condition' column.}
\item{facet}{Boolean, if TRUE, adds sample or condition names to peak tracks. Default=TRUE.}
\item{include.all.annotations}{Boolean, if TRUE, ALL annotations in annotation_bed will be plotted, even if they are outside of the gene region. May cause errors if some sites are on a different chromosome.}
\item{buffer}{Numeric, number of bp to add to each end of the gene region for plotting. Default=1000.}}
\examples{
# Basic peak track at Cox5b using 2 chromDicts to plot 2 tracks:
# makeTrack(DictList=list(chromDict1, chromDict2),
#           gtf=gtf,
#           gene="Cox5b",
#           rel=2,
#           meta = data.frame(v1=c("chromDict1", "chromDict2"),
#                             v2=c("treatment", "control")),
#           cols=c("grey", "skyblue"))
}
