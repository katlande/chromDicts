\name{makeLocusTrack}
\alias{makeLocusTrack}
\title{Make a peak track across a specified region.}
\usage{
makeLocusTrack(DictList, gtf, chr, start, end, annotation_bed=NULL, rel=NULL, cols=NULL,
               meta=NULL, facet=T, include.all.annotations=F, plot_title=NULL)
}
\description{
Creates peak tracks as grobs in ggplot using chromDicts at a specific gene.
}
\arguments{
\item{DictList}{list of chromDicts to use. One dict becomes one track.}
\item{gtf}{GTF containing genomic annotations. Designed for gtfs from ensembl.}
\item{chr}{chromosome of plotting locus.}
\item{start}{start position of plotting locus.}
\item{end}{end position of plotting locus.}
\item{annotation_bed}{Optional: bed-formatted data.frame containing loci to highlight.}
\item{rel}{Optional: numeric, the relative window ratio of the track windows to the annotation window. Larger value means the peak tracks take up a larger portion of vertical space. Will estimate based on the number of tracks/transcripts if not set.}
\item{cols}{Vector of colours to use for peak tracks.}
\item{meta}{Two column data.frame containing track names names(DictList) and Condition; will colour tracks by the 'Condition' column.}
\item{facet}{Boolean, if TRUE, adds sample or condition names to peak tracks. Default=TRUE.}
\item{include.all.annotations}{Boolean, if TRUE, ALL annotations in annotation_bed will be plotted, even if they are outside of the gene region. May cause errors if some sites are on a different chromosome.}
\item{plot_title}{Optional: character, plot title. Will use "chr:start-end" if not set.}}
\examples{
# Basic peak track :
#makeLocusTrack(dictList, gtf, chr="chr1", start=135824455, end=135907423,
#               plot_title="test")

# peak track with annotations, meta data, and custom colours:
#makeLocusTrack(dictList, gtf, chr="chr1", start=135824455, end=135907423,
#               plot_title="test",
#               meta=trackMeta,
#               annotation_bed = data.frame(c=c("chr1", "chr1"),
#                                          s=c(135828555,135854455),
#                                          e=c(135829655,135855455)),
#               cols=c("navy", "salmon"))
}
