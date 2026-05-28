#' @export
NormalizeBins <- function(input_file, output_file=NULL, constant_binsize=T, sample_names=NULL){
  res <- utils::read.delim(input_file, comment.char = "#", header=F)

  if(constant_binsize == T){
    bs <- as.numeric(res[[2]][2]-res[[2]][1])
    message(paste("Found a bin size of:"), bs)
  }

  # normalize:
  res[4:ncol(res)] <-
    apply(res[4:ncol(res)], 2, function(x){
      return((as.numeric(x)/sum(as.numeric(x)))*1e06)
    })

  # Add column names to res:
  if(is.null(sample_names)){
    colnames(res) <- c("chr", "start", "end", colnames(res)[4:ncol(res)])
  } else {
    colnames(res) <- c("chr", "start", "end", sample_names)
  }

  # add a single position for each bin at the centerpoint:
  if(constant_binsize==T){
    res$pos <- res$start + as.integer(bs*0.5)
  } else {
    res$pos <- res$start
  }

  # remove start and end from res:
  res <- res[c(1,ncol(res),4:(ncol(res)-1))]

  if(is.null(output_file)){
    return(res)
  } else {
    message("saving output file...")
    utils::write.table(x = res, file = output_file, sep="\t", quote=F, col.names = T, row.names = F)
    return(res)
  }
}
