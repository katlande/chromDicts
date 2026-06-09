#' @export
#' @importFrom data.table data.table
MakeGeneData <- function(chromDictList, c, s, e, force=TRUE){

  for(i in 1:length(chromDictList)){
    chromDict <- chromDictList[[i]]
    tmp <- chromDict[[c]]
    tmp <- stats::na.omit(tmp[ .( c(s:e) ) ])
    tmp$Set <- names(chromDictList)[i]
    if(exists("OUTPUTFILE")){
      if(force==TRUE){
        OUTPUTFILE <- rbind(OUTPUTFILE, setNames(tmp, colnames(OUTPUTFILE)))
      } else {
        OUTPUTFILE <- rbind(OUTPUTFILE, tmp)
      }
    } else {
      OUTPUTFILE <- tmp
    }
  }
  return(OUTPUTFILE)
}
