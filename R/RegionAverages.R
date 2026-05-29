#' @export
#' @importFrom magrittr %>%
RegionAverages <- function(DictList, regions, cols=NULL){
  apply(regions, 1, function(x){
    MakeGeneData(DictList, x[[1]], x[[2]], x[[3]]) %>%
      dplyr::group_by(Set) %>%
      dplyr::summarise(mean(value)) %>%
      as.data.frame()  %>%
      setNames(c("Set", "meanValue"))
  }) %>%
    do.call(rbind,.) -> g

  if(is.null(cols)){
    cols <- grDevices::colorRampPalette(c("#E7298A", "#D95F02", "#A6761D","#E6AB02", "#66A61E","#1B9E77","#666666","#7570B3"))(length(unique(g$Set)))
  }

  ggplot2::ggplot(g, ggplot2::aes(x=Set, fill=Set, y=meanValue))+
    ggplot2::geom_boxplot(outlier.alpha = 0)+
    ggplot2::geom_jitter(height=0, width=0.075, shape=21)+
    PCBS::Ol_Reliable()+
    PCBS::tilt()+
    ggplot2::scale_fill_manual(values=cols)+
    ggplot2::theme(legend.title = ggplot2::element_blank())+
    ggplot2::xlab("")+
    ggplot2::ylab("Mean Value per Region") -> gp

  return(gp)
}
