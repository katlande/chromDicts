#' @noRd
MakeGeneData <- function(chromDictList, c, s, e){

  for(i in 1:length(chromDictList)){
    chromDict <- chromDictList[[i]]
    tmp <- chromDict[[c]]
    tmp <- stats::na.omit(tmp[ .( c(s:e) ) ])
    tmp$Set <- names(chromDictList)[i]
    if(exists("OUTPUTFILE")){
      OUTPUTFILE <- rbind(OUTPUTFILE, tmp)
    } else {
      OUTPUTFILE <- tmp
    }
  }
  return(OUTPUTFILE)
}

#' @export
#' @importFrom magrittr %>%
makeGeneTrack <- function(DictList, gtf, gene, annotation_bed=NULL, rel=4, cols=NULL,
                      has.chr=T, meta=NULL, facet=T, include.all.annotations=F){

  message("Extracting gene information")
  tmp_gtf <- subset(gtf, V3 %in% c("gene", "transcript", "exon") & grepl(paste0(gene, ";"), V9))
  if(nrow(tmp_gtf)==0){
    warning(paste("No gene named", gene, "present in gtf! Maybe it has another name?"))
  }
  tmp_gtf$V9 <- gsub("\\;.*", "", gsub(".*transcript_id", "", tmp_gtf$V9))
  tmp_gtf$ypos <- as.numeric(factor(tmp_gtf$V9))

  message("Pulling values across gene body...")
  c <-ifelse(has.chr,  paste0("chr",tmp_gtf$V1[tmp_gtf$V3 == "gene"][1]),  tmp_gtf$V1[tmp_gtf$V3 == "gene"][1])
  s <- tmp_gtf$V4[tmp_gtf$V3 == "gene"][1]-1000
  e <- tmp_gtf$V5[tmp_gtf$V3 == "gene"][1]+1000

  if(! is.null(annotation_bed) & include.all.annotations==T){
    s <- min(c(min(annotation_bed[[2]]), s))
    e <- max(c(max(annotation_bed[[3]]), e))
  }

  tmp_gtf <- subset(tmp_gtf, V3 %in% c("transcript", "exon"))
  input <- MakeGeneData(DictList, c=c, s=s, e=e)

  if(!is.null(meta)){
    message("Adding supplied meta data...")
    colnames(meta)[1:2] <- c("Set", "Fill")
    input <- merge(input, meta, by="Set", all.x=T)
  }

  # gene annotation
  message("Creating annotation track...")
  arrow_ticks <- subset(tmp_gtf, V3 == "transcript") %>%
    dplyr::rowwise() %>%
    dplyr::reframe(x= seq(V4, V5, length.out = 7)[2:6], y= ypos, strand = V7)

  ggplot2::ggplot(tmp_gtf)+
    ggplot2::scale_x_continuous(paste("\nChromosome", tmp_gtf$V1[1], "(bp)"),
                                expand=c(0,0),
                                limits=c(min(input$pos), max(input$pos)))+
    ggplot2::geom_segment(data=subset(tmp_gtf, V3 == "transcript"),
                          mapping= ggplot2::aes(x=V4, xend=V5, y=ypos, yend=ypos))+
    ggplot2::geom_rect(data=subset(tmp_gtf, V3 == "exon"),
                       mapping= ggplot2::aes(xmin=V4, xmax=V5, ymax=ypos+0.25, ymin=ypos-0.3),  fill="black")+
    ggplot2::theme(panel.grid.major.y =  ggplot2::element_blank(),
          panel.grid.major.x =  ggplot2::element_line(color="grey", linewidth = 0.3),
          axis.text.y =  ggplot2::element_text(colour = "black", size=6),
          axis.text.x =  ggplot2::element_text(colour = "black"),
          panel.background =  ggplot2::element_rect(fill="white", color="black"),
          axis.title.y = ggplot2::element_blank(),
          plot.title= ggplot2::element_blank(),
          plot.margin =  ggplot2::margin(0, 0, 0.25, 0, "cm"))+
    ggplot2::scale_y_continuous(labels=levels(factor(tmp_gtf$V9)), expand=c(0.1,0.1),
                                breaks=unique(tmp_gtf$ypos)[order(unique(tmp_gtf$ypos))])+
    ggplot2::geom_segment(data = arrow_ticks,
                          mapping =  ggplot2::aes(x = x, y = y, yend = y,
                               xend = ifelse(strand == "+", x + diff(range(input$pos)) * 0.01,
                                             x - diff(range(input$pos)) * 0.01)),
                 arrow =  ggplot2::arrow(length = unit(3, "pt"), type = "closed"),
                 linewidth = 0.2, color = "black")-> annotgrob

  # make peaks
  message("Creating peak track(s)...")
  if(ncol(input) == 4 & facet==T){
    ggplot2::ggplot(input,  ggplot2::aes(x=pos))+
      ggplot2::geom_ribbon(ymin = 0,  ggplot2::aes(ymax=value))+
      ggplot2::facet_wrap(~Set, nrow=length(unique(input$Set))) -> g

  } else if(ncol(input) > 4 & facet==T){

    if(is.null(cols)){
      cols <- grDevices::colorRampPalette(c("skyblue", "darkgrey"))(length(unique(input$Fill)))
    }

    ggplot2::ggplot(input, ggplot2::aes(x=pos))+
      ggplot2::geom_ribbon(ymin = 0, ggplot2::aes(ymax=value, fill=Fill))+
      ggplot2::facet_wrap(~Set, nrow=length(unique(input$Set)))+
      ggplot2::scale_fill_manual(values=cols) -> g

  } else {
    ggplot2::ggplot(input, ggplot2::aes(x=pos))+
      ggplot2::geom_ribbon(ymin = 0, ggplot2::aes(ymax=value))-> g
  }

  g+
    PCBS::Ol_Reliable()+
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, face="bold", size=14),
                   axis.title.x = ggplot2::element_blank(), legend.position = "none")+
    ggplot2::ggtitle(gene)+
    ggplot2::scale_x_continuous(expand=c(0,0))+
    ggplot2::scale_y_continuous("Mean CPM/bin", expand=c(0,0)) -> trackgrob

  if(! is.null(annotation_bed)){
    # add annotations to the track:
    tmpdf <- annotation_bed[c(2:3)]
    tmpdf <- stats::setNames(tmpdf, c("xs", "xe"))
    # remove peak annotations outside of range
    tmpdf <- subset(tmpdf, !xs > e & !xe < s)
    tmpdf$xs <- ifelse(tmpdf$xs < s, s, tmpdf$xs)
    tmpdf$xe <- ifelse(tmpdf$xe > e, e, tmpdf$xe)
    #
    tmpdf$ys <- 0
    tmpdf$ye <- Inf
    trackgrob+
      ggplot2::annotate(geom="rect", xmin=tmpdf$xs, xmax=tmpdf$xe, ymin=tmpdf$ys,
                        ymax=tmpdf$ye, color=NA, fill="red", alpha=0.2) -> trackgrob
  }

  outgrob <- ggpubr::ggarrange(trackgrob, annotgrob, heights = c(rel,1), align="v", ncol=1, nrow=2)
  return(outgrob)
}
