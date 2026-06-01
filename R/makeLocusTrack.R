#' @export
#' @importFrom magrittr %>%
makeLocusTrack <- function(DictList, gtf, chr, start, end, annotation_bed=NULL, rel=NULL, cols=NULL,
                           meta=NULL, facet=T, include.all.annotations=F, plot_title=NULL){


  if(! is.null(annotation_bed)){

    if(any(gsub("chr", "", annotation_bed[[1]], ignore.case=T) != gsub("chr", "", chr, ignore.case = T))){

      message("Identified regions from other chromosomes in annotation file!")
      annotation_bed <- annotation_bed[gsub("chr", "", annotation_bed[[1]], ignore.case=T) %in% c(chr,  gsub("chr", "", chr, ignore.case = T)),]

      if(nrow(annotation_bed)==0){
        warning("No annotations found on input chromosome! Plotting without annotations...")
        annotation_bed <- NULL
      } else {
        message("Removed annotations from other chromosomes.")
      }
    }

    if(! is.null(annotation_bed) & include.all.annotations==T){
      start2 <- min(c(min(annotation_bed[[2]]), start))
      end2 <- max(c(max(annotation_bed[[3]]), end))
      if(! identical(c(min(c(min(annotation_bed[[2]]), start)), min(c(min(annotation_bed[[3]]), end))), c(start, end))){
        message("Expanding input locus boundaries to match annotation file...")
        start2 -> start
        end2 -> end
      }
    }
  }

  message("Extracting genomic information in locus...")
  tmp_gtf <- subset(gtf, V3 %in% c("transcript", "exon") &
                      V1 %in% c(chr,  gsub("chr", "", chr, ignore.case = T)) &
                      V5 > start &
                      V4 < end)

  # Determine whether or not to plot annotations based on input:
  if(nrow(tmp_gtf)==0){
    warning(paste0("No features present at: ", chr, ":", start, "-", end,  ". Plotting without genomic annotations..."))
    show.annotation <- FALSE
  } else {
    show.annotation <- TRUE
    # realign truncate annotations that overlap the locus boundary:
    tmp_gtf$V4 <- ifelse(tmp_gtf$V4 < start, start+1, tmp_gtf$V4)
    tmp_gtf$V5 <- ifelse(tmp_gtf$V5 > end, end-1, tmp_gtf$V5)

    tmp_gtf$V9 <- gsub("\\;.*", "", gsub(".*transcript_name", "", tmp_gtf$V9))
    tmp_gtf$ypos <- as.numeric(factor(tmp_gtf$V9))
  }

  if(is.null(plot_title)){
    plot_title <- paste0(chr, ":", start, "-", end)
  }

  # extract bin data across locus:
  input <- MakeGeneData(DictList, c=chr, s=start, e=end)

  # add meta data annotations:
  if(!is.null(meta)){
    message("Adding supplied meta data...")
    colnames(meta)[1:2] <- c("Set", "Fill")
    input <- merge(input, meta, by="Set", all.x=T)
  }

  # plot genomic annotation if applicable:
  if(show.annotation == T){
    arrow_ticks <- subset(tmp_gtf, V3 == "transcript") %>%
      dplyr::rowwise() %>%
      dplyr::reframe(x= seq(V4, V5, length.out = 7)[2:6], y= ypos, strand = V7)

    ggplot2::ggplot(tmp_gtf)+
      ggplot2::scale_x_continuous(paste("\nChromosome", tmp_gtf$V1[1], "(bp)"),
                                  expand=c(0,0),
                                  limits=c(start, end))+
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
  }

  # make peak tracks:
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
    ggplot2::ggtitle(plot_title)+
    ggplot2::scale_x_continuous(expand=c(0,0))+
    ggplot2::scale_y_continuous("Mean BPM", expand=c(0,0)) -> trackgrob

  if(! is.null(annotation_bed)){
    # add annotations to the track:
    tmpdf <- annotation_bed[c(2:3)]
    tmpdf <- stats::setNames(tmpdf, c("xs", "xe"))
    # remove peak annotations outside of range
    tmpdf <- subset(tmpdf, !xs > end & !xe < start)
    tmpdf$xs <- ifelse(tmpdf$xs < start, start, tmpdf$xs)
    tmpdf$xe <- ifelse(tmpdf$xe > end, end, tmpdf$xe)
    #
    tmpdf$ys <- 0
    tmpdf$ye <- Inf
    trackgrob+
      ggplot2::annotate(geom="rect", xmin=tmpdf$xs, xmax=tmpdf$xe, ymin=tmpdf$ys, ymax=tmpdf$ye,
               color=NA, fill="red", alpha=0.2) -> trackgrob
  }

  if(show.annotation == T){
    if(is.null(rel)){
      rel <- length(DictList)/ (log10(  ((nrow(tmp_gtf))/4)+1    ))
      message(paste0("Estimating plotting window... using rel=", formatC(rel, digits=3), "; you can increase rel to manually increase the relative size of the peak track window."))
    }
    outgrob <- ggpubr::ggarrange(trackgrob, annotgrob, heights = c(rel,1), align="v", ncol=1, nrow=2)
    return(outgrob)
  } else {
    return(trackgrob)
  }
}

