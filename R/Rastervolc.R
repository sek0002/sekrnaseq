#' Plotting RNA-seq data input from DE_analysis
#'
#' Using output file from DE_analysis to generate volcano, MA, correlation and MDS plots.
#'
#' @param data_in Output from DE_analysis
#' @param xaxis Input string for xaxis
#' @param yaxis Input string for yaxis
#' @param name_output Name to title plot
#' @param exp_name exp_name, also determine destination where file is saved
#'
#' @return ggplots of RNAseq data
#'
#' @import ggplot2
#' @import plotly
#' @import ggrepel
#' @import viridis
#' @import Glimma
#' @import ggrastr
#' @export

##Raster volc
rastervolc <- function(data_in,label_peaks,xaxis,yaxis,name_output,exp_name){

  # data_in$Threshold <- ifelse(data_in$FDR<=0.05|data_in$logFC <=1 & data_in$logFC >= -1 ,"nochange",
  #                               ifelse(data_in$logFC>1, "Open", "Close"))

  if(length(unique(data_in[,"Threshold"])) == 3){
    g1 <- ggplot(data = data_in, aes(x = data_in[,xaxis], y = -log10(data_in[,yaxis]))) +
      geom_point_rast(
        show.legend = FALSE,
        alpha = 1,
        size = 2,
        aes(col = data_in[,"Threshold"])
      ) +
      scale_colour_manual(values = c("Blue", "Grey", "Red")) +
      xlab("log2 fold change") +
      ylab("-log10(FDR)") +
      # ylim(0, 60) +
      # xlim(-3, 4) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      theme_linedraw() +
      theme(
        panel.grid = element_blank()
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # text = element_text(size = 20),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank()
      )


    sub<- data_in[data_in$SYMBOL %in% label_peaks,]
    out_vol<- g1+geom_text_repel(data=sub, aes(x= sub[,xaxis], y= -log10(sub[,yaxis]), label=SYMBOL),
                                 nudge_y = 0, nudge_x = 0)
    plot(out_vol)
    ggsave(file=paste0("Volcano/",exp_name,"/",name_output,".pdf", sep=""), plot=out_vol, units="in", width=7, height=6, useDingbats=FALSE)

  }else if(length(unique(data_in[,"Threshold"])) == 2){
    g1 <- ggplot(data = data_in, aes(x = data_in[,xaxis], y = -log10(data_in[,yaxis]))) +
      geom_point_rast(
        show.legend = FALSE,
        alpha = 1,
        size = 2,
        aes(col = data_in[,"Threshold"])
      ) +
      scale_colour_manual(values = c("Red", "Grey")) +
      xlab("log2 fold change") +
      ylab("-log10(FDR)") +
      # ylim(0, 60) +
      # xlim(-3, 4) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      theme_linedraw() +
      theme(
        panel.grid = element_blank()
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # text = element_text(size = 20),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank()
      )

    sub<- data_in[data_in$SYMBOL %in% label_peaks,]
    out_vol<- g1+geom_text_repel(data=sub, aes(x= sub[,xaxis], y= -log10(sub[,yaxis]), label=SYMBOL),
                                 nudge_y = 0, nudge_x = 0)
    plot(out_vol)
    ggsave(file=paste0("Volcano/",exp_name,"/",name_output,".pdf", sep=""), plot=out_vol, units="in", width=7, height=6, useDingbats=FALSE)
  }else{
    g1 <- ggplot(data = data_in, aes(x = data_in[,xaxis], y = -log10(data_in[,yaxis]))) +
      geom_point_rast(
        show.legend = FALSE,
        alpha = 1,
        size = 1,
        aes(col = data_in[,"Threshold"])
      ) +
      scale_colour_manual(values = c("Grey")) +
      xlab("log2 fold change") +
      ylab("-log10(FDR)") +
      # ylim(0, 60) +
      # xlim(-3, 4) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      theme_linedraw() +
      theme(
        panel.grid = element_blank()
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # text = element_text(size = 20),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank()
      )

    sub<- data_in[data_in$SYMBOL %in% label_peaks,]
    out_vol<- g1+geom_text_repel(data=sub, aes(x= sub[,xaxis], y= -log10(sub[,yaxis]), label=SYMBOL),
                                 nudge_y = 0, nudge_x = 0)
    plot(out_vol)
    ggsave(file=paste0("Volcano/",exp_name,"/",name_output,".pdf", sep=""), plot=out_vol, units="in", width=7, height=6, useDingbats=FALSE)

  }

}

