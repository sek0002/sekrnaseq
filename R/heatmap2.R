#' Heatmapscriptnew use heatmapscript if error
#'
#'
#' @param data CPM table with entrez or ensembl IDs (output of sekrnaseq::CPMlist function)

#' @return Heatmaps outputed to heatmap/exp_name folder
#'
#' @import ggplot2
#' @importFrom data.table setDT
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi mapIds
#' @import limma
#' @import ggplot2
#' @import plotly
#' @import ggrepel
#' @import viridis
#' @import Glimma
#' @import fgsea
#' @import superheat
#' @import tidyverse
#' @import heatmaply
#' @import dplyr
#' @import pheatmap
#' @export
heatmapscriptnew <- function(data=NULL,selcol=NULL, refno= 1, refgene= data, name = NULL, destination, clcol= TRUE){

  if (is.null(data)){
    stop('no reference gene list provided')}
  else if(is.null(selcol)){
    stop("select columns to plot")
  } else if(identical(refgene,data)){
    print("Showing all SYMBOLs")
  }
  if (refno == 1){
    print("Using first column for Symbols")
  }

  refgenelist <- refgene %>% na.omit()
  refgenelist <- refgenelist[,refno]
  subsetted_filein_avg <- data[which(data$SYMBOL %in% refgenelist),]
  subsetted_filein_unique_avg <- subsetted_filein_avg[!duplicated(subsetted_filein_avg$SYMBOL),]
  heat_map_matrix <- data.frame(subsetted_filein_unique_avg[,"SYMBOL"],subsetted_filein_unique_avg[,..selcol])
  print(nrow(heat_map_matrix))
  heat_map <- heat_map_matrix[,-1]
  rownames(heat_map) <- heat_map_matrix[,1]
  center_raw_mat_logged <- log2(heat_map+0.5)
  center_raw_mat_logged_matrix <- data.matrix(center_raw_mat_logged)

  set.seed(2019-03-21)
  heat<- pheatmap(center_raw_mat_logged, scale="row",  #double <<- spits asd out into global enviro
                  # treeheight_row = 0,   treeheight_col = 0,
                  file = paste(destination, "/",Sys.Date(),name, ".pdf", sep = ""),
                  fontsize_col = 7.5,
                  cluster_cols=clcol,
                  show_colnames = TRUE,
                  show_rownames = TRUE,
                  main = paste("" ,name, sep = ""),
                  cellwidth = 20,
                  # cellheight = 10,
                  color = colorRampPalette(c("#0099FF", "black", "#FEFF00"))(50), border_color = NA)
  return(heat)
}
