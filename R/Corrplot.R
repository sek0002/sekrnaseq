#' Interactive correlation plot
#'
#' using two data tables, generate interactive correlation plot
#'
#' @param frame1 1st data frame
#' @param frame2 2nd data frame
#' @param com_col Common column to match rows (ie. gene symbols)
#' @param FCcol label common fold change columns
#' @param FDRcol label common FDR/adj P value columns
#' @param namex name for x axis
#' @param namey name for y axis
#'
#' @return plotly interactive correlation plot of input dataframes
#'
#' @import ggplot2
#' @import plotly
#' @import ggrepel
#' @import viridis
#' @import Glimma
#' @import ggrastr
#' @export


##interactive correlation plots between two tables
corrplot <- function(frame1, frame2, com_col, FCcol, FDRcol, namex=NULL, namey=NULL) {
  colnames(frame1) <- paste("x", colnames(frame1), sep = "_")
  colnames(frame1)[names(frame1) == paste0("x_",com_col)]<- paste0(com_col)
  colnames(frame2) <- paste("y", colnames(frame2), sep = "_")
  colnames(frame2)[names(frame2) == paste0("y_",com_col)]<- paste0(com_col)
  frame1 <- na.omit(frame1)
  frame2 <- na.omit(frame2)
  filemerge <- Reduce(function(dtf1,dtf2)
    merge(dtf1,dtf2, by =com_col), list(frame1,frame2))
  corrmatrix <- filemerge[c(com_col, paste0("x_",FCcol,sep=""), paste0("y_",FCcol,sep=""), paste0("x_",FDRcol, sep=""), paste0("y_",FDRcol,sep=""))]

  corrmatrix["group"] <- "notsig"
  corrmatrix[which(corrmatrix[paste0("x_",FDRcol,sep="")] <=0.05), "group"] <- "FDRxsig"
  corrmatrix[which(corrmatrix[paste0("y_",FDRcol,sep="")] <=0.05), "group"] <- "FDRysig"
  corrmatrix<- within(corrmatrix, group[corrmatrix[paste0("x_",FDRcol, sep="")] <=0.05 & corrmatrix[paste0("y_",FDRcol, sep="")] <=0.05] <- "Overlap")

  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )

  x <- list(
    title = paste0(namex,sep=""),
    titlefont = f
    #range =c(-15,15)
  )

  y <- list(
    title = paste0(namey,sep=""),
    titlefont = f
    #range =c(-15,15)

  )

  p <- plot_ly(data = corrmatrix, x = corrmatrix[,2], y = corrmatrix[,3], text = corrmatrix$SYMBOL, color = corrmatrix$group) %>%
    layout(title ="Correlation") %>%
    layout(xaxis=x , yaxis=y)
  p
  # write.table(corrmatrix, file = paste0("DE/",namex,namey,".txt"), sep='\t',row.names = FALSE)

}
