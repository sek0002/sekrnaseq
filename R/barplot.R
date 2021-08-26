#' Barplot of scripts using fgsea algorithm
#'
#'
#' @param countdata Raw counts table with entrez or ensembl IDs (output from merged bam files)

#' @return Barplot from GSEA scrip
#' @import edgeR
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
#' @export


barplot_gsea <- function(height, x="NES", color='p.adjust', showCategory=8, font.size=12, title="", xlab=NULL, ylab=NULL,orderBy = "x", decreasing = TRUE, ...) {
  ## use *height* to satisy barplot generic definition
  ## actually here is an enrichResult object.
  object <- height

  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  }
  else if (x == "count" || x == "Count") {
    x <- "Count"
  }

  df <- fortify(object, showCategory=showCategory, by=x, ...)
  if (orderBy != "x" && !orderBy %in% colnames(df)) {
    message("wrong orderBy parameter; set to default `orderBy = \"x\"`")
    orderBy <- "x"
  }
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text = x)))
  }
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels = rev(unique(df$Description[idx])))


  if(colorBy %in% colnames(df)) {
    p <- ggplot(df, aes_string(x = "Description", y = x, fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE))
  } else {
    p <- ggplot(df, aes_string(x = "Description", y = x, fill = "Description")) +
      theme_dose(font.size) +
      theme(legend.position="none")
  }
  p + geom_bar(stat = "identity") + coord_flip() +
    ggtitle(title) + xlab(xlab) + ylab(ylab)
}
