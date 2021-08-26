#' DOTplot for GSEA scripts using fgsea algorithm
#'
#'
#' @param countdata Raw counts table with entrez or ensembl IDs (output from merged bam files)

#' @return Dotplot of results from GSEAscript
#'
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
#'
dotplot_gsea<- function (object, x = "geneRatio", color = "p.adjust", showCategory = 20, size = NULL, split = NULL, font.size = 12, title = "", orderBy = "x", decreasing = TRUE) {
  colorBy <- match.arg(color, c("pvalue", "p.adjust",
                                "qvalue","NES"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  }
  else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  }
  else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else if (is.null(size)) {
    size <- "Count"
  } else if(size == "p.adjust") {
    size <- "p.adjust"
  }
  df <- fortify(object, showCategory = showCategory, split = split)
  if (orderBy != "x" && !orderBy %in% colnames(df)) {
    message("wrong orderBy parameter; set to default `orderBy = \"x\"`")
    orderBy <- "x"
  }
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text = x)))
  }
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels = rev(unique(df$Description[idx])))

  if (size=="p.adjust"){
    ggplot(df, aes_string(x = x, y = "Description", size = -log10(df$p.adjust),
                          color = colorBy)) + geom_point() + scale_color_continuous(low = "blue",
                                                                                    high = "red", name = color, guide = guide_colorbar(reverse = TRUE)) +
      ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
      scale_size(range = c(3, 8))
  } else if (colorBy=="NES" & x=="p.adjust") {
    ggplot(df, aes_string(x = -log10(df$x), y = "Description", size = size,
                          color = colorBy)) + geom_point() + scale_color_continuous(high = "red",
                                                                                    low = "blue", name = color) + xlab(paste0("-log10(",x,")"))+
      ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
      scale_size(range = c(3, 8))
  }else if (is.null(color)) {
    ggplot(df, aes_string(x = df$x, y = "Description",size = -log10(df$p.adjust)
    )) + geom_point() + xlab(paste0(x))+
      ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
      scale_size(range = c(3, 8))
  } else {
    ggplot(df, aes_string(x = x, y = "Description", size = size,
                          color = colorBy)) + geom_point() + scale_color_continuous(low = "red",
                                                                                    high = "blue", name = color, guide = guide_colorbar(reverse = TRUE)) +
      ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
      scale_size(range = c(3, 8))
  }
}
