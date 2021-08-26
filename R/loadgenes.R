#' loadgenes from MsigDB genesets
#'
#'
#' @param path1 name of pathway
#' @param species a=mouse or b=human

#' @return Table of normalised and filtered DGElist objects for limma/voom and EdgeR
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


loadgenes <- function(path1,species){
  file_path <- switch(species,
                      a = "Genesets/mouse/",
                      b = "Genesets/human/")
  file_list <- list.files(path=file_path)
  dat <-vector(length = length(file_list))
  longlist <- c()
  for (i in 1:length(file_list)){
    listgenesets<- lapply(paste0(file_path,file_list[i]),function(x) get(load(x)))
    dat[i]<- listgenesets}

  for (i in 1:length(file_list)){
    longlist <- c(longlist, dat[[i]])
  }
  ann <- switch(species,
                a = AnnotationDbi::select(org.Mm.eg.db, keys=longlist[[path1]], columns=c("SYMBOL")),
                b = AnnotationDbi::select(org.Hs.eg.db, keys=longlist[[path1]], columns=c("SYMBOL")))
  return(ann)
}
