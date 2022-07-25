#' Easily convert human to mouse genes
#'
#'
#'
#' @param x input list of human genes
#'
#' @importFrom data.table setDT
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi mapIds
#' @import biomaRt
#' @return list of input genes with respective conversions
#' @export


convertHumanGeneList <- function(x){

  human = useEnsembl("ensembl",dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse = useEnsembl("ensembl",dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2)
  return(humanx)
}

