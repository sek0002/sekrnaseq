#' Easily convert between mouse, human and zebrafish genes
#'
#'
#'
#' @param x input list of genes
#'
#' @importFrom data.table setDT
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi mapIds
#' @import biomaRt
#' @return list of input genes with respective conversions
#' @export

convertMouseGeneList <- function(x){

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2)
  return(humanx)
}

#' @export
convertHumanGeneList <- function(x){

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2)
  return(humanx)
}

#' @export
convertZebratoHumanGeneList <- function(x){

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  zebra = useMart("ensembl", dataset = "drerio_gene_ensembl")
  genesV2 = getLDS(attributes = c("zfin_id_symbol","ensembl_gene_id"), filters = "ensembl_gene_id", values = x , mart = zebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2)
  return(humanx)
}
