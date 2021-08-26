#' GSEA_input scripts using fgsea algorithm
#'
#'
#' @param countdata Raw counts table with entrez or ensembl IDs (output from merged bam files)

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


gsea_input <- function(input,species,orderby="logFC",cat =NULL, subcat = NULL, nPerm=1000, minGSSize=2, maxGSSize=500, pvalueCutoff= 0.05, by="fgsea", plotpathway=NULL, expname, DEgroup, showCategory=20, customGS=NULL, plot=F, ...){
  if (is.null(input)){
    stop('no input gene list provided')}
  else if (is.null(cat)){
    print("H: Hallmarks, C1: Positional, C2: Curated, C3:Regulatory target, C4:Computational, C5: GOterm, C6: Oncogenic, C7: immunologic, more than one can be used... Visit: https://www.gsea-msigdb.org/gsea/msigdb/index.jsp")
  } else if (is.null(orderby)){
    print("ranked by logFC")
  }
  if (is.null(subcat)){
    print("Using entire geneset")
  }
  if (pvalueCutoff==1){
    print("output all enrichment results, limit to restrict")
  }
  if (by=="fgsea"){
    print("using fgsea, other options are DOSE")
  }
  if (is.null(plotpathway)){
    print("all pathways selected, Input plotpathway list for specific pathways")
  }


  if (is.null(customGS)){
    gseapathways <- switch(species,
                           human = msigdbr(species = "Homo sapiens"),
                           mouse = msigdbr(species = "Mus musculus"))
    if (is.null(subcat)) {
      gseapathways<-gseapathways%>% dplyr::filter(gs_cat == cat ) %>% dplyr::select(gs_name, gene_symbol)
    } else{
      gseapathways<-gseapathways%>% dplyr::filter(gs_cat == cat & gs_subcat %in% subcat) %>% dplyr::select(gs_name, gene_symbol)
    }} else {
      gseapathways<-customGS %>% dplyr::select(gs_name, gene_symbol)
    }


  gseaDat <- input

  ranks <- switch(orderby,
                  logFC = gseaDat$logFC,
                  pvalue = if(method=="edger"){
                    gseaDat$FDR
                  }else{gseaDat$adj.P.val}
  )


  names(ranks) <- as.character(gseaDat$SYMBOL)
  geneList <- sort(ranks,decreasing = T)
  barplot(geneList)
  head(geneList)

  results<-GSEA(geneList=geneList,
                TERM2GENE=gseapathways,
                nPerm = nPerm,
                verbose=TRUE,
                minGSSize =  minGSSize,
                maxGSSize = maxGSSize,
                pvalueCutoff = pvalueCutoff,
                by=by)

  results_all<-GSEA(geneList=geneList,
                    TERM2GENE=gseapathways,
                    nPerm = nPerm,
                    verbose=TRUE,
                    minGSSize =  minGSSize,
                    maxGSSize = maxGSSize,
                    pvalueCutoff = 1,
                    by=by)

  write.table(results_all, file=paste0("GSEA/",expname,"/",DEgroup,"_","GSEA_results_all.tsv"),row.names = FALSE, sep="\t")

  write.table(results, file=paste0("GSEA/",expname,"/",DEgroup,"_","GSEA_results.tsv"),row.names = FALSE, sep="\t")

  if (plot==T){
    if (is.null(plotpathway)){
      for (genesetID in results$ID) {
        plot<- gseaplot2(results, geneSetID = paste0(genesetID), title = paste0(genesetID))
        ggsave(file=paste0("GSEA/",expname,"/",DEgroup,"_",genesetID,".pdf", sep=""), plot=plot, units="in", width=7, height=6)
      }} else {
        for (genesetID in plotpathway)
          plot<- gseaplot2(results, geneSetID = genesetID, title = paste0(genesetID))
        ggsave(file=paste0("GSEA/",expname,"/",DEgroup,"_",genesetID,".pdf", sep=""), plot=test, units="in", width=7, height=6)
      }}
  return(results)

}
