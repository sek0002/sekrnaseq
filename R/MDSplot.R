#' Limma/voom and EdgeR differential analysis, MDS plot
#'
#' MDSplot with Limma/voom and EdgeR, with outputs.
#'
#' @param countdata Raw counts table with entrez or ensembl IDs (output from merged bam files)
#' @param sampleinfo Sampleinfo file (see template)
#' @param group_param select group_param column to compare from sampleinfo
#' @param thresh Minimum CPM threshold to exclude
#' @param keep Minimum number of samples below thresh before exclusion (value usually set to number of replciates per group)
#' @param species "human" or "mouse"
#' @param genekeytype "entrez" or "ensembl"
#' @param archived Should be set to FALSE, unless using archived Human Ensembl database
#' @param log Set to true to output normalised LogCPM
#' @param exp_name Name to call experiment
#' @param design Design matrix (use model.matrix to generate, see template)
#' @param DEgroup Specify comparison of interest for DE
#' @param cont.matrix Contrast matrix (use makeContrast to generate, see template)
#' @param plots Set to TRUE to generate TIFF and PDF volcano and MAplots in Volcano/exp_name
#' @param qlftest Default=TRUE, set to FALSE for edgeR LRT test
#' @param ntop default top 100, for MDSplot function only
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
#' @import ggbiplot
#' @import matrixStats
#' @export




MDSplot <- function(countdata,sampleinfo, thresh=0.5, keep, group_param, exp_name, design, DEgroup, species, archived="FALSE",genekeytype, plots=F, ntop=100) {
  myCPM <- cpm(countdata)
  threshold <- myCPM >= thresh
  keep_rep <- rowSums(threshold) >= keep
  counts.keep <- countdata[keep_rep,]
  y <- DGEList(counts.keep)
  logcounts <- cpm(y, log=TRUE)
  colSums(countdata)
  table(rowSums(threshold))
  dim(counts.keep)
  plot(myCPM[,1],countdata[,1], ylim=c(0,50), xlim=c(0,50),abline(v=thresh, h=10, col= "blue"))
  barplot(y$samples$lib.size,names=colnames(y),las=2, main = "Barplot of library sizes")


  ## normlalize for composition bias using EdgeR
  y1 <- calcNormFactors(y)

  if(archived=="FALSE") {
    godb_database<- switch(species, human = org.Hs.eg.db, mouse = org.Mm.eg.db)
    idtype <- switch(genekeytype, ensembl = "ENSEMBL", entrez = "ENTREZID")
    gene.ids <- mapIds(godb_database, keys=rownames(y1),
                       keytype=idtype, column="SYMBOL")
    y1$genes <- data.frame(keytype=rownames(y1), SYMBOL=gene.ids)} else {
      listMarts(host = 'http://grch37.ensembl.org')
      ensemblarchived <- useMart(host = 'http://grch37.ensembl.org', biomart= "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      res <- getBM(c("hgnc_symbol","ensembl_gene_id"), filters = "ensembl_gene_id", values = rownames(y1), mart = ensemblarchived)
      resunique<- res[!duplicated(res$ensembl_gene_id),]
      y1$genes <- switch(species, human = data.frame(keytype=rownames(y1), SYMBOL=resunique[,"hgnc_symbol"]), mouse = data.frame(keytype=rownames(y1), SYMBOL=resunique[,"mgi_symbol"]))

    }

  ### normalize using DEseq2 (DEseq)
  coldata<- data.frame(row.names=colnames(counts.keep))
  for (i in 2:ncol(sampleinfo)){
    coldata[colnames(sampleinfo[i])]<- factor(sampleinfo[,i])
  }

  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=eval(parse(text=paste0("~",group_param,sep=""))))
  dds <- DESeq(dds)
  rld <- rlogTransformation(dds)

  ### MDS plot of normalized data set to reactive parameter sampleinfo(column)
  col.status <- c(1:length(levels(sampleinfo[,group_param])))[sampleinfo[,group_param]]
  plotMDS(y1,col=col.status)

  ### DEseq normalization PCA plot Set to subset by reactive parameter sampleinfo(param_no.), ntop
  ntop <- ntop
  Pvars <- rowVars(assay(rld))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                        length(Pvars)))]
  PCA <- prcomp(t(assay(rld)[select, ]), scale = T)

  DEseqPCA<-ggbiplot(PCA, groups=coldata$param_1, obs.scale = 0.2,var.axes= FALSE, alpha = 0)+
    scale_colour_brewer(type="qual", palette=2)+
    geom_point(alpha=0.4, size=7, aes(col=groupDE))+
    theme_classic()+
    theme(legend.position = "bottom")+
    theme(text = element_text(size=15))
  DEseqPCA

  if(plots==T){
    ###intereactive MDS plot using Glimma/ EdgeR normalizatio
    labels <- paste(sampleinfo$Column)
    glMDSPlot(y1, labels=labels, groups=groupDE, folder=paste0("MDS/",exp_name))

    ggsave(file=paste0("MDS/",exp_name,"/PCA_select.svg"), units="in", plot=DEseqPCA, width=5, height=5)
  }
}
