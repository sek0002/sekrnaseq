#' Return table of filtered and TMM (edgeR) normalised CPMs
#'
#' This function takes raw countdata+sampleinfo+group_pram and filters on minimum CPM and sample replicates. The genekeytype, "ensembl" or "entrez" and species, "human" or "mouse" must be specified.
#' @param countdata Raw counts table with entrez or ensembl IDs (output from merged bam files)
#' @param sampleinfo Sampleinfo file (see template)
#' @param group_param select group_param column to compare from sampleinfo
#' @param thresh Minimum CPM threshold to exclude
#' @param keep Minimum number of samples below thresh before exclusion (value usually set to number of replciates per group)
#' @param species "human" or "mouse"
#' @param isEnsembl TRUE or FALSE, if false defaults to entrezid
#' @param log Set to true to output normalised LogCPM
#' @param human_ensembldb
#' @param mouse_ensembldb
#' @return Table of normalised and filtered CPMs with converted gene symbols
#'
#' @import edgeR
#' @importFrom data.table setDT
#' @import DESeq2
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import EnsDb.Mmusculus.v79
#' @import EnsDb.Hsapiens.v75
#' @importFrom AnnotationDbi mapIds
#' @export
### Normalized counts

CPMtable<- function(countdata,sampleinfo, thresh=0.5, keep, group_param, species, isEnsembl="FALSE" ,log=F, human_ensembldb= "EnsDb.Hsapiens.v75", mouse_ensembldb= "EnsDb.Mmusculus.v79"){
#
#   group_param1 <- eval(parse(text=paste0("~",group_param,sep="")))

  if (isEnsembl==F){
    print("Defaulting to EntrezID")}
  if (human_ensembldb=="EnsDb.Hsapiens.v75"){
    print("Defaulting to EnsDb.Hsapiens.v75, hg19, GRCh37")}
  if (mouse_ensembldb=="EnsDb.Mmusculus.v79"){
    print("Defaulting to EnsDb.Mmusculus.v79")}

  myCPM <- edgeR::cpm(countdata)
  threshold <- myCPM >= thresh
  keep_rep <- rowSums(threshold) >= keep
  counts.keep <- countdata[keep_rep,]
  y <- edgeR::DGEList(counts.keep)
  nonlogcounts <- data.frame(edgeR::cpm(y, log=log, normalized.lib.sizes=T))
  logcounts <- edgeR::cpm(y, log=T)
  logcounts1 <-data.frame(logcounts)
  colSums(countdata)
  table(rowSums(threshold))
  dim(counts.keep)
  plot(myCPM[,1],countdata[,1], ylim=c(0,50), xlim=c(0,50),abline(v=thresh, h=10, col= "blue"))
  barplot(y$samples$lib.size,names=colnames(y),las=2, main = "Barplot of library sizes")

  y1 <- calcNormFactors(y)
  ### normalize using DEseq2 (DEseq)
  coldata<- data.frame(row.names=colnames(counts.keep))
  for (i in 2:ncol(sampleinfo)){
    coldata[colnames(sampleinfo[i])]<- factor(sampleinfo[,i])
  }

  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=eval(parse(text=paste0("~",group_param,sep=""))))
  dds <- DESeq(dds)
  rld <- rlogTransformation(dds)
  normalized.counts <- as.data.frame(counts(dds, normalized=TRUE ))


  if(isEnsembl=="FALSE") {
    godb_database<- switch(species, human = org.Hs.eg.db, mouse = org.Mm.eg.db)
    idtype <- "ENTREZID"
    gene.ids <- mapIds(godb_database, keys=rownames(y1),
                       keytype=idtype, column="SYMBOL")
    y1$genes <- data.frame(keytype=rownames(y1), SYMBOL=gene.ids)} else {

      godb_database<- switch(species, human = get(human_ensembldb), mouse = get(mouse_ensembldb))
      idtype <- "ENSEMBL"
      gene.ids <- mapIds(godb_database, keys=rownames(y1),
                         keytype=idtype, column="SYMBOL")
      y1$genes <- data.frame(keytype=rownames(y1), SYMBOL=gene.ids)
      # listMarts(host = 'http://grch37.ensembl.org')
      # ensemblarchived <- useMart(host = 'http://grch37.ensembl.org', biomart= "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
      # ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      # res <- getBM(c("hgnc_symbol","ensembl_gene_id"), filters = "ensembl_gene_id", values = rownames(y1), mart = ensemblarchived)
      # resunique<- res[!duplicated(res$ensembl_gene_id),]
      # y1$genes <- switch(species, human = data.frame(keytype=rownames(y1), SYMBOL=resunique[,"hgnc_symbol"]), mouse = data.frame(keytype=rownames(y1), SYMBOL=resunique[,"mgi_symbol"]))
    }
  return(nonlogcounts)  }
