#' Limma/voom and EdgeR differential analysis
#'
#' DEanalysis with Limma/voom and EdgeR, with outputs.
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




DE_analysis_lrt <- function(countdata,sampleinfo, thresh=0.5, keep, group_param, exp_name, design, DEgroup, cont.matrix, species, genekeytype, archived="FALSE", plots=F) {
  list_of_folders<- c("DE","export","GSEA","Heatmaps","MDS","Volcano")

  for(x in list_of_folders){
    dir.create(path = paste0(x,"/",exp_name,"/interactive"), recursive = TRUE)
  }

  dir.create(path = paste0("GSEA/",exp_name,"/H"), recursive = TRUE)
  dir.create(path = paste0("GSEA/",exp_name,"/C2"), recursive = TRUE)
  dir.create(path = paste0("GSEA/",exp_name,"/C3"), recursive = TRUE)
  dir.create(path = paste0("GSEA/",exp_name,"/C5"), recursive = TRUE)
  dir.create(path = paste0("GSEA/",exp_name,"/C2_CGP"), recursive = TRUE)
  dir.create(path = paste0("GSEA/",exp_name,"/C3_TFT"), recursive = TRUE)
  dir.create(path = paste0("GSEA/",exp_name,"/C7"), recursive = TRUE)
  dir.create(path = paste0("GSEA/",exp_name,"/CHEA_TF"), recursive = TRUE)
  dir.create(path = paste0("GSEA/",exp_name,"/Lit_chipseq"), recursive = TRUE)


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

  ### voom/limma for all possible combinations
  v<- voom(y1,design,plot=TRUE)
  limmafit <- lmFit(v)
  limmafit.cont <- contrasts.fit(limmafit, cont.matrix)
  limmafit.cont <- eBayes(limmafit.cont)
  summa.limmafit <- decideTests(limmafit.cont, p.value=0.05)
  voom_file<- list()

  if(plots==T){
    for (x in colnames(cont.matrix)) {
      limma.res <- topTable(limmafit.cont, coef=x, sort.by='p', n=Inf)
      limma.res<- limma.res[order(as.numeric(row.names(limma.res))),]
      summa.limmafit.reorder <- summa.limmafit[order(as.numeric(row.names(summa.limmafit))),]
      limmafit.cont$limma.res <- limma.res
      summary(summa.limmafit)
      table(summa.limmafit)

      limmafit.cont$limma.res$Threshold <- ifelse(limmafit.cont$limma.res$adj.P.Val <=0.05 & limmafit.cont$limma.res$logFC >0,"Up", ifelse(limmafit.cont$limma.res$adj.P.Val <=0.05 & limmafit.cont$limma.res$logFC <0, "Down", "notsig"))

      limmafit.cont$limma.res$genes$SYMBOL <- limmafit.cont$limma.res$Gene.Name

      topgenes<- limma.res %>% dplyr::select("SYMBOL","adj.P.Val","logFC") %>% dplyr::filter(logFC >=1 & adj.P.Val<=0.05|logFC<=-1 & adj.P.Val<=0.05) %>% na.omit()

      maplotniceqlf(limmafit.cont$limma.res, xaxis= "AveExpr", yaxis= "logFC", paste0(x,"_voom_MA",sep=""), exp_name)
      nicevolcqlf(limmafit.cont$limma.res, xaxis= "logFC", yaxis= "adj.P.Val", paste0(x,"_Voom_volc",sep=""), exp_name)

      if(length(unique(limmafit.cont$limma.res$Threshold))==3){
        glXYPlot(x=limmafit.cont$limma.res$AveExpr, y=limmafit.cont$limma.res$logFC,
                 xlab="AveExpr", ylab="logFC", main=paste(x),groups=sampleinfo[,group_param], status=limmafit.cont$limma.res$Threshold,
                 anno=limmafit.cont$limma.res, side.main="SYMBOL", path= paste0("Volcano/",exp_name,"/interactive"),folder=paste0(x,"_MA"))

        glXYPlot(x=limmafit.cont$limma.res$logFC, y=-log10(limmafit.cont$limma.res$adj.P.Val),
                 xlab="logFC", ylab="-Log10FDR", main=paste(x),groups=sampleinfo[,group_param], status=limmafit.cont$limma.res$Threshold,
                 anno=limmafit.cont$limma.res, side.main="SYMBOL", path= paste0("Volcano/",exp_name,"/interactive"),folder=paste0(x,"_volc"))

      }else{

        glXYPlot(x=limmafit.cont$limma.res$AveExpr, y=limmafit.cont$limma.res$logFC,
                 xlab="AveExpr", ylab="logFC", main=paste(x),groups=sampleinfo[,group_param],
                 anno=limmafit.cont$limma.res, side.main="SYMBOL", path= paste0("Volcano/",exp_name,"/interactive"),folder=paste0(x,"_MA"))

        glXYPlot(x=limmafit.cont$limma.res$logFC, y=-log10(limmafit.cont$limma.res$adj.P.Val),
                 xlab="logFC", ylab="-Log10FDR", main=paste(x),groups=sampleinfo[,group_param],
                 anno=limmafit.cont$limma.res, side.main="SYMBOL", path= paste0("Volcano/",exp_name,"/interactive"),folder=paste0(x,"_volc"))

      }
      voom_file[x]<- list(limmafit.cont)

      write.table(limmafit.cont$limma.res, file=paste0("DE/",exp_name,"/voom_output_",x,".tsv",sep=""), sep="\t", row.names = F)
    }}else{
      for (x in colnames(cont.matrix)) {
        limma.res <- topTable(limmafit.cont, coef=x, sort.by='p', n=Inf)
        limma.res<- limma.res[order(as.numeric(row.names(limma.res))),]
        summa.limmafit.reorder <- summa.limmafit[order(as.numeric(row.names(summa.limmafit))),]
        limmafit.cont$limma.res <- limma.res
        summary(summa.limmafit)
        table(summa.limmafit)

        limmafit.cont$limma.res$Threshold <- ifelse(limmafit.cont$limma.res$adj.P.Val <=0.05 & limmafit.cont$limma.res$logFC >0,"Up", ifelse(limmafit.cont$limma.res$adj.P.Val <=0.05 & limmafit.cont$limma.res$logFC <0, "Down", "notsig"))

        limmafit.cont$limma.res$genes$SYMBOL <- limmafit.cont$limma.res$Gene.Name

        topgenes<- limma.res %>% dplyr::select("SYMBOL","adj.P.Val","logFC") %>% dplyr::filter(logFC >=1 & adj.P.Val<=0.05|logFC<=-1 & adj.P.Val<=0.05) %>% na.omit()

        voom_file[x]<- list(limmafit.cont)

        write.table(limmafit.cont$limma.res, file=paste0("DE/",exp_name,"/voom_output_",x,".tsv",sep=""), sep="\t", row.names = F)

      }
    }


  # ###EdgeR qlf
  y1 <- estimateDisp(y1, design, robust=TRUE)
  y1$common.dispersion
  plotBCV(y1)
  edgeRfit <- glmFit(y1, design) ### no robust= like DEGUST
  qlf_file<- list()

  if(plots==T){

    for (x in colnames(cont.matrix)){

      qlf <- glmLRT(edgeRfit, contrast=cont.matrix[,x])
      summary(decideTests(qlf))
      plotMD(qlf)

      EdgeR.res <- topTags(qlf, n=Inf)$table
      EdgeR.summa <- decideTests(qlf, p.value=0.05)
      EdgeR.summa.reorder<- EdgeR.summa[,1][match(rownames(EdgeR.res), rownames(EdgeR.summa))]
      EdgeR.res.sig <- topTags(qlf, n=Inf, p=0.05)$table
      qlf$EdgeR.res<- EdgeR.res

      # ## set lfc threshold value, FDR cutoff of 5%
      # tr <- glmTreat(edgeRfit, contrast=cont.matrix[,"hA1ROV_CtrlOV"], lfc=log2(1.2))
      # topTags(tr)
      #  summary(decideTests(tr))
      #  plotMD(tr)

      qlf$EdgeR.res$Threshold <- ifelse(qlf$EdgeR.res$FDR <=0.05 & qlf$EdgeR.res$logFC >0,"Up", ifelse(qlf$EdgeR.res$FDR <=0.05 & qlf$EdgeR.res$logFC <0, "Down", "notsig"))


      maplotniceqlf(qlf$EdgeR.res, xaxis= "logCPM", yaxis= "logFC", paste0(x,"_lrt_MA",sep=""), exp_name)

      nicevolcqlf(qlf$EdgeR.res,  xaxis= "logFC", yaxis= "FDR", paste0(x,"_lrt_volc",sep=""), exp_name)

      glMDPlot(qlf, main=paste(x), groups= sampleinfo[,group_param],
               status = EdgeR.summa.reorder, counts = y1$counts, side.main="SYMBOL",
               anno=qlf$EdgeR.res,path= paste0("Volcano/",exp_name,"/interactive"), folder=paste0(x,"_lrt_MA"))

      glXYPlot(x=qlf$EdgeR.res$logFC, y=-log10(qlf$EdgeR.res$FDR),
               xlab="logFC", ylab="-Log10FDR", main=paste(x),groups=sampleinfo[,group_param], status=EdgeR.summa.reorder,
               anno=qlf$EdgeR.res, side.main="SYMBOL", path= paste0("Volcano/",exp_name,"/interactive"),folder=paste0(x,"_lrt_volc"))


      qlf_file[x]<- list(qlf)

      write.table(qlf$EdgeR.res, file=paste0("DE/",exp_name,"/lrt_output_",x,".tsv",sep=""), sep="\t", row.names = F)
    }
  }else{
    for (x in colnames(cont.matrix)){

      qlf <- glmLRT(edgeRfit, contrast=cont.matrix[,x])
      summary(decideTests(qlf))
      plotMD(qlf)

      EdgeR.res <- topTags(qlf, n=Inf)$table
      EdgeR.summa <- decideTests(qlf, p.value=0.05)
      EdgeR.summa.reorder<- EdgeR.summa[,1][match(rownames(EdgeR.res), rownames(EdgeR.summa))]
      EdgeR.res.sig <- topTags(qlf, n=Inf, p=0.05)$table
      qlf$EdgeR.res<- EdgeR.res

      # ## set lfc threshold value, FDR cutoff of 5%
      # tr <- glmTreat(edgeRfit, contrast=cont.matrix[,"hA1ROV_CtrlOV"], lfc=log2(1.2))
      # topTags(tr)
      #  summary(decideTests(tr))
      #  plotMD(tr)

      qlf$EdgeR.res$Threshold <- ifelse(qlf$EdgeR.res$FDR <=0.05 & qlf$EdgeR.res$logFC >0,"Up", ifelse(qlf$EdgeR.res$FDR <=0.05 & qlf$EdgeR.res$logFC <0, "Down", "notsig"))
      qlf_file[x]<- list(qlf)

      write.table(qlf$EdgeR.res, file=paste0("DE/",exp_name,"/lrt_output_",x,".tsv",sep=""), sep="\t", row.names = F)
    }
  }

  ### MDS plot of normalized data set to reactive parameter sampleinfo(column)
  col.status <- c(1:length(levels(sampleinfo[,group_param])))[sampleinfo[,group_param]]
  plotMDS(y1,col=col.status)

  ###intereactive MDS plot using Glimma/ EdgeR normalizatio
  labels <- paste(sampleinfo$Column)
  glMDSPlot(y1, labels=labels, groups=groupDE, folder=paste0("MDS/",exp_name))

  ### DEseq normalization PCA plot Set to subset by reactive parameter sampleinfo(param_no.), ntop
  ntop <- 100
  Pvars <- rowVars(assay(rld))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                        length(Pvars)))]
  PCA <- prcomp(t(assay(rld)[select, ]), scale = F)

  DEseqPCA<-ggbiplot(PCA, groups=coldata$param_1, obs.scale = 0.2,var.axes= FALSE, alpha = 0)+
    scale_colour_brewer(type="qual", palette=2)+
    geom_point(alpha=0.4, size=5, aes(col=groupDE))+
    theme_classic()+
    theme(legend.position = "bottom")+
    theme(text = element_text(size=15))
  DEseqPCA
  ggsave(file=paste0("MDS/",exp_name,"/PCA.svg"), units="in", plot=DEseqPCA, width=5, height=5)



  # Plot the heatmap
  logcountsnorm <- cpm(y1, log=TRUE)
  var_genes <- apply(logcountsnorm, 1, var)
  select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
  highly_variable_lcpm <- logcountsnorm[select_var,]
  select_lowvar <- names(sort(var_genes, decreasing=FALSE))[1:500]
  lowly_variable_lcpm <- logcountsnorm[select_lowvar,]

  mypalette <- brewer.pal(11,"RdYlBu")
  morecols <- colorRampPalette(mypalette)
  mypalette1 <- brewer.pal(11,"PiYG")
  morecols1 <- colorRampPalette(mypalette1)
  heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",scale="row")
  heatmap.2(lowly_variable_lcpm,col=rev(morecols1(50)),trace="none", main="Top 500 least variable genes across samples",scale="row")

  ###RLE plots
  gene_medians <- apply(logcounts, 1, median)
  # Calculate the deviation by subtracting the medians from the log counts
  rle <- sweep(logcounts, 1, STATS=gene_medians, FUN="-")
  # Get normalised log counts
  norm_counts <- cpm(y1, normalized.lib.sizes=TRUE, log=TRUE)
  # Calculate the gene medians for the normalised log counts
  norm_gene_medians <- apply(norm_counts, 1, median)
  # Calculate the deviation by subtracting the medians from the log counts
  norm_rle <- sweep(norm_counts, 1, STATS=norm_gene_medians, FUN="-")

  # Plot the boxplots
  par(mfrow=c(1,2))
  boxplot(rle, outline=FALSE, xlab="", ylab="Relative log expression", las=2)
  abline(h=0, col="red", lty="dotted")
  title("RLE plot (unnormalised)")
  boxplot(norm_rle, outline=FALSE, xlab="", ylab="Relative log expression", las=2)
  abline(h=0, col="red", lty="dotted")
  title("RLE plot (normalised)")

  ### Plotting normalized figures (x= reactive col number)
  par(mfrow=c(1,2))
  plotMD(y,column = 3,main = "Unormalized")
  abline(h=0,col="grey")
  plotMD(y1,column = 3, main="TMM Unormalized")
  abline(h=0,col="grey")

  ### LogCPMs expression box-whisker
  boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2, main="Boxplots of logCPMs (unnormalised)")
  abline(h=median(logcounts),col="blue")
  boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
  abline(h=median(v$E),col="blue")

  data_out<- list(voom_file=voom_file,qlf_file=qlf_file)
  return(data_out)
}
