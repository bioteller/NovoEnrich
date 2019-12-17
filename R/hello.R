#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
ne_env <- function(){
library(clusterProfiler,quietly = T)
library("data.table",quietly = T)
library(stringr,quietly = T)
library(ggplot2,quietly = T)
library(dplyr,quietly = T)
library(tidyverse,quietly = T)
}
testGO <- "../NoveSmart_US/NovoSmart_US.nr/data/classification/GO_classification.xls"
testKG <- "../NoveSmart_US/NovoSmart_US.nr/data/classification/KEGG_classification.xls"
DGE <- "../NoveSmart_US/NovoSmart_US.nr/data/DGE/ToDP00vsToDV00.DEG.xls"

# gset
ne_gset <- function(classificationfile) {
  GO_class <- fread(classificationfile,sep = "\t")
  TERM2GENE <- list()
  TERM2NAME <- list()
  for (i in unique(GO_class$`GO Term (Lev1)`)){
    temp <- GO_class[GO_class$`GO Term (Lev1)` == i,c("GO ID (Lev4)","Gene List")]
    m <- lapply(unlist(temp[,"Gene List"]),FUN = function(x) unlist(str_split(x,pattern=",")))
    names(m) <- temp$`GO ID (Lev4)`
    TERM2GENE[[i]] <- data.frame(Term = substr(names(unlist(m)),1,10),Gene = unlist(m))
    TERM2NAME[[i]] <- GO_class[GO_class$`GO Term (Lev1)` == i,c("GO ID (Lev4)","GO Term (Lev4)")]
  }
  return(list(TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME))
}

ne_kset <- function(classificationfile){
  KEGG_class <- fread(classificationfile,sep = "\t")
  m <- lapply(unlist(KEGG_class[,"Gene IDs"]),FUN = function(x) unlist(str_split(x,pattern=",")))
  names(m) <- KEGG_class$`Pathway ID`
  TERM2GENE <- data.frame(Term = substr(names(unlist(m)),1,7),Gene = unlist(m))
  rownames(TERM2GENE) <- seq(1,nrow(TERM2GENE),1)
  TERM2NAME <- KEGG_class[,c("Pathway ID","KEGG Pathway")]

  return(list(TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME))
}

# read dge
ne_dge  <- function(DGE_file,lg2fc=0,include="ALL",pval=1,padj=1) {
  dge <- fread(DGE_file,sep = "\t",data.table = F)
  dge <- dge[which(dge$pval < pval & dge$padj < padj),]
  gene_list <- unlist(log2((dge[,2]+1)/(dge[,3]+1)))
  names(gene_list) <- dge$gene_id
  gene_list <- gene_list[order(gene_list,decreasing = T)]
  if (include == "ALL") {
    gene_list <- gene_list[which(abs(gene_list) > lg2fc)]
  }else if (include =="UP") {
    gene_list <- gene_list[which(gene_list > lg2fc)]
  }else if (include == "DOWN"){
    gene_list <- gene_list[which(gene_list < -lg2fc)]
  }
  return(gene_list)
}

# enrich
ne_enrichGO <- function(gset,dge){
  enrichr_GO <- list()
  TERM2GENE <- gset$TERM2GENE
  TERM2NAME <- gset$TERM2NAME
  for (i in names(TERM2GENE)){
    enrichr_GO[[i]] <- enricher(gene = names(dge),TERM2GENE = TERM2GENE[[i]],TERM2NAME = TERM2NAME[[i]],pvalueCutoff = 1,qvalueCutoff = 1)
    enrichr_GO[[i]]@result$ontology <- i
    enrichr_GO[[i]]@result$gene_ratio <- unlist(lapply(enrichr_GO[[i]]@result$GeneRatio, function(x) eval(parse(text=x))))
  }

  return(enrichr_GO)
}

ne_enrichKEGG <- function(kset,dge){
  TERM2GENE <- kset$TERM2GENE
  TERM2NAME <- kset$TERM2NAME
  enrichr_KEGG <- enricher(gene = names(dge),TERM2GENE = TERM2GENE,TERM2NAME = TERM2NAME,pvalueCutoff = 1,qvalueCutoff = 1)
  enrichr_KEGG@result$ontology <- "KEGG"
  enrichr_KEGG@result$gene_ratio <- unlist(lapply(enrichr_KEGG@result$GeneRatio, function(x) eval(parse(text=x))))
  return(enrichr_KEGG)
}


# graph

Goplot_prepare <- function(dataset,
                           ont=c("Biological Process","Cellular Component","Molecular Function"),
                           n = c(10,10,10),s_name="p.adjust",desc=T){
  tmp = data.frame()
  dataset <- as.data.frame(dataset)
  dataset <- dataset[order(dataset[,s_name],decreasing = desc),]
  for (i in ont){
    tmp <- rbind(tmp,dataset[dataset[,"ontology"] ==i,][0:n[which(ont == i)],])
  }
  return(tmp)
}

ne_goplot <- function(dataset,
                      title="",
                      ont=c("Biological Process","Cellular Component","Molecular Function"),
                      y_name="p.adjust",
                      x.size=7,
                      # x_name="Description",
                      sep=T,
                      n = c(10,10,10),
                      fill="ontology",
                      lab_fix=10,
                      desc=F){
  d <- rbindlist(lapply(dataset,function(x) as.data.frame(x)))
  d <- Goplot_prepare(d,n=n,ont=ont,desc=desc,s_name = y_name)
  d <- d[order(d$p.adjust,decreasing = T),]
  d$Description2 <- d$Description
  if( lab_fix !=0 ){
    d$Description2 <- ifelse(nchar(d$Description) > lab_fix, paste0(str_sub(d$Description,1,lab_fix),"..."),d$Description)}

  d$Description <- factor(d$Description, levels = rev(d$Description))
  jkl <- sym(y_name)
  d %>%
    mutate(name=fct_reorder(Description,!!jkl,.desc = T)) %>%
    group_by(ontology) %>%
    ggplot(aes(x=Description,y=!!jkl,fill=!!sym(fill))) +
    geom_bar(stat = "identity",width = 0.8) +
    scale_x_discrete(breaks=d$Description, labels = d$Description2)+
    ggtitle(title)+
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1,size = x.size),
          axis.line = element_line(),
          strip.text =  element_text(face = "bold",size = 7,colour = "white"),
          strip.background = element_rect(fill = rgb(0.2,0.2,0.2,0.8),colour = "black"),
          plot.margin = unit(c(1,1,1,1),units = "cm"),
          legend.position="right",
          plot.title = element_text(hjust = 0.5)) -> ph
          if (sep) { ph <- ph + facet_grid(.~ontology,scales="free_x",space = "free_x")}
  return(ph)
}
ne_plot <- function(dataset,
                      title="",
                      y_name="p.adjust",
                      # x_name="Description",
                      x.size=8,
                      n = 10,
                      fill="p.adjust",
                      lab_fix=10,
                      desc=F){
  if (class(dataset) == "list"){dataset <- rbindlist(lapply(dataset,function(x) as.data.frame(x)))}
  d <- dataset[order(dataset$p.adjust,decreasing = T),]
  d$Description2 <- d$Description
  if( lab_fix !=0 ){
    d$Description2 <- ifelse(nchar(d$Description) > lab_fix, paste0(str_sub(d$Description,1,lab_fix),"..."),d$Description)}

  d$Description <- factor(d$Description, levels = rev(d$Description))
  jkl <- sym(y_name)
  d %>%
    mutate(name=fct_reorder(Description,!!jkl,.desc = T)) %>%
    top_n(n) %>%
    ggplot(aes(x=Description,y=!!jkl,fill=!!sym(fill))) +
    geom_bar(stat = "identity",width = 0.8) +
    scale_x_discrete(breaks=d$Description, labels = d$Description2)+
    ggtitle(title)+
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1,size = x.size),
          axis.line = element_line(),
          strip.text =  element_text(face = "bold",size = 7,colour = "white"),
          strip.background = element_rect(fill = rgb(0.2,0.2,0.2,0.8),colour = "black"),
          plot.margin = unit(c(1,1,1,1),units = "cm"),
          legend.position="right",
          plot.title = element_text(hjust = 0.5)) -> ph
  return(ph)
}


if(F) {
  testGO <- "../NoveSmart_US/NovoSmart_US.nr/data/classification/GO_classification.xls"
  testKG <- "../NoveSmart_US/NovoSmart_US.nr/data/classification/KEGG_classification.xls"
  DGE <- "../NoveSmart_US/NovoSmart_US.nr/data/DGE/ToDP00vsToDV00.DEG.xls"

  packageVersion("NovoEnrich")

  gset <- ne_gset(classificationfile = testGO)
  gde <- ne_dge(DGE_file = DGE,lg2fc = 1,)
  dge <- gde
  kset <- ne_kset(classificationfile = testKG)

  GOR <- ne_enrichGO(gset = gset,dge = gde)
  KOR <- ne_enrichKEGG(kset = kset,dge = gde)
barplot(GOR$`Molecular Function`,)
  ne_plot(KOR,title = "KEGG enrichment",n=30,lab_fix = 30,y_name = "Count",fill="Count")
  ne_goplot(GOR,title = "GO Enrichment",n=c(10,10,8),fill="ontology",lab_fix = 0)

  GSEAOR <- GSEA(gde,TERM2GENE = gset$TERM2GENE[[2]],TERM2NAME = gset$TERM2NAME[[2]],pvalueCutoff = 1)
  colnames(GSEAOR@result)
  ridgeplot(GSEAOR)
  gseaplot(GSEAOR,geneSetID = "GO:0046483")
  otype(GSEAOR)
  #GSEAOR@result
  }
if (require("NovoEnrich")) {
  devtools::install_github("bioteller/NovoEnrich")
  library(NovoEnrich)
}


