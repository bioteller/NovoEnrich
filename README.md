# NovoEnrich V0.99 alpha #
#### NovoEnrich for TR project ####

### Author

**Xiaobo Li** lixiaobo@novogene.com

### Description
GO, KEGG and other enrichment analysis for TR project. Version 0.5. 2019-12-3

Examples
## NovoEnrich ##


### HOW TO INSTALL **NovoEnrich** ###
```R
devtools::install_github("bioteller/NovoEnrich")
library(NovoEnrich)
```

### HOW TO  do enrich analysis via NovoEnrich ###
```R
rm(list = ls())
library(NovoEnrich)

testGO <- "../NoveSmart_US/NovoSmart_US.nr/data/classification/GO_classification.xls"
testKG <- "../NoveSmart_US/NovoSmart_US.nr/data/classification/KEGG_classification.xls"
DGE <- "../NoveSmart_US/NovoSmart_US.nr/data/DGE/ToDP00vsToDV00.DEG.xls"

packageVersion("NovoEnrich")

gset <- ne_gset(classificationfile = testGO)
gde <- ne_dge(DGE_file = DGE)
kset <- ne_kset(classificationfile = testKG)

GOR <- ne_enrichGO(gset = gset,dge = gde)
KOR <- ne_enrichKEGG(kset = kset,dge = gde)

ne_plot(KOR,title = "KEGG enrichment",n=25,lab_fix = 30)
ne_goplot(GOR,title = "GO Enrichment",n=c(10,10,3),fill="p.adjust")
```
> NovoEnrich 0.90 aftersale of TR
