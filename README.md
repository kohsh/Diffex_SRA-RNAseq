# ✒️ SRA-Diffex-analysis

my Guthub repo include pipeline for running Differential Gene Expression (Diffex) analysis on [SRA](https://www.ncbi.nlm.nih.gov/sra) fastq files using [edgeR](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) tool. The pipeline also include all the required pipeline and their dependent packages for sample/ gene Quality checks, Normalization, PCA analysis, ANOVA analysis, Batch effect corrections,topGO and visulaization which are all written in an [R](https://github.com/rstudio/rstudio) programming language. 

## ⚙️ Technologies & Tools

![](https://img.shields.io/badge/Code-RScript-informational?style=flat&logo=<#FF6000>&logoColor=white&color=2bbc8a)
![](https://img.shields.io/badge/Tools-Rstudio-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)
![](https://img.shields.io/badge/Tools-GitHub-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)
![](https://img.shields.io/badge/Tools-SRAtoolkit-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)

```R
#!/usr/bin/Rscript

# load packages
library("XML")
library("xml2")
library("data.table")
library("VennDiagram")
library("RUVSeq")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library("org.Hs.eg.db")
library("topGO")
library("ggbiplot")
library("reshape")
library("gplots")
library("calibrate")
library("biomaRt")
library("sva")
library("ggplot2")
library("corrplot")
library("gage")
library("Rgraphviz")
library("pathview")
library("keep")
library("bootstrap")
library("edgeR")
library("dplyr")
library("shiny")
library("DT")
library("rafalib")
library("preprocessCore")
library("impute")
library("matrixStats")
library("splines")
library("foreach")
library("doParallel")
library("fastcluster")
library("dynamicTreeCut")
library("survival")
library("GO.db")
library("WGCNA")
library("flashClust")
library("bigmemory")
library("DESeq2")
require("ggpubr")
require("tidyverse")
require("Hmisc")
library("ggstatsplot")
library("hrbrthemes")
library("dslabs")
data(gapminder)
library("psych")
library("EnhancedVolcano")
library("viridis")


# Set working directory
setwd(path.expand("/Path/To/Data") )

# Read in the data
covar <- read.delim("meta_df.txt", sep='\t', header=TRUE, stringsAsFactors=FALSE)
reads <- read.delim("reads.txt", sep='\t', header=TRUE, stringsAsFactors=FALSE)

# Filtering out lowly expressed genes based on sample groups
group <- paste(c(rep("AD_PBMC",22), rep("PD_PBMC",46)))
y <- DGEList(counts= reads, group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE] 
y <- calcNormFactors(y,method = "TMM")
reads_norm <- cpm(y, normalized.lib.sizes = T)
write.csv(reads_norm,"Normalised_reads.csv")

# Unsupervised clustering to identify outliers
datExpr <- t(reads_norm)
datTraits <- dplyr::select(covar, c('Sample', 'Project_Numeric', 'Strategy_Numeric', 'Condition', 'Age', 'Sex', 'Braak', 'RIN'))
head(datTraits)
rownames(covar) = datTraits$Sample
# show that row names agree
table(rownames(covar)==rownames(datExpr)) #should return TRUE if datasets align correctly.
datTraits$Sample <- as.numeric(factor(datTraits$Sample))
datTraits$Project_Numeric <- as.numeric(factor(datTraits$Project_Numeric))
datTraits$Strategy_Numeric <- as.numeric(factor(datTraits$Strategy_Numeric))
datTraits$Condition <- as.numeric(factor(datTraits$Condition))
datTraits$Age <- as.numeric(factor(datTraits$Age))
datTraits$Sex <- as.numeric(factor(datTraits$Sex))
datTraits$Braak <- as.numeric(factor(datTraits$Braak))
datTraits$RIN <- as.numeric(factor(datTraits$RIN))
str(datTraits)

# Sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(t(datExpr),type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)

# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 #often -2.5

# the color vector indicates outlyingness (red)
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(datTraits,signed=FALSE))
str(traitColors)
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor,traitColors)
pdf("SampleDendro.pdf", width = 20, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree,groupLabels=names(datColors),colors=datColors,main="Sample Dendrogram and Trait Heatmap")
dev.off()

# Remove outlying samples from expression and trait data
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
# the following 2 lines differ from what is written in the book
datExpr=datExpr [!remove.samples,]
datTraits=datTraits[!remove.samples,]
# Recompute the sample network among the remaining samples
A=adjacency(t(datExpr),type="distance")
# Let's recompute the Z.k values of outlyingness
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)


# PCA analaysis(variation between data records(finding important genes))
rsem.pca <- prcomp(data.frame(t(reads_norm)), scale = TRUE)

pdf(file = "PCA_Condition_outliers.pdf")
ggbiplot(rsem.pca, choices = 1:2, obs.scale = 1, var.scale = 1, group = as.factor(covar$Condition), 
ellipse = TRUE, circle = TRUE, varname.size = 2, labels.size=3, var.axes = F, labels=covar$Sample) +
scale_color_manual(name="Condition", values=c("orange", "purple", "green")) +  
scale_shape_manual(name="Condition", values=c(17:19)) +
geom_point(aes(colour=as.factor(covar$Condition), shape=as.factor(covar$Condition)), size = 3) +
theme(legend.direction ="horizontal", 
      legend.position = "top")
dev.off()

pdf(file = "PCA_Project_outliers.pdf", width = 10, height = 8)
ggbiplot(rsem.pca, choices = 1:2, obs.scale = 1, var.scale = 1, group = as.factor(covar$Project_Numeric), 
ellipse = TRUE, circle = TRUE, varname.size = 2, labels.size=3, var.axes = F, labels=covar$Sample) +
scale_color_manual(name="Project_Numeric", values=c("orange", "purple", "green", "yellow", "red", "pink", "blue", "gray")) + 
scale_shape_manual(name="Project_Numeric", values=c(17:24)) +
geom_point(aes(colour=as.factor(covar$Project_Numeric), shape=as.factor(covar$Project_Numeric)), size = 3) +
theme(legend.direction ="horizontal", 
      legend.position = "top")
dev.off()

pdf(file = "PCA_Starategy_outliers.pdf")
ggbiplot(rsem.pca, choices = 1:2, obs.scale = 1, var.scale = 1, group = as.factor(covar$Strategy_Numeric), 
ellipse = TRUE, circle = TRUE, varname.size = 2, labels.size=3, var.axes = F, labels=covar$Strategy) +
scale_color_manual(name="Strategy_Numeric", values=c("orange", "purple", "green", "pink")) + 
scale_shape_manual(name="Strategy_Numeric", values=c(17:20)) +
geom_point(aes(colour=as.factor(covar$Strategy_Numeric), shape=as.factor(covar$Strategy_Numeric)), size = 3) +
theme(legend.direction ="horizontal", 
      legend.position = "top")
dev.off()

# ANOVA analysis
z=y$counts
m=melt(data.frame(z))
colnames(m) <- c("sample_ID","counts")
dim(m)
age <- rep(as.numeric(covar$Age, each=nrow(z)))
age
sex <- rep(as.numeric(covar$Sex, each=nrow(z)))
sex
dis <- rep(as.numeric(covar$Condition, each=nrow(z)))
dis
study <- rep(as.numeric(covar$Project_Numeric, each=nrow(z)))
study
strategy <- rep(as.numeric(covar$Strategy_Numeric, each=nrow(z)))
strategy
braak <- rep(as.numeric(covar$Braak, each=nrow(z)))
braak
RIN <- rep(as.numeric(covar$RIN, each=nrow(z)))
RIN
batch1 <- rep(as.numeric(covar$RUVgW_1, each=nrow(z)))
batch1
batch2 <- rep(as.numeric(covar$RUVgW_2, each=nrow(z)))
batch2
batch3 <- rep(as.numeric(covar$RUVrW_1, each=nrow(z)))
batch3
batch4 <- rep(as.numeric(covar$RUVrW_2, each=nrow(z)))
batch4
matrix <- data.frame(m, Condition=dis, Age=age, Sex=sex, Project_Numeric=study, Strategy_Numeric=strategy)
fit1 <- lm(counts ~ Condition + Age + Sex + Project_Numeric + Strategy_Numeric + braak + RIN, data=matrix) 
summary(fit1)
a <- anova(fit1)
pdf(file="Anovar.pdf")
maxval = 100
nfac <- length(a[,1])-1
barplot(100*a$"Sum Sq"[1:nfac]/sum(a$"Sum Sq"[1:nfac]),names.arg=rownames(a[1:nfac,]),ylim=c(0,maxval),las=3)
dev.off()

# limma removeBatchEffect to visualise effect of removing batch using a linear model
reads_forlim <- data.frame(reads_norm)
lim_batch <- covar$Project_Numeric
reads_lim <- removeBatchEffect(reads_forlim, batch=lim_batch)

rsem.pca <- prcomp(t(reads_lim), scale = TRUE)
# Batch effect pca 
pdf(file = "PCA_limma_logcorrect2_study.pdf", width = 15, height = 12)
ggbiplot(rsem.pca, choices = 1:2, obs.scale = 1, var.scale = 1,groups = as.factor(covar$Project_Numeric), ellipse = TRUE, circle = TRUE, var.axes = F,labels=covar$Sample) +
  scale_color_discrete(name="Group", h=c(0,360) + 20) +
  theme(legend.direction = 'vertical', legend.position = 'right')
dev.off()

# limma removeBatchEffect to visualise effect of removing batch using a linear model
reads_forlim <- data.frame(reads_norm)
lim_batch <- covar$Strategy_Numeric
reads_lim <- removeBatchEffect(reads_forlim, batch=lim_batch)

pdf(file = "PCA_limma_logcorrect2_Strategy.pdf", width = 15, height = 12)
ggbiplot(rsem.pca, choices = 1:2, obs.scale = 1, var.scale = 1,groups = as.factor(covar$Strategy_Numeric), ellipse = TRUE, circle = TRUE, var.axes = F,labels=covar$Sample) +
  scale_color_discrete(name="Group", h=c(0,360) + 20) +
  theme(legend.direction = 'vertical', legend.position = 'right')
dev.off()

# Create the design matrix for the groups
design <- model.matrix( ~0+group+covar$Sex)
dim(design)

# Plot the estimated dispersions using BCV
y <- estimateGLMCommonDisp(y,design,verbose=T)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)
pdf(file="BCV.pdf")
plotBCV(y, xlab="Average log CPM", ylab="Biological coefficient of variation", pch=16, cex=0.2, col.common="red", col.trend="blue", col.tagwise="black")
dev.off()

# Make your contrasts (compare one group with the other)
Contrast <- glmLRT(fit, contrast=c(-1, 1, 0))

# Set or adjust Genes and Mutations according to their P-values 

summary(dt_Contrast <- decideTestsDGE(Contrast, p=1, adjust="BH"))
#UP/Down/Nosignificant genes
isDE_Contrast <- as.logical(dt_Contrast)
#give a vector of TRUE and FALSE
DE_Contrast <- rownames(y)[isDE_Contrast]
Contrast$table$FDR <- p.adjust(Contrast$table$PValue, method="BH")
#P-value adjustment (Given a set of p-values)
DEG_Contrast <- Contrast$table[Contrast$table$FDR <1,]
DEG_Contrast=DEG_Contrast[DEG_Contrast$FDR > 0,]
DEG_Contrast=DEG_Contrast[DEG_Contrast$PValue < 1,]

# P.Value histogram
pdf("Contrast_pvalues.pdf")
hist(Contrast$table$PValue, breaks=seq(0, 1, 0.05))
dev.off()

# Produce a relative log expression (RLE) plot
pdf("Contrast_lrt_RLE.pdf", 
width = 60, height = 15); par(cex = 3); par(mar = c(4,4,4,4))
plotRLE(Contrast$fitted.values, outline=FALSE, ylim=c(-6, 6), col=colors(), las=2)
dev.off()

# Build a biomaRt query
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="http://apr2022.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
##mouse <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="http://may2021.archive.ensembl.org", path="/biomart/martservice", dataset="mmusculus_gene_ensembl")
gene <- gsub("\\..*", "", rownames(DEG_Contrast)) 
results_DEG_Contrast <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id", "entrezgene_id", "gene_biotype","transcript_biotype","description", "chromosome_name", "start_position", "end_position", "strand","go_id"), filters="ensembl_transcript_id", values = gene, mart = mart, useCache = FALSE)
## results_DEG_Contrast <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end'), filters="ensembl_gene_id", values = gene, mart = mart, useCache = FALSE)

idx_DEG_Contrast=match(gsub("\\..*", "", rownames(DEG_Contrast)), (results_DEG_Contrast$ensembl_transcript_id)) 
DEG_Contrast$hgnc_symbol=results_DEG_Contrast$hgnc_symbol[idx_DEG_Contrast]
DEG_Contrast$ensembl_gene_id=results_DEG_Contrast$ensembl_gene_id[idx_DEG_Contrast]
DEG_Contrast$ensembl_transcript_id=results_DEG_Contrast$ensembl_transcript_id[idx_DEG_Contrast]
DEG_Contrast$entrezgene_id=results_DEG_Contrast$entrezgene_id [idx_DEG_Contrast]
DEG_Contrast$gene_biotype=results_DEG_Contrast$gene_biotype[idx_DEG_Contrast]
DEG_Contrast$transcript_biotype=results_DEG_Contrast$transcript_biotype[idx_DEG_Contrast]
DEG_Contrast$description=results_DEG_Contrast$description[idx_DEG_Contrast]
DEG_Contrast$chromosome_name=results_DEG_Contrast$chromosome_name[idx_DEG_Contrast]
DEG_Contrast$start_position=results_DEG_Contrast$start_position[idx_DEG_Contrast]
DEG_Contrast$end_position=results_DEG_Contrast$end_position[idx_DEG_Contrast]
DEG_Contrast$strand=results_DEG_Contrast$strand[idx_DEG_Contrast]
write.csv(DEG_Contrast, file="Contrast_DEG.csv")

top_Contrast=rownames(DEG_Contrast)
Contrast_cpm <- cpm(y)[top_Contrast, ]
write.csv(Contrast_cpm, file="Contrast_DEG_cpm.csv")

# Plot heatmap
unique_DEG_Contrast=unique(DEG_Contrast$entrezgene_id)
z_score_Contrast <- ((Contrast_cpm - rowMeans(Contrast_cpm))/apply(Contrast_cpm, 1, sd))
z_score_Contrast[z_score_Contrast > 3] = 2
z_score_Contrast[z_score_Contrast < -3] = -2
z_score_Contrast_gene <- z_score_Contrast[1:140,]
col_groups <- substr(colnames(z_score_Contrast), 1, 2)
z_score_Contrast[,col_groups == "1"] <- z_score_Contrast[,col_groups == "1"] * 5
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(z_score_Contrast)
mat_colors <- list(group = brewer.pal(2, "Pastel1"))
names(mat_colors$group) <- unique(col_groups)
pdf(file = "heatmap_plot11.pdf", width = 12, height = 6);
pheatmap(z_score_Contrast_gene, color= inferno(10), show_colnames = T,
annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = FALSE,
  fontsize          = 5,
  main              = "Default Heatmap", clustering_method = 'ward.D2'
)
dev.off()

# MA plot
pdf("Contrast_MA_plot.pdf")
plotSmear(Contrast, de.tags=DE_Contrast)
abline(h=c(-1,1), col="blue")
dev.off()

# GO Enrichment analysis
# BP
anno_bp <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "ensembl")
str(anno_bp)
allGenes_bp <- unique(unlist(anno_bp, use.names =TRUE))
str(allGenes_bp)
uniqueDEG_Contrast <- unique(results_DEG_Contrast$ensembl_gene_id)
head(uniqueDEG_Contrast)
geneList_Contrast <- factor(as.integer(allGenes_bp %in% uniqueDEG_Contrast))
names(geneList_Contrast) <- allGenes_bp
str(geneList_Contrast)
GOdata_bp_Contrast <- new("topGOdata", ontology = "BP", allGenes = geneList_Contrast, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl")
resultFisher_bp_Contrast <- runTest(GOdata_bp_Contrast, algorithm = "classic", statistic = "fisher")
resultKS_bp_Contrast <- runTest(GOdata_bp_Contrast, algorithm = "classic", statistic = "ks")
resultKS.elim_bp_Contrast <- runTest(GOdata_bp_Contrast, algorithm = "elim", statistic = "ks")
pdf("Contrast_gotree_BP.pdf")
showSigOfNodes(GOdata_bp_Contrast, score(resultFisher_bp_Contrast), firstSigNodes = 10, useInfo = "all")
dev.off()
allRes_BP_Contrast <- GenTable(GOdata_bp_Contrast, classicFisher = resultFisher_bp_Contrast, classicKS = resultKS_bp_Contrast, elimKS = resultKS.elim_bp_Contrast, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30)
for( i in 1:nrow(allRes_BP_Contrast)){ allRes_BP_Contrast$hgnc_symbol[i]=paste(as.vector(results_DEG_Contrast$hgnc_symbol)[which(grepl(allRes_BP_Contrast$GO.ID[i],results_DEG_Contrast$go_id))],collapse=",")}
for (i in 1:nrow(allRes_BP_Contrast)){ allRes_BP_Contrast$ensg[i]=paste(as.vector(results_DEG_Contrast$ensembl_gene_id)[which(grepl(allRes_BP_Contrast$GO.ID[i],results_DEG_Contrast$go_id))],collapse=",")}
write.csv(allRes_BP_Contrast, file="Contrast_toptree_BP.csv")

# CC
anno_cc <- annFUN.org("CC", mapping = "org.Hs.eg.db", ID = "ensembl")
allGenes_cc <- unique(unlist(anno_cc))
uniqueDEG_Contrast <- unique(results_DEG_Contrast$ensembl_gene_id)
geneList_Contrast <- factor(as.integer(allGenes_cc %in% uniqueDEG_Contrast))
names(geneList_Contrast) <- allGenes_cc
GOdata_cc_Contrast <- new("topGOdata", ontology = "CC", allGenes = geneList_Contrast, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl")
resultFisher_cc_Contrast <- runTest(GOdata_cc_Contrast, algorithm = "classic", statistic = "fisher")
resultKS_cc_Contrast <- runTest(GOdata_cc_Contrast, algorithm = "classic", statistic = "ks")
resultKS.elim_cc_Contrast <- runTest(GOdata_cc_Contrast, algorithm = "elim", statistic = "ks")
pdf("Contrast_gotree_CC.pdf")
showSigOfNodes(GOdata_cc_Contrast, score(resultFisher_cc_Contrast), firstSigNodes = 10, useInfo = "all")
dev.off()
allRes_cc_Contrast <- GenTable(GOdata_cc_Contrast, classicFisher = resultFisher_cc_Contrast, classicKS = resultKS_cc_Contrast, elimKS = resultKS.elim_cc_Contrast, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30)
for (i in 1:nrow(allRes_cc_Contrast)){ allRes_cc_Contrast$hgnc_symbol[i]=paste(as.vector(results_DEG_Contrast$hgnc_symbol)[which(grepl(allRes_cc_Contrast$GO.ID[i],results_DEG_Contrast$go_id))],collapse=",")}
for (i in 1:nrow(allRes_cc_Contrast)){ allRes_cc_Contrast$ensg[i]=paste(as.vector(results_DEG_Contrast$ensembl_gene_id)[which(grepl(allRes_cc_Contrast$GO.ID[i],results_DEG_Contrast$go_id))],collapse=",")}
write.csv(allRes_cc_Contrast, file="Contrast_toptree_CC.csv")

# MF
anno_mf <- annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "ensembl")
allGenes_mf <- unique(unlist(anno_mf))
uniqueDEG_Contrast <- unique(results_DEG_Contrast$ensembl_gene_id)
geneList_Contrast <- factor(as.integer(allGenes_mf %in% uniqueDEG_Contrast))
names(geneList_Contrast) <- allGenes_mf
GOdata_mf_Contrast <- new("topGOdata", ontology = "MF", allGenes = geneList_Contrast, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl")
resultFisher_mf_Contrast <- runTest(GOdata_mf_Contrast, algorithm = "classic", statistic = "fisher")
resultKS_mf_Contrast <- runTest(GOdata_mf_Contrast, algorithm = "classic", statistic = "ks")
resultKS.elim_mf_Contrast <- runTest(GOdata_mf_Contrast, algorithm = "elim", statistic = "ks")
pdf("Contrast_gotree_MF.pdf")
showSigOfNodes(GOdata_mf_Contrast, score(resultFisher_mf_Contrast), firstSigNodes = 10, useInfo = "all")
dev.off()
allRes_mf_Contrast <- GenTable(GOdata_mf_Contrast, classicFisher = resultFisher_mf_Contrast, classicKS = resultKS_mf_Contrast, elimKS = resultKS.elim_mf_Contrast, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30)
for (i in 1:nrow(allRes_mf_Contrast)){ allRes_mf_Contrast$hgnc_symbol[i]=paste(as.vector(results_DEG_Contrast$hgnc_symbol)[which(grepl(allRes_mf_Contrast$GO.ID[i],results_DEG_Contrast$go_id))],collapse=",")}
for (i in 1:nrow(allRes_mf_Contrast)){ allRes_mf_Contrast$ensg[i]=paste(as.vector(results_DEG_Contrast$ensembl_gene_id)[which(grepl(allRes_mf_Contrast$GO.ID[i],results_DEG_Contrast$go_id))],collapse=",")}
write.csv(allRes_mf_Contrast, file="Contrast_toptree_MF.csv")

# Unique genes from BP/CC/MF analysis
unique_DEG_Contrast=unique(DEG_Contrast$entrezgene_id)
go_DEG_Contrast=goana(unique_DEG_Contrast)
top_BP_DEG_Contrast=topGO(go_DEG_Contrast,ont="BP",n=30)
for (i in 1:nrow(top_BP_DEG_Contrast)){ top_BP_DEG_Contrast$hgnc_symbol[i]=paste(as.vector(results_DEG_Contrast$hgnc_symbol)[which(grepl(row.names(top_BP_DEG_Contrast)[i],results_DEG_Contrast$go_id))],collapse=",")}
for (i in 1:nrow(top_BP_DEG_Contrast)){ top_BP_DEG_Contrast$ensg[i]=paste(as.vector(results_DEG_Contrast$ensembl_gene_id)[which(grepl(row.names(top_BP_DEG_Contrast)[i],results_DEG_Contrast$go_id))],collapse=",")}
write.csv(top_BP_DEG_Contrast, file="Contrast_top1_BP.csv")
top_MF_DEG_Contrast=topGO(go_DEG_Contrast,ont="MF",n=30)
for (i in 1:nrow(top_MF_DEG_Contrast)){ top_MF_DEG_Contrast$hgnc_symbol[i]=paste(as.vector(results_DEG_Contrast$hgnc_symbol)[which(grepl(row.names(top_MF_DEG_Contrast)[i],results_DEG_Contrast$go_id))],collapse=",")}
for (i in 1:nrow(top_MF_DEG_Contrast)){ top_MF_DEG_Contrast$ensg[i]=paste(as.vector(results_DEG_Contrast$ensembl_gene_id)[which(grepl(row.names(top_MF_DEG_Contrast)[i],results_DEG_Contrast$go_id))],collapse=",")}
write.csv(top_MF_DEG_Contrast, file="Contrast_top1_MF.csv")
top_CC_DEG_Contrast=topGO(go_DEG_Contrast,ont="CC",n=30)
for (i in 1:nrow(top_CC_DEG_Contrast)){ top_CC_DEG_Contrast$hgnc_symbol[i]=paste(as.vector(results_DEG_Contrast$hgnc_symbol)[which(grepl(row.names(top_CC_DEG_Contrast)[i],results_DEG_Contrast$go_id))],collapse=",")}
for (i in 1:nrow(top_CC_DEG_Contrast)){ top_CC_DEG_Contrast$ensg[i]=paste(as.vector(results_DEG_Contrast$ensembl_gene_id)[which(grepl(row.names(top_CC_DEG_Contrast)[i],results_DEG_Contrast$go_id))],collapse=",")}
write.csv(top_CC_DEG_Contrast, file="Contrast_top1_CC.csv")

# KEGG Pathways
kg.hsa <- kegg.gsets(species = "hsa", id.type = "kegg", check.new=TRUE)
kegg.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
data(egSymb)
kegg.gs.sym <- lapply(kegg.gs, eg2sym)
lapply(kegg.gs.sym[1:3],head)
all_ensg_mart <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "entrezgene_id"), filters ="ensembl_transcript_id", values = sub("\\..*","",rownames(Contrast)), mart = mart, useCache = FALSE)
all_ensg_mart= all_ensg_mart[!is.na(all_ensg_mart$entrezgene_id),]
all_ensg_mart = all_ensg_mart[!duplicated(all_ensg_mart$entrezgene_id), ]
idx_Contrast=match(all_ensg_mart$ensembl_transcript_id, sub("\\..*","",rownames(Contrast$table)))
exp.fc=Contrast$table$logFC[idx_Contrast]
names(exp.fc)=all_ensg_mart$entrezgene_id [idx_Contrast]
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
head(fc.kegg.p $greater[, 1:5], 4)
head(fc.kegg.p $less[, 1:5], 4)
write.csv(rbind(fc.kegg.p$greater, fc.kegg.p$less), file ="Contrast_fc.kegg.p.csv")
fc.kegg.sig=sigGeneSet(fc.kegg.p)
write.csv(rbind(fc.kegg.sig, fc.kegg.sig), file="Contrast.kegg.sig.csv")
fc.kegg.p.2p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL, same.dir = F)
write.csv(rbind(fc.kegg.p.2p$greater, fc.kegg.p.2p$less), file="Contrast_fc.kegg.p.2p.csv")
fc.kegg.2p.sig=sigGeneSet(fc.kegg.p.2p)
write.csv(rbind(fc.kegg.2p.sig, fc.kegg.2p.sig), file="Contrast.kegg.2d.sig.csv")
```
