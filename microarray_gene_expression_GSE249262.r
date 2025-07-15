# This script analyzes the microarray gene expression dataset GSE249262.
# It performs data preprocessing, normalization, and dimensionality reduction using UMAP.
################################################################
# Install necessary libraries
package_list <- c("umap","BiocManager")
for (pkg in package_list) {
  if (!require(pkg, quietly = T))
    install.packages(pkg, dependencies = T)
    library(pkg, quietly = TRUE)
}

bioconductor_packages <- c("GEOquery","limma")
for (packg in bioconductor_packages) {
  if (!require(packg, quietly = T))
    BiocManager::install(packg, force = TRUE)
    library(packg, quietly = TRUE)
}

##### Differential expression analysis with limma

# load series and platform data from GEO
gset <- getGEO("GSE249262", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL23159", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- make.names(gset@phenoData@data[["timepoint:ch1"]])
sel <- which(gsms !="cell.line.control")
gsms <- gsms[sel]

# filter out excludes samples (marked as "cell line control")
gset <- gset[ ,sel]

# assign samples to groups and set up design matrix
gsms <- factor(gsms)
levels(gsms) <- make.names(c("Tumor-Base","Control","Tumor-Week 10","Tumor-Week 4"))
gset$group <- gsms
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gsms)

# skip missing values
gset <- gset[complete.cases(exprs(gset)), ] 

# log2 transformation
ex <- exprs(gset)
exprs(gset) <- log2(ex)

# calculate precision weights and show plot of mean-variance trend by groups
v <- voomaByGroup(gset, group=gset$group, design, plot=T, cex=0.2, pch=".", col=1:nlevels(gset$group))

# attach gene annotation
v$genes <- fData(gset)

# fit linear model
fit <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
gsms_levels <- levels(gsms)
cts <- c(paste(gsms_levels[2],"-",gsms_levels[1],sep=""), paste(gsms_levels[2],"-",gsms_levels[4],sep=""), paste(gsms_levels[2],"-",gsms_levels[3],sep=""))
contrast.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, sort.by="B", adjust.method = "BH", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","SPOT_ID","SPOT_ID.1"))
write.table(tT, file="top_250_DE_gene.tsv", row.names=F, sep="\t")

