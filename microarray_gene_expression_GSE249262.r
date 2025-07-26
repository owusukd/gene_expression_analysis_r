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
groups <- make.names(c("Tumor-Base","Control","Tumor-Week 10","Tumor-Week 4"))
levels(gsms) <- groups
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

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, sort.by="B", adjust.method = "BH", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="BH", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change) for each contrast
par(mfrow=c(2,2), mar=c(4,4,2,1))
for (ct in 1:ncol(fit2$coefficients)) {
  volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20, highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
}

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
par(mfrow=c(2,2), mar=c(4,4,2,1))
for (ct in 1:ncol(fit2$coefficients)) {
  plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
  abline(h=0)
}

# General expression data analysis
# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gsms)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE249262", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gsms[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE249262", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gsms, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gsms, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gsms), pch=20, col=1:nlevels(gsms), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)


