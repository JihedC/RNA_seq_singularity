library(Seurat)
library(readr)
library(dplyr)
library(patchwork)

out <- read.csv("MORC3_heterozygous.csv", header = TRUE, row.names = 1)
head(out)

mydata <- CreateSeuratObject(counts =  t(out), project = "Morc3_heterozygous")

mydata
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^mt-")
head(mydata@meta.data)

pdf(file = "/exports/humgen/jihed/scTE/QC_Morc3_heterozygous.pdf")
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000)

mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mydata), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "/exports/humgen/jihed/scTE/Variable_feature_heterozygous.pdf")

plot1 + plot2

dev.off()

all.genes <- rownames(mydata)
mydata <- ScaleData(mydata, features = all.genes)

mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))

pdf(file = "/exports/humgen/jihed/scTE/elbow_plot_heterozygous.pdf")
ElbowPlot(mydata)
dev.off()

mydata <- RunUMAP(mydata, dims = 1:20)
pdf(file = "/exports/humgen/jihed/scTE/UMAP_heterozygous.pdf")
DimPlot(mydata, reduction = "umap")
dev.off()

hetero.markers <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(hetero.markers, file="/exports/humgen/jihed/scTE/markers_hetero.txt", sep="\t", quote=F)

saveRDS(mydat, file="/exports/humgen/jihed/scTE/heterozygous.rds")