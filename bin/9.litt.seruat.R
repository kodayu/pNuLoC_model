library(stringdist)
library(seqinr)
library(dplyr)
library(gridExtra)
#simmethod <- "osa"
library(Seurat, lib.loc = "/home/yukai6/R/x86_64-conda-linux-gnu-library/4.0")
allmethods <- c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex")
for (simmethod in allmethods) {
  ##tSNE
  pbmc <- readRDS(paste("../9.litingting/Seurat.", simmethod, ".rds", sep = ""))
  pbmc <- ScaleData(pbmc, do.scale = F, do.center = F)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 10000)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  ElbowPlot(pbmc, ndims = 50)
  pbmc <- FindNeighbors(pbmc, dims = 1:20)
  pbmc <- FindClusters(pbmc, resolution = 0.8)
  pbmc <- RunUMAP(pbmc, dims = 1:20)
  pbmc <- RunTSNE(pbmc, dims = 1:20, check_duplicates = FALSE)
  saveRDS(pbmc, file = paste("../9.litingting/Seurat.", simmethod, ".cluster.rds", sep = ""))
  p1 <- DimPlot(pbmc, reduction = "tsne",pt.size=0.5,group.by = "PredClass", 
                label = T, cols = c("gray", "green", "red"))
  p2 <- DimPlot(pbmc, reduction = "tsne",pt.size=0.5,group.by = "ValidClass", 
                label = T, cols = c("gray", "green", "red"))
  pdf(paste("../9.litingting/Seurat.", simmethod, ".cluster.pdf", sep = ""), width = 10, height = 4)
  grid.arrange(p1, p2, nrow = 1, ncol = 2)
  dev.off()
}

