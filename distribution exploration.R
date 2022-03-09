#accession is GSE84465
setwd("E:/Bioinformatics/data")
a <- read.table("./rawdata/GSE84465_GBM_All_data.csv.gz")
head(rownames(a))
tail(rownames(a),10)
a <- a[1:(nrow(a)-5),]#remove the last five rows

b <- read.table("./SraRunTable.txt",
                sep = ",", header = T)
table(b$Patient_ID) # 4 human primary GBM samples
table(b$TISSUE) # 2 cell types: tumor cores and peripheral regions
table(b$TISSUE,b$Patient_ID)

# tumor and peripheral information of groups
head(colnames(a))
head(b$plate_id)
head(b$Well)
b.group <- b[,c("plate_id","Well","TISSUE","Patient_ID")]
paste0("X",b.group$plate_id[1],".",b.group$Well[1])
b.group$sample <- paste0("X",b.group$plate_id,".",b.group$Well)
head(b.group)
identical(colnames(a),b.group$sample)

# select tumor cell
index <- which(b.group$TISSUE=="Tumor")
length(index)#2343
group <- b.group[index,] 
head(group)

alt <- a[,index] 
dim(alt)# 23460  2343
identical(colnames(alt),group$sample)

#qsmooth normalization
matrix <- alt[,c(1,2)]
factor <- as.factor(alt[1,])
factor
library(qsmooth)
alt_qs <- qsmooth(object = alt,group_factor = factor)get
qs.data <- alt_qs@qsmoothData
sessionInfo()

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)

save(qs.data,file = "./daqs.data.csv")
qs.data1 <- CreateSeuratObject(counts = qs.data)
qs.data2 <- FindVariableFeatures(qs.data1,
                            selection.method = "vst", nfeatures = 2000)

qs3 <- ScaleData(qs.data2)
qs4 <- RunPCA(object = qs3, pc.genes = VariableFeatures(sce))

DimHeatmap(qs4, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(qs4)

qs4 <- FindNeighbors(qs4, dims = 1:18)
qs4 <- FindClusters(qs4, resolution = 0.9)
table(qs4@meta.data$RNA_snn_res.0.9)

#tSNE
set.seed(123)
qs4 <- RunTSNE(object = qs4, dims = 1:18, do.fast = TRUE)
DimPlot(qs4,reduction = "tsne",label=T)
LabelClusters(DimPlot(qs4, reduction = "tsne"),id = 'ident')

VlnPlot(qs4, features = c("MS4A1", "CD79A"))#is this plot okay?
info <- qs4@meta.data

info <- read.table("./info0.csv")

#check distribution
qs4.0 <- info[info$V6 == "0",]#cluster 0
qs4.0 <- qs4.0$V1
a.0 <- a[,qs4.0]
a.0.avg <- data.frame(gene = rownames(a.0), means = rowMeans(a.0[,-1]))
a.0.avg <- a.0.avg$means
a.0.avg <- a.0.avg[a.0.avg != 0]
a.0.avg.log <- log2(a.0.avg)

qs4.1 <- info[info$V6 == "1",]#cluster 1
qs4.1 <- qs4.1$V1
a.1 <- a[,qs4.1]
a.1.avg <- data.frame(gene = rownames(a.1), means = rowMeans(a.1[,-1]))
a.1.avg <- a.1.avg$means
a.1.avg <- a.1.avg[a.1.avg != 0]
a.1.avg.log <- log2(a.1.avg)

qs4.2 <- info[info$V6 == "2",]#cluster 2
qs4.2 <- qs4.2$V1
a.2 <- a[,qs4.2]
a.2.avg <- data.frame(gene = rownames(a.2), means = rowMeans(a.2[,-1]))
a.2.avg <- a.2.avg$means
a.2.avg <- a.2.avg[a.2.avg != 0]
a.2.avg.log <- log2(a.2.avg)

qs4.3 <- info[info$V6 == "3",]#cluster 3
qs4.3 <- qs4.3$V1
a.3 <- a[,qs4.3]
a.3.avg <- data.frame(gene = rownames(a.3), means = rowMeans(a.3[,-1]))
a.3.avg <- a.3.avg$means
a.3.avg <- a.3.avg[a.3.avg != 0]
a.3.avg.log <- log2(a.3.avg)

qs4.4 <- info[info$V6 == "4",]#cluster 4
qs4.4 <- qs4.4$V1
a.4 <- a[,qs4.4]
a.4.avg <- data.frame(gene = rownames(a.4), means = rowMeans(a.4[,-1]))
a.4.avg <- a.4.avg$means
a.4.avg <- a.4.avg[a.4.avg != 0]
a.4.avg.log <- log2(a.4.avg)

qs4.5 <- info[info$V6 == "5",]#cluster 5
qs4.5 <- qs4.5$V1
a.5 <- a[,qs4.5]
a.5.avg <- data.frame(gene = rownames(a.5), means = rowMeans(a.5[,-1]))
a.5.avg <- a.5.avg$means
a.5.avg <- a.5.avg[a.5.avg != 0]
a.5.avg.log <- log2(a.5.avg)

qs4.6 <- info[info$V6 == "6",]#cluster 6
qs4.6 <- qs4.6$V1
a.6 <- a[,qs4.6]
a.6.avg <- data.frame(gene = rownames(a.6), means = rowMeans(a.6[,-1]))
a.6.avg <- a.6.avg$means
a.6.avg <- a.6.avg[a.6.avg != 0]
a.6.avg.log <- log2(a.6.avg)

qs4.7 <- info[info$V6 == "7",]#cluster 7
qs4.7 <- qs4.7$V1
a.7 <- a[,qs4.7]
a.7.avg <- data.frame(gene = rownames(a.7), means = rowMeans(a.7[,-1]))
a.7.avg <- a.7.avg$means
a.7.avg <- a.7.avg[a.7.avg != 0]
a.7.avg.log <- log2(a.7.avg)

qs4.8 <- info[info$V6 == "8",]#cluster 8
qs4.8 <- qs4.8$V1
a.8 <- a[,qs4.8]
a.8.avg <- data.frame(gene = rownames(a.8), means = rowMeans(a.8[,-1]))
a.8.avg <- a.8.avg$means
a.8.avg <- a.8.avg[a.8.avg != 0]
a.8.avg.log <- log2(a.8.avg)

qs4.9 <- info[info$V6 == "9",]#cluster 9
qs4.9 <- qs4.9$V1
a.9 <- a[,qs4.9]
a.9.avg <- data.frame(gene = rownames(a.9), means = rowMeans(a.9[,-1]))
a.9.avg <- a.9.avg$means
a.9.avg <- a.9.avg[a.9.avg != 0]
a.9.avg.log <- log2(a.9.avg)

qs4.10 <- info[info$V6 == "10",]#cluster 10
qs4.10 <- qs4.10$V1
a.10 <- a[,qs4.10]
a.10.avg <- data.frame(gene = rownames(a.10), means = rowMeans(a.10[,-1]))
a.10.avg <- a.10.avg$means
a.10.avg <- a.10.avg[a.10.avg != 0]
a.10.avg.log <- log2(a.10.avg)

qs4.11 <- info[info$V6 == "11",]#cluster 11
qs4.11 <- qs4.11$V1
a.11 <- a[,qs4.11]
a.11.avg <- data.frame(gene = rownames(a.11), means = rowMeans(a.11[,-1]))
a.11.avg <- a.11.avg$means
a.11.avg <- a.11.avg[a.11.avg != 0]
a.11.avg.log <- log2(a.11.avg)

qs4.12 <- info[info$V6 == "12",]#cluster 12
qs4.12 <- qs4.12$V1
a.12 <- a[,qs4.12]
a.12.avg <- data.frame(gene = rownames(a.12), means = rowMeans(a.12[,-1]))
a.12.avg <- a.12.avg$means
a.12.avg <- a.12.avg[a.12.avg != 0]
a.12.avg.log <- log2(a.12.avg)

qs4.13 <- info[info$V6 == "13",]#cluster 13
qs4.13 <- qs4.13$V1
a.13 <- a[,qs4.13]
a.13.avg <- data.frame(gene = rownames(a.13), means = rowMeans(a.13[,-1]))
a.13.avg <- a.13.avg$means
a.13.avg <- a.13.avg[a.13.avg != 0]
a.13.avg.log <- log2(a.13.avg)

qs4.14 <- info[info$V6 == "14",]#cluster 14
qs4.14 <- qs4.14$V1
a.14 <- a[,qs4.14]
a.14.avg <- data.frame(gene = rownames(a.14), means = rowMeans(a.14[,-1]))
a.14.avg <- a.14.avg$means
a.14.avg <- a.14.avg[a.14.avg != 0]
a.14.avg.log <- log2(a.14.avg)

qs4.15 <- info[info$V6 == "15",]#cluster 15
qs4.15 <- qs4.15$V1
a.15 <- a[,qs4.15]
a.15.avg <- data.frame(gene = rownames(a.15), means = rowMeans(a.15[,-1]))
a.15.avg <- a.15.avg$means
a.15.avg <- a.15.avg[a.15.avg != 0]
a.15.avg.log <- log2(a.15.avg)

qs4.16 <- info[info$V6 == "16",]#cluster 16
qs4.16 <- qs4.16$V1
a.16 <- a[,qs4.16]
a.16.avg <- data.frame(gene = rownames(a.16), means = rowMeans(a.16[,-1]))
a.16.avg <- a.16.avg$means
a.16.avg <- a.16.avg[a.16.avg != 0]
a.16.avg.log <- log2(a.16.avg)

qs4.17 <- info[info$V6 == "17",]#cluster 17
qs4.17 <- qs4.17$V1
a.17 <- a[,qs4.17]
a.17.avg <- data.frame(gene = rownames(a.17), means = rowMeans(a.17[,-1]))
a.17.avg <- a.17.avg$means
a.17.avg <- a.17.avg[a.17.avg != 0]
a.17.avg.log <- log2(a.17.avg)

boxplot(a.0.avg.log,	a.1.avg.log,	a.2.avg.log,	a.3.avg.log,	a.4.avg.log,	a.5.avg.log,	a.6.avg.log,	a.7.avg.log,	a.8.avg.log,	a.9.avg.log,	a.10.avg.log,	a.11.avg.log,	a.12.avg.log,	a.13.avg.log,	a.14.avg.log,	a.15.avg.log,	a.16.avg.log,	a.17.avg.log,
        xlab = "cell clusters", ylab = "log2(expression level)", frame.plot = F, col = "darkgrey", border="blue")
axis(side = 1, at = c(1:18), labels = c(1:18))
ggplot(a.0.avg.log,	a.1.avg.log,	a.2.avg.log,	a.3.avg.log,	a.4.avg.log,	a.5.avg.log,	a.6.avg.log,	a.7.avg.log,	a.8.avg.log,	a.9.avg.log,	a.10.avg.log,	a.11.avg.log,	a.12.avg.log,	a.13.avg.log,	a.14.avg.log,	a.15.avg.log,	a.16.avg.log,	a.17.avg.log
       +geom_boxplot())
