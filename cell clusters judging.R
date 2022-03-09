#accession is GSE84465
a <- read.table("C:/Users/Jekyll/Desktop/R/GSE84465_GBM_All_data.csv")
#rownames are symbol IDs
#colnames are samples
head(rownames(a))
tail(rownames(a),10)
a <- a[1:(nrow(a)-5),]#remove the last five rows

#primitive counts data
#3,589 cells of 4 human primary GBM samples, accession number GSE84465
#2,343 cells from tumor cores and 1,246 cells from peripheral regions
b <- read.table("C:/Users/Jekyll/Desktop/R/SraRunTable.txt",
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

a.filt <- a[,index] 
dim(a.filt)# 23460  2343
identical(colnames(a.filt),group$sample)

#qsmooth normalization
matrix <- a.filt[,c(1,2)]
factor <- as.factor(a.filt[1,])
factor
library(qsmooth)
a.filt_qs <- qsmooth(object = a.filt,group_factor = factor)
qs.data <- a.filt_qs@qsmoothData
sessionInfo()

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)

save(qs.data,file = "C:/Users/Jekyll/Desktop/R/qs.data.csv")
qs.data1 <- CreateSeuratObject(counts = qs.data)
qs.data2 <- FindVariableFeatures(qs.data1,
                            selection.method = "vst", nfeatures = 2000)
# 中心化，为下一步PCA做准备
qs3 <- ScaleData(qs.data2)
qs4 <- RunPCA(object = qs3, pc.genes = VariableFeatures(sce))
#PC_ 1 
#Positive:  C1QB, SAT1, CD74, HLA-DRB1, HLA-DRA, SLC11A1, RGS1, IFI30, FTL, GPR183,OLR1, CD83, VSIG4, ACSL1, CD14, ALOX5AP, PLEK, CTSB, APOE, STAB1,ZNF331, FCGBP, NFKBID, TFRC, ITGAX, C3, CCL3, PLIN2, RANBP2, CCL4 
#Negative:  CLU, CNN3, PTPRZ1, GPM6B, TUBA1A, CPE, BAALC, PON2, DDR1, MAP1B,MT3, F3, GPM6A, PTN, GFAP, TSPAN7, GATM, NRCAM, TRIM9, ITGB8 
#SYT11, MAGED1, CKB, ITM2C, BCAN, GPRC5B, PMP2, GAP43, FHL1, AQP4 
#PC_ 2 
#Positive:  PLP1, CLDND1, CNDP1, CNP, PLLP, ENPP2, QDPR, MOG, SPOCK3, KLK6 UGT8, STMN1, EDIL3, OMG, MAL, C11orf9, CARNS1, ERMN, CNTN2, TMEM144 ,ANLN, C8orf46, FXYD6, MOBP, OPALIN, SEPT4, SH3GL3, DNAJC6, TF, SCD5 
#Negative:  MT2A, IGFBP7, C1R, CHI3L1, CAV1, ALDOA, C8orf4, FLNA, SPARC, GAPDH, SERPING1, SHC1, IGFBP2, SOD2, CA12, ACTB, CSRP2, SPOCD1, CHPF, CYR61, CXCL14, TNC, JUNB, C1S, TAGLN, ANXA1, WWTR1, OSMR, MGST1, ENO2 
#PC_ 3 
#Positive:  CLDND1, PLP1, CNDP1, ENPP2, SEPT4, KLK6, MOG, QDPR, EDIL3, PTGDS, SPOCK3, CARNS1, MAL, ERMN, C11orf9, CNTN2, TMEM144, SLC44A1, MOBP, OPALIN, APLP1, TF, MBP, EFHD1, ANLN, NID1, UGT8, CAPN3, SH3GL3, PPAP2C 
#Negative:  TSPAN31, GPM6A, CDK4, TSFM, GNS, TBK1, SRGAP1, CTDSP2, DYRK2, CCND2, PMP2, MEST, LHFPL3, PTPRZ1, MEG3, TRIM9, SLC1A2, CLU, DTNA, CNR1, MDM2, METTL1, FXYD6, AVIL, GFAP, MAP2, BCAN, DBI, SCARA3, BAALC 
#PC_ 4 
#Positive:  NID1, TSPAN31, CDK4, TSFM, TBK1, IFITM1, MEG3, DYRK2, MYO1B, FRZB, UACA, TIMP3, GNS, PDLIM1, LAMA4, COL6A1, GGT5, CD248, MDM2, ISLR, LHFPL3, IGFBP4, FN1, CTDSP2, DCN, COL18A1, TSC22D1, CYP1B1, FBLN2, SRGAP1 
#Negative:  PON2, NDRG2, SLC1A3, TTYH1, FAM107A, DDR1, LRIG1, MLC1, AQP4, CST3, GPM6B, GATM, PAQR4, GPRC5B, METRN, PLEKHB1, LPL, FADS2, HOPX, PLP1, MT3, AIF1L, ID3, AHCYL1, ADCYAP1R1, NCAN, GPR37L1, TRIM47, MBP, CXCL14 
#PC_ 5 
#Positive:  TSPAN31, GNS, CTDSP2, SRGAP1, TBK1, MEG3, TSFM, CNR1, DYRK2, GLUL, GAP43, MAP1B, AVIL, CDK4, MDM2, SCARA3, CLU, LHFPL3, METTL1, METTL21B, BAALC, CRYAB, TUBB3, APOL4, IFI6, GFAP, PMP2, GALNTL2, BCL6, CSF1 
#Negative:  DBI, ATP5A1, HMGN2, HES6, C1orf61, CD9, KCNE1L, TUBA1B, MSMO1, AQP4, TUSC3, PDLIM1, PTN, HIST1H4C, BCAN, PRCP, S100B, HSPB1, RARRES2, CA2, AGT, CD34, DDIT3, PSAT1, TM4SF1, KPNA2, NREP, LIMA1, STMN1, TIMP4 

#前12个主成分
DimHeatmap(qs4, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(qs4)
#判断最终选取的主成分数，此处判断18个
qs4 <- FindNeighbors(qs4, dims = 1:18)
qs4 <- FindClusters(qs4, resolution = 0.9)
#Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#Number of nodes: 2343
#Number of edges: 71361
#Maximum modularity in 10 random starts: 0.8583
#Number of communities: 18
#Elapsed time: 0 seconds
table(qs4@meta.data$RNA_snn_res.0.9)
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
#390 371 284 169 141 127 121 109 104  83  82  78  70  53  45  43  42  31 
#分成18个clusters

#tSNE可视化
set.seed(123)
qs4 <- RunTSNE(object = qs4, dims = 1:18, do.fast = TRUE)
DimPlot(qs4,reduction = "tsne",label=T)

#寻找每个cluster的高变代表基因，并选取前5个，进行可视化
table(qs4@meta.data$seurat_clusters) 
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
#390 371 284 169 141 127 121 109 104  83  82  78  70  53  45  43  42  31
p <- list()
for( i in unique(qs4@meta.data$seurat_clusters) ){
  markers_df <- FindMarkers(object = qs4, ident.1 = i, min.pct = 0.25)
  print(x = head(markers_df))
  markers_genes =  rownames(head(x = markers_df, n = 4))
  p1 <- VlnPlot(object = qs4, features =markers_genes,log =T,ncol = 2)
  p[[i]][[1]] <- p1
  p2 <- FeaturePlot(object = qs4, features=markers_genes,ncol = 2)
  p[[i]][[2]] <- p2
}
###############p_val avg_log2FC pct.1 pct.2    p_val_adj
#CYBA   6.184804e-76       -Inf 0.996     1 1.450955e-71
#LAPTM5 3.865581e-75       -Inf 0.996     1 9.068652e-71
#SRGN   5.824499e-75       -Inf 0.996     1 1.366427e-70
#TYROBP 1.917306e-74       -Inf 0.996     1 4.498001e-70
#FCER1G 1.883614e-72       -Inf 0.996     1 4.418958e-68
#BASP1  3.873381e-72       -Inf 0.996     1 9.086953e-68
#SOX10  4.636147e-25        Inf     1     1 1.087640e-20
#PCSK1N 1.051096e-21  -151.0312     1     1 2.465870e-17
#C1QL1  2.027363e-21   525.0891     1     1 4.756194e-17
#FGF12  1.486649e-19        Inf     1     1 3.487678e-15
#NKAIN4 1.892244e-19        Inf     1     1 4.439204e-15
#LUZP2  4.637767e-19       -Inf     1     1 1.088020e-14

p[[1]][[2]]
p[[1]][[1]]


