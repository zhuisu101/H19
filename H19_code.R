library(Seurat)
library(dplyr)
setwd('/media/sdb/ljt/GSE213740_RAW')
AD1 = Read10X(data.dir ="AD1") 
AD1 <- CreateSeuratObject(counts = AD1, project = "AD1",min.cells = 3, min.features = 200) 

folders=list.files('./', pattern='[123456]$')
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder)
})

scRNA <- merge(sceList[[1]], 
               y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                     sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]]), 
               add.cell.ids = folders)
scRNA

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

scRNA <- subset(scRNA, 
                subset = nCount_RNA > 200 &   
                  nFeature_RNA < 5000 &       
                  nFeature_RNA > 600 &            
                  percent.mt < 30)       

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

class(scRNA@meta.data) 
a = scRNA@meta.data
table(scRNA@meta.data$orig.ident)
a$group = c(rep('AD',61777),rep('N',35965))
scRNA = AddMetaData(object = scRNA,metadata=a)
table(scRNA@meta.data$group)

scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA), 10); top10 

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scRNA); plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE); plot2

scRNA <- ScaleData(scRNA, features = rownames(scRNA))              

scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE) 
scRNA <- JackStraw(scRNA, num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA, dims = 1:20) 
JackStrawPlot(scRNA, dims = 1:20)

ElbowPlot(scRNA)

library(clustree)
scRNA <- FindNeighbors(scRNA, dims = 1:15) #dims = 1:15 
scRNA <- FindClusters(object = scRNA,
                      resolution = c(seq(.1,1.6,.2)))
clustree(scRNA@meta.data, prefix = "RNA_snn_res.")  

scRNA <- RunUMAP(scRNA, dims = 1:15) 

scRNA <- RunTSNE(scRNA, dims = 1:15) 

scRNA <- FindNeighbors(scRNA, dims = 1:15) #dims = 1:15 
scRNA <- FindClusters(scRNA, resolution = 0.1) 

scRNA2 <- FindClusters(scRNA, resolution = 0.3) 

head(scRNA2@meta.data)

library(cowplot)
library(stringr)
scRNA2$patient=str_replace(scRNA$orig.ident,"_.*$","")
p1 <- DimPlot(scRNA2, reduction = "tsne", group.by = "patient", pt.size=0.5)
p2 <- DimPlot(scRNA2, reduction = "tsne", group.by = "ident",   pt.size=0.5, label = TRUE,label.size = 5, repel = TRUE);p2 
p3 <- DimPlot(scRNA2, reduction = "umap", group.by = "patient", pt.size=0.5);p3
p4 <- DimPlot(scRNA2, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE);p4

scRNA.markers <- FindAllMarkers(scRNA2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scRNA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

top5 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

####ComplexHeatmap####
library(ComplexHeatmap)
expression_data = GetAssayData(scRNA,slot = 'counts') 
expression_data = log2(expression_data+1)
expression_data[1:4,1:4]

cluster_info = sort(scRNA2$seurat_clusters)

expression_data <- as.matrix(expression_data[top5$gene, names(cluster_info)])

library(devtools)
devtools::install_github("caleblareau/BuenColors")

library(BuenColors)
col <- jdb_color_maps[1:21] # 要与cluster_info的level一致
names(col) <- levels(cluster_info)

Heatmap(expression_data,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_split = cluster_info)

#cell annotation
VlnPlot(scRNA, features = c("CD3D", "CD3E","GZMA","NKG7"))   #T cell   
VlnPlot(scRNA, features = c('CD19', 'CD79A', 'CD79B', 'MS4A1'))  # B cell  
VlnPlot(scRNA, features = c('CST3', 'LYZ', 'FCGR3B', 'CSF3R'))  # Neutrophil  
VlnPlot(scRNA, features = c('FGF7', 'MME',"DCN","LUM","FBN1","COL1A1",'COL3A1'))  # Fibroblasts  
VlnPlot(scRNA, features = c('CD68','CD163', 'S100A8', 'S100A9'))  # Monocytes 
VlnPlot(scRNA, features = c('TAGLN', 'MYH11', 'CNN1','ACTG2','ACTA2')) #SMC 
VlnPlot(scRNA, features = c("VWF","CLDN5","ICAM2","CDH5")) #Endothelial
VlnPlot(scRNA, features = c("CCR7","IFITM1")) #DC  

marker_genes = c("CD3D", "CD3E","GZMA",'NKG7',                     
                 'CD19', 'CD79A', 'CD79B', 'MS4A1',                   
                 'CST3', 'LYZ', 'FCGR3B', 'CSF3R',          
                 'FGF7', 'MME',"DCN","LUM","FBN1","COL1A1",'COL3A1',             
                 'CD68','CD163', 'S100A8', 'S100A9',              
                 'TAGLN', 'MYH11', 'CNN1','ACTG2','ACTA2',                       
                 "VWF","CLDN5","ICAM2","CDH5",
                 "CCR7","IFITM1")                        

# 1. DotPlot
DotPlot(scRNA2,features = marker_genes)+coord_flip()
# 2. VlnPlot
library(MySeuratWrappers) 
VlnPlot(scRNA, features = marker_genes,stacked=T,pt.size =0.25) 

new.cluster.ids <- c("0"="SMC",
                     "1"="Monocytes", 
                     "2"="Fibroblasts", 
                     "3"= "Monocytes", 
                     "4"= "T cell", 
                     "5"= "Endothelial",
                     "6"= "SMC", 
                     "7"= "Fibroblasts", 
                     "8"= "SMC",
                     '9'='Neutrophil',
                     '10'='Fibroblasts',
                     '11'='Endothelial',
                     '12'='Endothelial',
                     '13'='Endothelial',
                     '14'='Fibroblasts',
                     '15'='Monocytes',
                     '16'='Neutrophil',
                     '17'='Monocytes',
                     '18'='SMC',
                     '19'='B cell',
                     '20'='SMC') 
names(new.cluster.ids) <- levels(scRNA2$seurat_clusters)

scRNA2 <- RenameIdents(scRNA2, new.cluster.ids)

scRNA2[['cell_type']] = unname(new.cluster.ids[scRNA2@meta.data$seurat_clusters])

VSMC = scRNA2[, Idents(scRNA2) %in% c( "SMC" )]
df.data = VSMC@assays$RNA@counts %>% as.data.frame()
head(rownames(df.data))

#####substract H19+ & H19- cells####
H19 <- df.data[c("H19"),] %>% as.data.frame() %>% t() %>% as.data.frame()
H19$var <- ifelse(H19$H19 > 0, "H19+", "H19-")
H19$id<-rownames(H19)

cell <- as.data.frame(VSMC@meta.data)
cell$id<-rownames(cell)
cell<-cell[,c("id","orig.ident")]

count<-merge(H19,cell,by="id")
table(count$orig.ident)
VSMC@meta.data$H19 <- count$var 
table(VSMC@meta.data$H19)

DefaultAssay(VSMC) <- "RNA"
all.genes <- rownames(VSMC)
VSMC <- ScaleData(VSMC, features = all.genes)

Idents(VSMC) <- "H19" 
VSMC@active.assay
head(VSMC@active.ident)

#
DEG_AD_Con <- FindMarkers(VSMC, min.pct = 0.25, logfc.threshold = 0.06,
                          group.by = "H19",
                          ident.1 ="H19+",
                          ident.2="H19-")
write.csv(DEG_AD_Con,file = "DEG_smc_AD_Con.CSV")

library(clusterProfiler)
library(GEOquery)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

DEG <- read.csv("DEG_smc_AD_Con.CSV")
head(DEG)

OrgDb = "org.Hs.eg.db"
gene_convert <- bitr(DEG$X, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = OrgDb)
DEG = DEG%>% inner_join(gene_convert,by=c("X"="SYMBOL"))
head(gene_convert)
head(DEG)

ont = "BP"# posiible value: BP, CC, MF, all
go.results <- enrichGO(DEG$ENTREZID, keyType="ENTREZID",ont="BP",OrgD = OrgDb, readable = FALSE)
head(go.results)
go.results <- enrichGO(DEG$ENTREZID, keyType="ENTREZID",ont="BP",OrgD = OrgDb, readable = TRUE)
head(go.results)
write.csv(go.results,file = "go.results.CSV")

go----------------------------
p1 <-dotplot(go.results,label_format=10000)
p1

p2 <-barplot(go.results, showCategory = 10,label_format=10000)
p2

p3 <-emapplot(pairwise_termsim(go.results))
p3

organism = "hsa"

#kegg 
kegg.results <- enrichKEGG(DEG$ENTREZID, organism = organism)
head(kegg.results)

kegg.results <- setReadable(kegg.results, OrgDb = OrgDb, keyType='ENTREZID')
head(kegg.results)
write.csv(kegg.results,file = "kegg.results.CSV")

#
P4 <-dotplot(kegg.results,label_format=10000)
P4

p5 <-barplot(kegg.results, showCategory = 10,label_format=10000)
p5

p6<-emapplot(pairwise_termsim(kegg.results))
p6


####scMetabolism####
library(scMetabolism)
countexp.Seurat<-sc.metabolism.Seurat(obj = VSMC, method = "VISION", 
                                      imputation = F, ncores = 2, metabolism.type = "KEGG")
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
metabolism.matrix

#DimPlot
DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Glycolysis / Gluconeogenesis", 
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, size = 1)

#Dotplot
input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:30]

DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "group", 
                   norm = "y")

input.pathway1 = c('Glycolysis / Gluconeogenesis','Citrate cycle (TCA cycle)','Pyruvate metabolism','Oxidative phosphorylation',
                   'Fatty acid elongation','Fatty acid biosynthesis','Fatty acid degradation','Synthesis and degradation of ketone bodies',
                   'Glycerolipid metabolism','Galactose metabolism','Glycerolipid metabolism','Arachidonic acid metabolism',
                   'Butanoate metabolism','Ascorbate and aldarate metabolism','Arachidonic acid metabolism')

input.pathway2 = c('Alanine, aspartate and glutamate metabolism','D-Glutamine and D-glutamate metabolism','Arginine biosynthesis',
                   'Arginine and proline metabolism')
DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway2, 
                   phenotype = "group", 
                   norm = "y")


cairo_pdf('scMetabolism_VSMC.pdf',width=7,height =5.5,family="Arial")
DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway1, 
                   phenotype = "group", 
                   norm = "y")
dev.off()