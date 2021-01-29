### muxiaoqin 20210125 Seurt的学习
setwd("C:/Wanglab/Mouxiaoqin/单细胞/result")
library(dplyr)
library(Seurat)
library(patchwork)
raw_data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
# 为17055个细胞中21027个基因的count值矩阵
dim(raw_data)
# [1] 21027 17055
# 创建seurat对象
# 基因至少在min.cells个细胞中表达
# 每个细胞中至少表达min.genes个基因
sc.object <- CreateSeuratObject(counts = raw_data, project = "macaca_brain", min.cells = 3, min.features = 200 )
dim(sc.object)
# [1] 17764 17055

#画出Count和feature数的分布图
FeatureScatter(object = sc.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sc.object[["percent.mt"]] <- PercentageFeatureSet(sc.object, pattern = "^ENSMFAG")
##sc.object[["percent.mt"]] <- PercentageFeatureSet(sc.object, pattern = "^MT-")
#人和鼠的线粒体基因才是"^MT-"
VlnPlot(sc.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

plot1 <- FeatureScatter(sc.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
## 因为是单细胞核测序，所有没有mt ×因为之前的步骤错了才没有！！！
#根据feature图删除离群值,选择feature大于200小于4500的细胞。我瞎选的。如何确定是离群值？？？看图自己选吧。
sc.object <- subset(sc.object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5)
dim(sc.object)
# [1] 17764 16292

#对数据进行归一化处理
sc.object <- NormalizeData(sc.object, normalization.method = "LogNormalize", scale.factor = 10000)
#分析差异表达的基因，即在数据集中找出在一些细胞中高表达一些细胞中低表达的基因,选取了方差最大的2000个feature做下游的分析
sc.object <- FindVariableFeatures(sc.object, selection.method = "vst", nfeatures = 2000)
#方差最大的十个基因
top10 <- head(VariableFeatures(sc.object), 10)
# 画出这些差异表达基因
plot1 <- VariableFeaturePlot(sc.object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot2
#scale data线性转化 
#scale就是将数据归一化为均值为0方差为1的数据
all.genes <- rownames(sc.object)
sc.object <- ScaleData(sc.object, features = all.genes)

###去批次
library(harmony)
sce@meta.data$group <- substring(sce@assays$RNA@counts@Dimnames[[2]],first=18,last=10000000L)
sce <- sce %>% 
  RunHarmony("group", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sce, 'harmony')

#pca降维
sc.object <- RunPCA(sc.object, features = VariableFeatures(object = sc.object))
print(sc.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc.object, dims = 1:2, reduction = "pca")
DimPlot(sc.object, reduction = "pca")

#如何决定数据集的维度如何选择PCs
#在选择几个PC进行下游的分析时，我们可以通过热图查看这些方差的来源
DimHeatmap(sc.object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sc.object, dims = 1:20, cells = 500, balanced = TRUE)
saveRDS(sc.object, file = "C:/Wanglab/Mouxiaoqin/单细胞/result/sc.object.rds")

#用JackStraw随机产生了一个数据集，又跑了一次PCA，构建了一个空集对照，发现显著的PC多集中在低P值的位置
nu.object <- JackStraw(sc.object, num.replicate = 100)
nu.object1 <- ScoreJackStraw(nu.object, dims = 1:20)
JackStrawPlot(nu.object1, dims = 1:20) + ggtitle("JackStraw")
#这需要很长的时间，可以用ElowPlot代替
library(ggplot2)
ElbowPlot(sc.object)+NoLegend() + ggtitle("ElbowPlot")
#用三个方法选择PCs
#1.观察PCs热图找到异质性的来源，这个热图怎么看
#2.用随机的对照数据集统计，耗时而且无法确定合适的cutoff
#3.ElowPlot直接选择PC的拐点
#选择不同的PCs会对下游的结果产生巨大的影响

#去批次后融合还不如没有去批次的，就不去了
#library(harmony)
#library(Rcpp)
#sce1 <- sc.object
#sce1@meta.data$group <- substring(sce@assays$RNA@counts@Dimnames[[2]],first=18,last=10000000L)
#sce <- sce %>% 
#  RunHarmony("group", plot_convergence = TRUE)
  #harmony_embeddings <- Embeddings(sce, 'harmony')
#options(repr.plot.height = 5, repr.plot.width = 12)
#p1 <- DimPlot(object = sce, reduction = "harmony", pt.size = .1, group.by = "group")
#p2 <- VlnPlot(object = sce, features = "harmony_1", group.by = "group",  pt.size = .1)
#p3 <- DimPlot(sce1, group.by = "group", reduction = "pca", pt.size = .1)


#细胞聚类分析
#Seurat采用基于图的聚类方法，简而言之这些方法将单元格嵌入图结构中，
#例如K近邻（KNN）图，在具有相似特征表达模式的单元格之间绘制边，
#然后尝试将该图划分为高度互连的“准斜体”或“社区”。

#首先基于PCA空间中的欧式距离构建一个KNN图
#然后细化他们的边缘权重（相似性
#运用了一个模块优化方法Louvain algorithm (default) or SLM
#选了15个PCs
sc.object <- readRDS("C:/Wanglab/Mouxiaoqin/单细胞/result/sc.object.rds")
sc.object1 <- sc.object
sc.object1 <- FindNeighbors(sc.object1, reduction='pca', dims = 1:6)
sc.object1 <- FindClusters(sc.object1, resolution = 0.5)
#可以看一下前5个细胞属于哪一类
head(Idents(sc.object12), 5)

#非线性降维
#例如t-SNE和umap这样的非线性降维方法会考虑数据的各个方面
#将他们降维到一个低维空间里
#上面基于图像的聚类方法应该和这些低维图共定位
sc.object1 <- RunUMAP(sc.object1, dims = 1:9)
p1 <-DimPlot(sc.object1, reduction = "umap")
p2 <-DimPlot(sc.object1, reduction = "umap")
p1
p2
dev.off()
saveRDS(sc.object1, file = "C:/Wanglab/Mouxiaoqin/单细胞/result/sce1.rds")
#找到差异表达基因
#对每个cluster定义一个阳性marker和阴性marker
#一个cluster和所有细胞对比找出差异表达基因
#比如cluster1的marker基因，min.pct参数设置两组细胞中被检测的feature所占的比例
#cluster1.markers <- FindMarkers(sc.object, ident.1 = 1, min.pct = 0.25)
#head(cluster1.markers, n = 5)
#找到cluster5中能与cluster0到3区分开的基因
#cluster5.markers <- FindMarkers(sc.object, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)
# 找到每个cluster中相比其他剩余所有cluster的marker基因，且只保留positive
sc.object1 <- readRDS("C:/Wanglab/Mouxiaoqin/单细胞/result/sc.object.rds")
library(dplyr)
library(Seurat)
library(patchwork)
sc.object.markers <- FindAllMarkers(sc.object1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TEST <- sc.object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#用一些参数对差异表达基因进行test，比如ROC test回归每个marker的统计力，从0到1的数值

###细胞类型注释
if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("SingleR")
browseVignettes("SingleR")
library(SingleR)
library(tibble)
hpca.se <- HumanPrimaryCellAtlasData()##运行了好久
bpe.se <- BlueprintEncodeData()
seurat.obj <- sc.object1
seurat.obj@meta.data$cell.type <- Idents(seurat.obj)
test <- as.SingleCellExperiment(seurat.obj)
Anno <- SingleR(test = test,
                ref = list(HP = hpca.se , BP = bpe.se),
                labels = list(hpca.se$label.main , bpe.se$label.main),
                method = "cluster",
                cluster = test$cell.type)

Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj <- RenameIdents(seurat.obj, new.cluster.ids)
DimPlot(seurat.obj, reduction = "umap")


############没有什么用啊
#用软件注释
library(scCATCH)
sce_markers <- findmarkergenes(sc.object1,species = 'Human',cluster = 'All',match_CellMatch = TRUE,cancer = NULL,tissue = 'Brain',cell_min_pct = 0.25,logfc = 0.25,pvalue = 0.05)
clu_ann2 <- scCATCH(sce_markers$clu_markers,species = 'Human',cancer = NULL,tissue = 'Brain')
write.table(clu_ann, file = "C:/Wanglab/Mouxiaoqin/单细胞/result/clu_ann.txt", sep="\t")

                       
###手动注释
####画叠在一起的vlnplot
library(patchwork)
library(ggplot2)
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
StackedVlnPlot(sc.object1, features = finally_marker[18:19], pt.size = 0)+ scale_x_continuous(breaks = 0:17,labels=1:18)

ppt_marker <- c('SLC17A7','GAD1','MAG','A2M','SLC1A2','SATB2','GAD2','TAC1','PCDH8','DRD2','ADORA2A','PENK','PCP4',
                'NECAB2','LMO7','CALB1','PDGFRA','CSPG4','GJA1','MBP','MOBP','MOG','RELN','AIF1','CX3CR1','PTPRC',
                'HLA−DRA','TIAM1','CADPS2','CACNA2D1','PHLDB2','KCNH8','PDE1A','GRB14','ND6','PVALB','VIP','SST','LAMP5','COL9A1','SCL1A2','STAB1')
markers_exist_test <- intersect(x=sc.object1@assays[["RNA"]]@counts@Dimnames[[1]], y = ppt_marker)

finally_marker<-c('SLC17A7','CACNA2D1','SATB2','GAD2','PVALB','CADPS2','MBP',
'MOBP','MOG','MAG','PCP4','SLC1A2','PHLDB2','KCNH8','PDE1A','ND6','PDGFRA','COL9A1','A2M')
VlnPlot(sc.object, features = markers_exist, pt.size = 0)
FeaturePlot(sc.object1, features = markers1)

final.object <- RenameIdents(object = sce1, 
                             "0" = "EXC1",
                             "1" = "EXC2",
                             "2" = "EXC3",
                             "3" = "EXC4",
                             "4" = "EXC5",
                             "5" = "EXC4",
                             "6" = "EXC7",
                             "7" = "Int",
                             "8" = "EXC8",
                             "9" = "OLIG",
                             "10" = "Purkinje cells",
                             "11" = "Astro",
                             "12" = "EXC9",
                             "13" = "EXC10",
                             "14" = "EXC11",
                             "15" = "EXC12",
                             "16" = "OPC",
                             "17" = "Endo")

DimPlot(final.object, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
#用VlnPlot展示这个基因在所有cluster中的表达概率分布
#featurePlot可以在tSNE或这PCA上标记feature
#Doheatmap可以画出指定基因的表达热图
top5 <- sc.object.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(sc.object, features = top5$gene)
top10 <- sc.object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(top10, file = "./top10marker", sep = "\t")

#指定细胞类型去识别cluster
#需要凭借背景知识判断marker之后属于哪类细胞
#常用的注释类型包有singleR，数据库有CellMarker
#markers根据已知文献查询得到脑组织常见细胞类型的makers
#为了方便改名，另用一个final.object记载改名后的结果
# 改成缩写，或者在注释出现问题时方便修改

#####不太理解这个clustree
library(clustree)
library(ggraph)
intergrated=FindClusters(sc.object1,resolution = c(seq(0,1,.2)))
clustree(intergrated@meta.data,prefix = "RNA_snn_res.") 

intergrated <- AddMetaData(intergrated,intergrated@reductions$umap@cell.embeddings,
                           col.name = c("UMAP_1","UMAP_2"))
intergrated <- AddMetaData(intergrated,intergrated@reductions$pca@cell.embeddings,
                           col.name = colnames(intergrated@reductions$pca@cell.embeddings))
clustree_overlay(intergrated, prefix = "RNA_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2")
