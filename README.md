## 单细胞seurat方法用R学习篇
### 数据说明
> 数据：`crab-eating macaque single nuclei 10X V3 `食蟹猴单核测序 这个是single nuclei 测序 可以对premRNA进行定量

> 单核和单细胞的区别？？？

> 数据位置：
```
data_VISp="/picb/neurosys/LJ/BigData/GLHu/20191218/"
data_VISp2="/picb/neurosys/LJ/BigData/GLHu/20191225/"
```
> 不是把所有日期的数据放到一起，而是每个单独跑，最后去批次。

> 食蟹猴参考基因组数据： 

`/data/neurosys-svr2/Junjiem/DecisionTree/monkey_selfsequencing/ref/macaca_fascicularis` 

> 此外（也可以在[ensemble](http://ensemblgenomes.org/)里面找）怎么找？？？

### 学习资源
> [seurat](https://satijalab.org/seurat/)用于后续分析
> [cellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)下机数据处理的官方软件
> [周周笔记](https://github.com/small-west/single_cell_RNA_seq/blob/main/10xgenomics_data_preprocesing.md)里是鼠的数据做的单细胞
> [微信单细胞测序流程文章](https://mp.weixin.qq.com/s/jXxoRHC1FcHQMGbgPCADoA)
> 先学cellRanger，run的时候抽时间看文献，那三篇文献就是些review。这三篇文献如何链接进来？？？
####multiqc的安装
```
conda config --add channels bioconda
onda config --add channels conda-forge
conda create -n py3.5 python=3.5
conda env list
conda activate py3.5
conda install multiqc
multiqc -help
multiqc ./
```
参照：[multiqc的安装](https://www.jianshu.com/p/4783ffbb1347)

###前处理--cellranger
####制作参考基因组
```
(base) -bash-4.2$ 
nohup wget ftp://ftp.ensembl.org/pub/release-102/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_5.0.dna_sm.toplevel.fa.gz &
nohup wget ftp://ftp.ensembl.org/pub/release-102/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_5.0.102.chr.gtf.gz &
gunzip Macaca_fascicularis.Macaca_fascicularis_5.0.dna_sm.toplevel.fa.gz
gunzip Macaca_fascicularis.Macaca_fascicularis_5.0.102.chr.gtf.gz
## 为什么是用dna_sm这个，与dna的有什么区别？
cellranger mkgtf Macaca_fascicularis.Macaca_fascicularis_5.0.102.chr.pre_mrna.gtf Macaca_fascicularis.Macaca_fascicularis_5.0.102.chr.pre_mrna.filtered.gtf --attribute=gene_biotype:protein_coding
cellranger mkref --genome=Macaca_fascicularis.Macaca_fascicularis_5.0.102.chr.pre_mrna.filtered --fasta=Macaca_fascicularis.Macaca_fascicularis_5.0.dna_sm.toplevel.fa --genes=Macaca_fascicularis.Macaca_fascicularis_5.0.102.chr.pre_mrna.filtered.gtf
```
####count
```
refer_path=/picb/neurosys/mouxiaoqin/BigData/ref/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_5.0.102.chr.pre_mrna.filtered
##1225的样本
fastq_path1225=/picb/neurosys/mouxiaoqin/BigData/GLHu/20191225/
cellranger count --id=19R139 \
--transcriptome=$refer_path \
--fastqs=$fastq_path1225 \
--sample=19R139-2-1,19R139-2-2,19R139-2-3,19R139-2-4

##1218的样本
fastq_path1218=/picb/neurosys/mouxiaoqin/BigData/GLHu/20191218/
cellranger count --id=19X118 \
--transcriptome=$refer_path \
--fastqs=$fastq_path1218 \
--sample=19X118-2-1,19X118-2-2,19X118-2-3,19X118-2-4
```
####aggr
```


```
