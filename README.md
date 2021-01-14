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


