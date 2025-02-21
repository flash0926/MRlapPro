# MRlapPro
- MRlapPro包是基于MRlap和GenomicSEM包的魔改版
1. 增加了读取文件的速度
2. 将部分代码修改成效率更高的替代版本
3. 修改了GenomicSEM中标准化数据时，无用反复读写文件带来的内存消耗，时间消耗
4. 将必要的文件直接封装到了包里
5. 删除了例子里的BMI和SBP 减少安装MRlapPro包的时间
6. 兼容更多的GWAS格式数据

- 耗时上，可以大概减少原MRlap包的90%无用耗时，主要目的是为了支持大量的文件批量分析。

- 修改或删除的源码部分已经在代码中标注

- MRlapPro包遵循GPL-2.0协议
- 原MRlap包：https://github.com/n-mounier/MRlap

## 安装

You can install the current version of `MRlap` with:

``` r
# Directly install the package from github
# install.packages("remotes")
# if use MRlap internal testing data
remotes::install_github("n-mounier/MRlap")
# install MRlapPro
remotes::install_github("flash0926/MRlapPro")
# install MRlapPro force update
remotes::install_github("flash0926/MRlapPro", force=T)
```

## 例子
``` r
BMI <- system.file("data/", "BMI_Data.tsv.gz", package="MRlap")
SBP <- system.file("data/", "SBP_Data.tsv.gz", package="MRlap")
A = MRlapPro::MRlap(exposure = BMI,
          exposure_name = "BMI_100Ksample",
          outcome = SBP,
          outcome_name = "SBP_100Ksample",
          ld = system.file("Data/eur_w_ld_chr", package="MRlapPro"),
          hm3 = system.file("Data/eur_w_ld_chr", "w_hm3.snplist", package="MRlapPro"))
A
```

## 特别说明
因为代码需要使用到IEU的clump功能
- 所以需要按照IEU提供的文档，进行Token配置，否则无法使用其功能。
- 也可以使用本地clump功能，查看包体参数说明传入plink以及bfile的文件。
可以下载网页中提供的LD.zip文件 解压到自己的工作目录使用
https://kimfiles.oss-cn-beijing.aliyuncs.com/LD.zip


还有一个问题
MRlap的流程中的mr（ivw）分析，snp的流转流程是

暴露筛阈值和结局的snp取交集后，再做的clump

而我们做MR的标准流程，应该是暴露筛阈值，先clump，再和结局的snp取交集

这里就会有一个问题，就是MRlap中的snp，暴露和结局先取交集，会导致snp减少，然后clump时会有较大的差异（snp会不同）

这里的差异，就会导致ivw的分析结果，和我们自己做的不一样

所以，如果纠结的同学，请自己修改mr部分的流程源码，在run_MR.R内


## 引用

查看原MRlap的引用说明
https://github.com/n-mounier/MRlap/blob/master/README.md
 
也可以直接引用MRlapPro包
https://github.com/flash0926/MRlapPro


## MendelR 搭配用法
可以跳转下方链接查看：
https://flash0926.yuque.com/org-wiki-flash0926-kivyu0/gpa1ys/ggcx9qoyhzsr1dgc

支持多个文件传入，根据步骤准备好环境，然后直接调用代码分析出结果