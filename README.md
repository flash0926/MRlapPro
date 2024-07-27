# MRlapPro
- MRlapPro包是基于MRlap和GenomicSEM包的魔改版
1. 增加了读取文件的速度
2. 将部分代码修改成效率更高的替代版本
3. 修改了GenomicSEM中标准化数据时，无用反复读写文件带来的内存消耗，时间消耗
4. 将必要的文件直接封装到了包里
5. 删除了例子里的BMI和SBP 减少安装MRlapPro包的时间

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
```


## 引用

查看原MRlap的引用说明
https://github.com/n-mounier/MRlap/blob/master/README.md
 
也可以直接引用MRlapPro包
https://github.com/flash0926/MRlapPro
