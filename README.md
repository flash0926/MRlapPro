# MRlapPro
 A MRlapPro package
 使用的MRlapPro包是基于MRlap和GenomicSEM包的魔改版
1. 增加了读取文件的速度
2. 将部分代码修改成效率更高的替代版本
3. 修改了GenomicSEM无用反复读写文件带来的内存消耗，时间消耗
4. 将必要的文件直接封装到了包里
    耗时上，可以大概减少原MRlap包的90%无用耗时，主要目的是为了支持大量的文件批量分析。
修改或删除的源码部分已经在代码中标注

MRlapPro包遵循GPL-2.0协议
原MRlap包：https://github.com/n-mounier/MRlap

## Installation

You can install the current version of `MRlap` with:

``` r
# Directly install the package from github
# install.packages("remotes")
remotes::install_github("flash0926/MRlapPro")
```

## 例子
`` r
BMI <- system.file("data/", "BMI_Data.tsv.gz", package="MRlap")
SBP <- system.file("data/", "SBP_Data.tsv.gz", package="MRlap")
A = MRlapPro::MRlap(exposure = BMI,
          exposure_name = "BMI_100Ksample",
          outcome = SBP,
          outcome_name = "SBP_100Ksample",
          ld = system.file("Data/eur_w_ld_chr", package="MRlapPro"),
          hm3 = system.file("Data/eur_w_ld_chr", "w_hm3.snplist", package="MRlap"))
``


## 引用

查看原MRlap的引用说明
https://github.com/n-mounier/MRlap/blob/master/README.md
 
也可以直接引用MRlapPro包
https://github.com/flash0926/MRlapPro
