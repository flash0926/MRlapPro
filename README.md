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