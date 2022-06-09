# 染色质开放区域识别 Open Chromatin Regions Dection

文件说明：

```
.
├── README.md                   本文件
├── build
│   └── EMPTY                   build结果        
├── cmd
│   ├── parser                  从PaGenBase中提取高表达基因的TSS区域
│   ├── pick_none_tss.py        低表达基因的TSS区域的仿真生成脚本
│   ├── plot_acc.py             绘制Accurency的分类条形图
│   ├── plot_features.py        绘制曲线特征的核密度图
│   └── subsample.py            对cfDNA测序数据进行随机采样
├── data
│   ├── GRCh37.gene.bed         全基因组TSS的位置信息
│   └── HK.all.txt.bed          从PaGenBase中提取持家基因TSS区域
├── entry.sh                    特征提取+模型训练的运行脚本
├── entry_fractions.sh          引入不同覆盖度下的特征提取+模型训练的运行脚本
├── extract_features.py         特征提取主程序
├── predict.py                  模型预测主程序
└── utils
    ├── Config.py               配置文件
    ├── Peak.py                 波峰类
    ├── draw_utils.py           绘图库
    ├── feat_utils.py           特征提取库
    └── log.py                  日志库
```

<details>
  <summary>Click to expand!</summary>
  
```
github.com/AlDanial/cloc v 1.90  T=0.03 s (604.1 files/s, 66522.1 lines/s)
-------------------------------------------------------------------------------
Language                     files          blank        comment           code
-------------------------------------------------------------------------------
Python                          13            282            212           1261
Bourne Shell                     3             18              5             65
Markdown                         1              2              0             27
-------------------------------------------------------------------------------
SUM:                            17            302            217           1353
-------------------------------------------------------------------------------
```
</details>