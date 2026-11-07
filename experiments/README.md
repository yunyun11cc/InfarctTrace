# Experiments

每个目录对应一个明确的实验问题，编号按开展顺序递增。

| 实验 | 问题 | 状态 |
|---|---|---|
| [01_fibroblast_screening](01_fibroblast_screening/README.md) | 心梗后 D7 的成纤维细胞表达变化能否在独立队列复现？ | 已完成探索性基线，含中文教程与已执行 notebook |

实验目录保存 `README.md`、`protocol.md`、`run.py`，有教学需求时增加 `tutorial.md` 和 `notebook.ipynb`。设计写明输入、比较、阈值和解释边界；结果输出到项目根目录 `outputs/<实验目录名>/`。后续实验直接引用 `datasets/` 下的源数据，不各自复制下载文件。

当前代码仍集中在具体实验中。只有多个实验实际共用某段计算代码时，再提取公共模块。
