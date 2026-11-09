# GSE252027 教学笔记

[阶段 01 数据与差异分析.ipynb](01_gse252027_data_and_differential.ipynb) 按“样本纳入 → 长度缩放计数 → 文库大小因子 → 总体时间检验 → 逐时间点比较 → D28 敏感性 → GO/KEGG/GSEA → 六张图”的顺序讲解，边读输出边对照 `experiments/01_preprocessing_differential/` 的代码。先阅读 [GSE252027 数据手册](../datasets/GSE252027/README.md)，再按 [阶段 01 运行说明](../experiments/01_preprocessing_differential/README.md) 生成输出。

笔记只读取已有阶段 01 输出，不重新拟合模型，按顺序执行即可跑完。样本取舍与统计口径见 [阶段方案](../experiments/01_preprocessing_differential/protocol.md)，公式与真实算例见 [数学原理](../docs/math/01_preprocessing_differential.md)。

使用项目 `.venv` 内核。正式工程入口在 `experiments/`。

在项目根目录启动 JupyterLab。

```powershell
uv run --locked jupyter lab
```
