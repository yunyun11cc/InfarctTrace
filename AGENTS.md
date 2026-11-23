# InfarctTrace 工程约定

## 环境与入口

- 使用仓库根目录的 uv 项目和 `.venv`。依赖由 `pyproject.toml` 与 `uv.lock` 管理，按 `.python-version` 选择 Python。
- 从仓库根目录运行 `uv sync --locked` 和 `uv run --locked python experiments/<阶段>/main.py`。通过各阶段的 `--steps` 选择复跑范围。

## 实验与产物

- `experiments/` 保存分析实现、输入资产与实验文档。`datasets/` 保存源数据与数据说明，`docs/outputs.md` 是工程产物索引。
- `outputs/` 的一级目录与实际产生产物的实验阶段同名。阶段内按工作流归档，使用 `matrices/`、`tables/`、`figures/`、`reports/` 表达产物类型。
- 阶段 01 使用 `01_preprocessing_differential/`。阶段 02 使用 `02_feature_selection/<mode>/`，在对应模式下组织方法结果、重采样、EFS 聚合和生物学互证。阶段 03 使用 `03_external_validation/`，将队列清洗、表达与 ROC 评估、单细胞分析和跨队列汇总分别归档。
- `outputs/` 保存当前工程可追溯的矩阵、分析表、程序绘制的科研图和运行报告。组会展示、手工排版、重复图表以及已停用的产物通过 Git 历史查阅。实验说明和结果解读放在 `experiments/<阶段>/docs/`。
- 用阶段配置确定默认输出位置与跨阶段输入。移动产物时同步更新生产者、消费者、命令、索引和来源记录。临时复跑使用入口的输出覆盖参数或工作区外的临时目录。

## 验证

- 目录迁移保留分析表内容、基因和样本顺序、阈值、随机种子及冻结候选表校验值。核对迁移文件的哈希，并从标准入口验证新路径读写。
- 优先复用已保存的矩阵和统计表复核绘图及下游步骤。改变计算逻辑时用同一输入和环境复算对应结果。
