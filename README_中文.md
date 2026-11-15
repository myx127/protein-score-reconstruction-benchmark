# 本稿分析代码：本地使用与公开发布

本目录已整理为可独立上传 GitHub 的代码仓库，代码包版本为 **1.0.0**。对应论文为 *Benchmarking Protein Network Score Reconstruction under Structured Feature Missingness*。其中 CRNP 的算法版本仍为 **0.2.0**。

**公开内容包括代码、独立模拟示例、蛋白与覆盖元数据、论文汇总结果和图件；不包含 HMACS 或 Emory 参与者级定量数据、真实参与者评分或用真实数据拟合的模型。**

## 先运行什么

请在本目录打开终端。推荐 Python 3.11，另需 R；具体版本见 [环境安装说明](docs/INSTALL.md)。无需获得 HMACS 数据即可完成以下运行。

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m pip install --no-deps -e .
Rscript R/install_dependencies.R
python -B -m unittest discover -s tests -v
python -B examples/run_demo.py --bootstrap 50
```

`Rscript` 不在 PATH 时，先设置 `$env:RSCRIPT` 为实际 Rscript.exe 路径；已有依赖库可通过 `R_LIBS` 指定。必须使用本仓库的 CRNP R 实现，程序不会改用 Python 近似替代。

默认输出为 `outputs/synthetic_demo/`，含两个虚拟批次、三个虚拟模块、四类掩码的 **24 个恢复评价单元**，以及 **35 名虚拟参与者、六个分支、18 个单元**的配对检测示例。所有数值和编号独立生成。默认 50 次重采样用于演示流程，不代表论文的正式推断，也不应引用示例效应量。

可选运行：

```powershell
python -B examples/run_demo.py --output outputs/demo_with_timing --bootstrap 50 --with-timing
python -B scripts/check_r_contracts.py
python -B scripts/redraw_paper_figures.py
python -B tools/verify_release.py
```

第一条额外运行 24 个工作负载、72 次预热、360 条正式计时；不同机器耗时会不同。绘图入口只读论文汇总表，不需要参与者数据。每次分析请选择新的输出目录，以免覆盖已有完整运行。

## 获批数据后如何运行

1. 向通讯作者曾燕老师 **zengyan68@wust.edu.cn** 申请 HMACS 数据。具体审批和使用条件由队列及所在机构决定。
2. 将获批定量矩阵和本地配置放在本目录之外。恢复分析需要两份冻结标准化批次矩阵；Emory 分支另外需要 HMACS 的源定义矩阵及允许使用的六份检测输入。
3. 根据 [输入格式](docs/INPUT_SCHEMAS.md) 修改 `config/` 示例的外部副本，再执行分析入口。

```powershell
python -B scripts/run_benchmark.py recovery --config ../governed_inputs/hmacs.local.json --private-workdir ../private_runs/recovery_run1
python -B scripts/run_benchmark.py timing --config ../governed_inputs/hmacs.local.json --private-workdir ../private_runs/timing_run1
python -B scripts/run_benchmark.py emory --config ../governed_inputs/emory.local.json --private-workdir ../private_runs/emory_run1
```

路径相对于各 JSON 配置文件所在目录解析。论文元数据默认从代码包读取，无需移动到受控目录。**定量输入与所有真实数据运行输出必须在公开目录外**；程序会检查这一边界。配置示例未绑定您的电脑路径。

## 已保留的科学规则

- 完整及直接评分按相应蛋白权重绝对值之和归一化，残差含分母变化的影响。
- 来源残差岭回归仅用来源批次训练；保留原调参网格、单标准误规则和截距处理。
- CRNP 采用冻结 R 代码，恢复分析输出通过来源完整评分的均值与标准差逆变换；不利用目标完整评分校准。
- Within-CF 使用目标训练折的完整评分，是较强信息条件下的诊断，不能写成部署模型或理论上界。
- 先平均各方向的 36 个单元，再等权平均两个方向；相对改善与信息捕获比例在汇总风险后计算。
- HMACS 正式推断为固定模型、预处理及预测后的 2,000 次目标参与者配对 bootstrap；Emory 为 10,000 次。
- HMACS 的原始直接评分与 Emory 的来源校准、分支标准化直接评分不同，解释时应保留这一差别。

详细依据见 [方法](docs/METHODS.md) 和 [复现范围与验证](docs/REPRODUCIBILITY.md)。本包从冻结预处理矩阵开始，不重新完成原始质谱处理或联合模块发现。本轮采用独立模拟数据验证入口和数值契约，**没有重新拟合真实 HMACS/Emory 数据来声称正式结果全部复现**。

## 上传 GitHub

可以将本目录作为仓库根目录，按 [上传指南](docs/GITHUB_UPLOAD.md) 操作。本轮未创建公共仓库、未上传文件。`.gitignore` 已排除运行目录、环境和本地配置；它不能替代上传前检查。

`validation/` 保存本轮实际运行的脱敏验证摘要。`SHA256SUMS.json` 用于核对发布文件。仓库创建后，再把真实 URL 加入 README、CITATION.cff 和论文的 Code availability；DOI 应在真实归档成功后填写。
