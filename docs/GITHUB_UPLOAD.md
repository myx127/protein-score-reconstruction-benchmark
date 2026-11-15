# 上传与版本归档

本代码包对应公开仓库 [myx127/protein-score-reconstruction-benchmark](https://github.com/myx127/protein-score-reconstruction-benchmark)，论文引用版本为 [v1.0.0](https://github.com/myx127/protein-score-reconstruction-benchmark/releases/tag/v1.0.0)。本目录就是仓库根目录。以下保留首次发布与后续更新的操作说明；现有仓库更新时直接提交、推送，不需再次创建同名仓库。

## 推荐上传顺序

1. 在本目录运行 `python -B tools/verify_release.py`，阅读 `validation/VALIDATION_SUMMARY.json` 和 `docs/DATA_ACCESS.md`。本包的可追溯清单是 `SHA256SUMS.json`。
2. 在 GitHub 创建空仓库，填写实际账户和仓库名称。因为本地已有 README、许可和忽略规则，不要再自动生成这些文件；创建流程见 [GitHub 官方说明](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository)。
3. 使用 GitHub Desktop 把本目录添加为本地仓库，或按下面的 Git 命令操作。提交前查看文件清单；应只有代码、元数据、汇总表、图件和文档。定量数据、本地环境、患者模型、访问凭据及真实数据运行输出不应出现。
4. 首次提交后推送到您的仓库，再在网页检查 README、文件清单和数据访问声明。创建 `v1.0.0` release，写明这是代码包版本，CRNP 算法版本为 0.2.0；见 [GitHub release 官方说明](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)。
5. 需要可引用 DOI 时，可在 Zenodo 连接 GitHub 并为该仓库启用归档，然后发布新的正式 release；启用步骤见 [Zenodo 官方指南](https://help.zenodo.org/docs/github/enable-repository/)。应以实际归档成功后显示的标识为准。
6. 将实际仓库 URL、release 标识和真实 DOI 回填 README、CITATION.cff 和稿件的 Code availability。内容修改后重建清单，保证提交内容与清单一致。

## 命令行方式

以下由您在仓库根目录执行。先在 GitHub 创建空仓库，并把下面的示例 URL 改成实际值。不要把密码或 token 写入脚本。

```powershell
git init -b main
git add .
git status --short
git diff --cached --stat
git commit -m "Prepare reproducible benchmark code v1.0.0"
git remote add origin https://github.com/YOUR_ACCOUNT/protein-score-reconstruction-benchmark.git
git push -u origin main
```

`YOUR_ACCOUNT` 只是示例占位，不是已经创建的仓库。`.gitignore` 在 Git 上传时生效；浏览器手动拖拽文件并不等于经过 Git 忽略规则。不要从受控运行目录直接拖拽结果进入公共仓库。

## 修改后如何重新核对

```powershell
python -B -m unittest discover -s tests -v
python -B tools/verify_release.py --write-report --write-manifest
python -B tools/verify_release.py
```

如果修改了计算入口或算法，还应运行相应的模拟/数值检查并更新验证摘要；修改 README 中的 URL 不需要重复完整的科学实验。`verify_release.py` 只做文件、字段、文本和哈希检查，不自动判断所有输出是否符合数据使用条件。

`.gitattributes` 保留文件字节，避免不同系统的 Git 自动换行转换破坏冻结 R 文件及清单哈希。生成独立压缩包时，可执行 `python -B tools/build_archive.py --output ../protein-score-reconstruction-benchmark-v1.0.0.zip`；输出必须位于仓库外且尚不存在，压缩包只收录通过清单核对的发布文件。

## Code availability 英文模板

> The analysis code, frozen feature/coverage metadata, aggregate results, and an independent synthetic-data demonstration are available at [actual repository URL], release [actual version/DOI]. HMACS participant-level data are not publicly distributed and may be requested from the corresponding author, Yan Zeng (zengyan68@wust.edu.cn), subject to applicable cohort and institutional approvals and data-use conditions. Third-party data remain subject to their original access conditions.

未真实公开前，不应把模板改写成“已公开”或提供虚构链接。代码开源与参与者数据受控访问可以并行成立。
