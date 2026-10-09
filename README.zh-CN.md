# RepBase 转 RepeatMasker 适配器

[English](README.md) · **简体中文**

将 RepBase EMBL 记录转换为 RepeatMasker 兼容的 FASTA 重复序列库，提供可复现的分类映射、序列规范化、重复处理，以及可追溯的转换报告。

**测试环境：** RepBase 31.09 · RepeatMasker 4.2.4 · RMBlast 2.17.1

> 本项目是独立的转换工具，**并非** GIRI 官方的 RepBase RepeatMasker Edition 发行版。本仓库**不**分发任何 RepBase 序列数据。

## 概述

近期的 RepBase 发行版以 EMBL 格式提供重复序列及其注释。这些记录必须先转换为 RepeatMasker 兼容的 FASTA 头部格式，才能通过 `-lib` 选项使用。

本工具提供：

- 将 RepBase EMBL `.ref` 记录转换为 `>sequence_id#class/subclass` 形式的 FASTA。
- 使用改编自 [Terrier](https://github.com/rbturnbull/terrier) 项目的 RepBase→RepeatMasker 映射表进行分类。
- 大小写不敏感的关键词匹配；无法判定分类时保留 `Unknown`。
- 重复标识符的解析，并留下可追溯的重命名记录。
- 非 IUPAC 序列字符的规范化。
- 生成本地转换与审计报告。

原始的 RepBase 文件不会被修改。

## 环境要求

- **Python ≥ 3.11**（使用标准库 `tomllib`；已在 3.13.13 测试）。无需任何第三方 Python 包。
- 在相应 GIRI 许可下获得的 RepBase EMBL 记录。
- [RepeatMasker](https://www.repeatmasker.org/)（已在 4.2.4 测试）。
- [RMBlast](https://www.repeatmasker.org/RMBlast.html)，用于重复注释（已在 2.17.1 测试）。
- BLAST+ 的 `makeblastdb`，用于可选的数据库索引校验。

无需 Terrier 神经网络推理，也无需安装 REPET。转换器只使用**分类映射表**，不使用 Terrier 的预测模型。

## 使用方法

### 1. 获取 RepBase EMBL 数据

从官方发布渠道下载并解压所需的 RepBase 版本。请将解压后的源文件放在本仓库之外的目录中。

本文测试的转换使用 RepBase 31.09 的 39 个主 `.ref` 文件，**不包含** `appendix/` 子目录（其中大量记录 ID 与主文件重复）。转换器只读取顶层 `*.ref` 文件。

### 2. 转换数据库

```bash
python repbase_to_rm.py \
    --input-dir /path/to/RepBase31.09.embl \
    --mapping   mappings/repbase-to-repeatmasker.toml \
    --output    output/repbase_31.09.rm.fa \
    --out-dir   output/reports
```

审计报告（`id_mapping.tsv`、`unmapped_records.tsv`、`base_normalization.tsv`、`duplicate_report.tsv`）会写到输出库旁边，或在使用 `--out-dir` 时写入该目录。全部选项见 `python repbase_to_rm.py --help`。

### 3. 校验 FASTA 库

```bash
makeblastdb -in output/repbase_31.09.rm.fa -dbtype nucl -parse_seqids
```

构建完成的数据库应具有唯一的序列标识符和合法的核苷酸序列。

可选：将分类结果与参考库交叉验证（例如官方 2018 RepeatMasker Edition）：

```bash
python scripts/cross_validate.py \
    --old-lib /path/to/repbase_2018.lib \
    --new-lib output/repbase_31.09.rm.fa
```

### 4. 运行 RepeatMasker

```bash
RepeatMasker -e rmblast \
    -lib output/repbase_31.09.rm.fa \
    -xsmall -gff \
    -dir rm_results \
    genome.fa
```

转换后的库通过 `-lib` 显式指定，不会替换 RepeatMasker 已安装的参考库。用户可根据物种、研究问题和重复类别对结果注释进行过滤。

### 5. 测试

合成数据测试（不需要 RepBase 序列，仅用标准库）：

```bash
python tests/test_conversion.py
# 或：python -m unittest discover -s tests -v
```

## 转换策略

**分类。** 取每条记录 `KW` 字段中第一个大小写不敏感匹配的关键词，映射到 RepeatMasker 分类体系。没有识别到分类的条目保留为 `Unknown`，不做猜测。

**重复标识符。** 共享同一标识符且规范化后序列相同的记录会被合并；共享同一标识符但序列不同的记录会获得唯一后缀（`_dup2`、`_dup3`……）。重复比较在序列规范化之后进行。

**序列规范化。** 非 IUPAC 的 `X` 与 `O` 转换为 `N`，`U` 转换为 `T`；合法的 IUPAC 简并字符原样保留。

**标识符净化。** RepeatMasker FASTA 标识符中非法的字符（`:`、`(`、`)`、`@`）替换为 `_`。

**可追溯性。** 转换会在本地保留报告，关联原始标识符、输出标识符、分类判定与序列规范化事件。

## 验证

适配器于 2026-10-09 使用 RepBase 31.09 与 RepeatMasker 4.2.4 完成评估。

| 指标 | 结果 |
| --- | --- |
| 输入 EMBL 记录 | 126,885 |
| 输出 FASTA 记录 | 124,930 |
| 已分类记录 | 124,067 |
| 未分类记录 | 863（0.69%） |
| 分类覆盖率 | 99.31% |
| 与 2018 RepeatMasker Edition 的顶层类一致率（共有名称） | 98.10% |
| `makeblastdb -parse_seqids` | 通过 |
| RepeatMasker 功能测试 | 通过 |
| 与官方 FASTA 的不同序列（规范化后） | 124,924 = 124,924（集合一致） |
| 真实基因组比较（NTN 甲虫，chr2 前 10 Mb） | 31.09 = 18.08% 掩蔽，2018 = 14.66%；+3.42pp，主要落在已分类 TE |

功能测试确认了 Gypsy、L1、Helitron、Penelope 及一条未分类重复的注释与分类输出；随机片段对照未产生任何重复命中。

这些测试确立了转换兼容性与基本注释功能，**不**代表全基因组注释的灵敏度或特异度，也不代表优于精选的 2018 RepeatMasker Edition。

完整验证细节（含真实基因组比较）见 [`docs/VALIDATION.md`](docs/VALIDATION.md)。

## 局限

- 输出是自定义库，不是官方精选的 RepBase RepeatMasker Edition。
- 分类基于显式关键词映射，而非基于序列的预测。
- 2018 RepeatMasker Edition 中的部分亚类区分在本库中被归到更粗的超家族层级（如 `DNA/TcMar-Tc1` → `DNA/TcMar`）。
- 完整转换库不会自动限制到特定分类群。
- 转换库包含非 TE 的重复类别（卫星、结构 RNA）；下游应用应自行选取合适的类别。
- 有 863 条输出记录（0.69%）因 `KW` 中无可映射的超家族而保持 `Unknown`；它们被保留而非丢弃。

## 数据可用性与许可

RepBase 由遗传信息研究所（GIRI）维护，受其数据使用与再分发条款约束。GIRI 的学术用户条款限制对数据库、其组成部分及衍生材料的再分发。

> **2026 年更新 —— RepBase 将转为 CC0。** 2026 年 6 月，GIRI 宣布 RepBase 将以 **CC0 公有领域许可**发布，作为将 RepBase 与 Dfam 统一为单一开放框架的努力的一部分（Storer 等，*Mobile DNA* 17:17，2026；DOI [10.1186/s13100-026-00409-9](https://doi.org/10.1186/s13100-026-00409-9)）。这是一项**分多次发行推进的过渡**：该公告不意味着当前 Dfam 4.0 已完整收录 RepBase 31.09，也不意味着现有的各个 RepBase 发行版今天即可再分发。在过渡覆盖你所使用的具体发行版与日期之前，请仍按适用于该版本的条款从 GIRI 获取 RepBase。本仓库依然不分发任何 RepBase 数据。

本仓库包含转换软件与文档，但**不**分发 RepBase 源记录、转换后的共识序列或源自序列的测试数据。用户须自行在适用许可下获取 RepBase。

分类映射改编自 Terrier 项目（Apache-2.0）。见 [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)。

本仓库的原创代码以 MIT 许可发布，见 [`LICENSE`](LICENSE)。该许可不延伸到 RepBase 数据。

## 参考文献与致谢

- **RepBase：** Genetic Information Research Institute — <https://www.girinst.org/>
- **RepeatMasker：** <https://www.repeatmasker.org/>
- **Terrier：** <https://github.com/rbturnbull/terrier>
- **Terrier 论文：** Turnbull et al. (2025). *Terrier: a deep learning repeat classifier*. Briefings in Bioinformatics, 26(4), bbaf442. <https://doi.org/10.1093/bib/bbaf442>
- **RepBase–Dfam 统一（2026）：** Storer, J. M., Hubley, R. M., Rosen, J. B., Wheeler, T. J., & Smit, A. F. A. *Unifying Repbase and Dfam: a new open foundation for transposable element research*. Mobile DNA, 17, 17. <https://doi.org/10.1186/s13100-026-00409-9>（PMID 42343465）。GIRI 声明：<http://www.girinst.org/repbase/repbase_and_dfam.html>

## 项目状态

转换流程已通过 RepBase 31.09 的格式验证与 RepeatMasker 功能测试，作为独立开发、可复现的转换工具分发。下游基因组注释表现需按具体应用单独评估。

**Release v1.0.0 —— 已完成 RepBase 31.09 与 RepeatMasker 4.2.4 的转换验证与功能兼容性验证。** 此处 “Validated” 指通过本项目自身的测试（转换完整性、RepeatMasker 兼容性、真实基因组行为比较），**不**代表 GIRI 官方认证，也不代表全基因组注释准确率保证。
