# SolubleMPNN → colabfold-msa → OpenDDE 批量设计流水线

EcBfr v8-0-0 的可溶性设计流水线: 迭代调用 SolubleMPNN 生成一条设计序列,
经 colabfold-msa 获得多序列比对, OpenDDE 预测结构, 一致性检查
(ptm + backbone RMSD) 通过后记入 `report/summary.csv`, 直到凑齐目标数量。

来源: `zsqlab08_EcBfr-Mirror/pipelines/pipeline-4.0.1/2_SolubleMPNN/batch/scripts`
(另有 `manual/` 目录对应人工单条流程, 不在本目录内)。

## 流程总览 (5 步)

1. **SolubleMPNN 生成序列** — 固定 FIXED_RESIDUES, 每个 candidate 一条设计序列
   (seed = SEED_BASE + n)
2. **MSA 准备** — 取设计序列, 按模板切分位置切成 3 片段
   (1-79 : 80-85 : 86-158, linker 片段直接用设计序列自身残基), 提交
   colabfold-msa (`--pair-mode paired`, 不用 template); `paired.a3m` 作为
   预测 MSA, 并提取第一条 (query) 为 `single.a3m`
3. **OpenDDE 预测** — 由 `opendde_template.json` 生成 `opendde.json`
   (sequence + `pairedMsaPath`=single.a3m + `unpairedMsaPath`=paired.a3m),
   运行 `opendde pred`
4. **一致性检查** — `ptm >= 0.8` 且对齐 `node/v8-0-0_mutate_conserved.pdb`
   后 backbone RMSD `<= 1.5` (排除 linker 80-85)
5. **记录** — 通过则追加 `report/summary.csv` (`candidate-8.n.0,<序列>`),
   随后 `summary_to_fasta.py` 生成 fasta 并出 seqlogo 图

## 目录约定

```
batch/
├── scripts/              本目录
│   ├── main.sh           入口: 跑主循环 + 生成 seqlogo
│   ├── main_v8-0-0.sh    主驱动 (5 步循环, 幂等可断点续跑)
│   ├── prepare_msa_input.py    取设计序列并切 3 片段
│   ├── extract_single_a3m.py   从 paired.a3m 提取 query
│   ├── make_opendde_json.py    生成 opendde.json
│   ├── check_candidate.py      ptm / RMSD 一致性检查
│   ├── summary_to_fasta.py     summary.csv → fasta
│   ├── opendde_template.json   OpenDDE 输入模板
│   └── renumber_protein_structure_file.py
│                             (symlink → snippet/renumber_protein_structure_files.py)
├── node/                 设计用参考结构 v8-0-0_mutate_conserved.pdb
├── data/v8-0-0/          批量输出根 (colabfold-msa 与 opendde 都按记录名建子目录)
│   └── candidate-8.n.0/  单个 candidate: seqs/, msa_input.fa, paired/, single.a3m,
│                         opendde.json, seed_101/predictions/
└── report/               summary.csv, v8-0-0.log, summary.fa, summary_seqlogo_{a,b}.png
```

colabfold-msa 与 opendde 都是按记录名/作业名自动建子目录的批量工具, 统一以
`BATCH_DIR=data/v8-0-0` 为输出根, 每个 candidate 的结果落在
`CAND_DIR=${BATCH_DIR}/candidate-8.n.0` 下。

## 依赖

三个 micromamba 环境 (驱动内硬编码绝对路径):

- **ProteinMPNN** `/opt/envs/ProteinMPNN` — `repo/LigandMPNN/run.py` +
  权重 `models/solublempnn_v_48_020.pt` (soluble_mpnn 模型)
- **BioRazer** `/opt/envs/BioRazer` — `biorazer colabfold-msa` / `biorazer plot-seqlogo`
- **OpenDDE** `/opt/envs/OpenDDE` — `opendde pred`; 需要
  `OPENDDE_ROOT_DIR` 与 `PYTORCH_CUDA_ALLOC_CONF` (与 env 的
  `activate.d/env_vars.sh` 一致, 驱动内已 export, 默认值可被覆盖)

## 运行

在 `batch/` 目录下执行:

```sh
bash scripts/main.sh
```

前台运行, stderr 合并进 stdout 并 tee 到 `report/v8-0-0.log`; 完成后生成
`summary.fa` 与两张 seqlogo 图。

### 可调参数 (环境变量)

| 变量 | 默认 | 含义 |
|---|---|---|
| `SEED_BASE` | 101 | SolubleMPNN / OpenDDE 的 seed 基数, 实际 seed = SEED_BASE + n |
| `TARGET_CANDIDATES` | 3 | 需要的可用 candidate 数量 |
| `MAX_ITERATIONS` | 50 | 最大迭代次数 (n 上限) |
| `PTM_CUTOFF` | 0.8 | 一致性检查 ptm 阈值 |
| `RMSD_CUTOFF` | 1.5 | 一致性检查 backbone RMSD 阈值 (Å) |
| `EXCLUDE_LINKER` | 80-85 | RMSD 计算排除的 linker 残基 |

示例: `TARGET_CANDIDATES=5 bash scripts/main.sh`

## 各步骤细节

### Step 1 — SolubleMPNN

```sh
FIXED_RESIDUES="A1 A2 A3 ... A179"   # 全长编号 (chain A), 见 main_v8-0-0.sh 顶部
/opt/envs/ProteinMPNN/bin/python -u /opt/envs/ProteinMPNN/repo/LigandMPNN/run.py \
    --pdb_path node/v8-0-0_mutate_conserved.pdb \
    --out_folder data/v8-0-0/candidate-8.n.0 \
    --fixed_residues "$FIXED_RESIDUES" \
    --omit_AA C \
    --model_type soluble_mpnn \
    --checkpoint_soluble_mpnn /opt/envs/ProteinMPNN/models/solublempnn_v_48_020.pt \
    --seed $((SEED_BASE + n)) \
    --batch 1 --number_of_batches 1
```

- 输出在 `CAND_DIR/seqs/*.fa`: 第 1 条记录是参考序列, **第 2 条是第一条设计序列**。
- 设计序列若已存在于 `report/summary.csv`, 跳过该迭代 (换 seed 继续)。

### Step 2 — MSA 准备

```sh
biorazer colabfold-msa \
    -i data/v8-0-0/candidate-8.n.0/msa_input.fa \
    -o data/v8-0-0 \
    --pair-mode paired \
    --no-fetch-templates
```

- `prepare_msa_input.py` 取 seqs/*.fa 第 2 条, 校验长度与参考一致后按
  `frag1_len=79, mid_len=6` 切成 `frag1:mid:frag2` 单行 fasta
  (片段间用 `:` 分隔), 并在 stdout 输出 `DESIGNED` / `NORMALIZED` (供上层解析)。
- **不用 template**: 此前对比过是否使用 template 结果相近, 为简化流程不加。
- MSA 完成后, `extract_single_a3m.py` 把 `paired/paired.a3m` 的第一条
  (query) 提取为 `single.a3m` —— OpenDDE 的 `pairedMsaPath` 需要它。
- `make_opendde_json.py` 生成 `opendde.json`:
  - `sequence` = NORMALIZED (与 MSA query 一致)
  - `pairedMsaPath` = **single.a3m** (绝对路径)
  - `unpairedMsaPath` = **paired/paired.a3m** (绝对路径)
  - ⚠️ 不要把 `paired.a3m` 放进 `pairedMsaPath`

### Step 3 — OpenDDE 预测

```sh
opendde pred \
    -i data/v8-0-0/candidate-8.n.0/opendde.json \
    -o data/v8-0-0 \
    --need_atom_confidence True \
    --seeds 101
```

输出落在 `CAND_DIR/seed_101/predictions/`。

### Step 4 — 一致性检查 (`check_candidate.py`)

1. 读 `seed_*/predictions` 下 `*_summary_confidence_sample_*.json` 的 ptm;
2. 用 `renumber_protein_structure_file.py` 把 sample 的 `.cif` 重编号为 `.pdb`
   (`--renumber-res 1_1,86_107`, 与参考结构编号一致);
3. biotite 叠合对齐 `node/v8-0-0_mutate_conserved.pdb`, 计算 backbone
   (N/CA/C/O) RMSD, 默认排除 linker 80-85;
4. 按 ptm 降序检查, **第一个同时满足** `ptm >= PTM_CUTOFF` 且
   `rmsd <= RMSD_CUTOFF` 的 sample → `VERDICT ACCEPT`。

### Step 5 — 记录与 seqlogo

```sh
python scripts/summary_to_fasta.py report/summary.csv   # → report/summary.fa
biorazer plot-seqlogo \
    -i report/summary.fa \
    -o report/summary_seqlogo_a.png \
    --renumber-res 0_1,85_107 \
    --res-id-range 1_85 \
    --first-tick-id 1 \
    --mark-res-ids 6,7,9,10,...,178 \
    --entropy-cutoff 0
```

- 图 a 覆盖全长编号 1-85 (frag1 + linker), 图 b 覆盖 107-179 (frag2)。
- `MARK_RES_IDS` = FIXED_RESIDUES 在 1-179 中的补集 (再排除不存在的 86-106),
  即**全部可突变(设计)位点**。

## 编号说明

- 设计序列共 158 aa: frag1(1-79) : linker(80-85) : frag2(86-158)。
- 全长 v8-0-0 编号为 1-179, 其中 **86-106 区段在设计中被删去**, 不存在。
- 重编号映射: `1_1,86_107` (设计 1 → 全长 1, 设计 86 → 全长 107)。

## 批处理语义 (main_v8-0-0.sh)

- 迭代 n 从 **2** 开始: 8.0.0 是原始设计, 8.1.4 之前已被用于质粒编号, 避免混淆。
- 每个 candidate 各步骤**幂等** (产物已存在则跳过):
  `seqs/*.fa` → 跳过 SolubleMPNN; `paired/paired.a3m` → 跳过 colabfold-msa;
  `seed_101/predictions` → 跳过 OpenDDE; summary.csv 已有该 candidate → 跳过检查。
- `trap INT/TERM`: 收到 Ctrl+C / kill 立即退出循环——否则被中断的子进程只会让
  该次迭代"失败", 主循环会当成失败迭代继续跑下一个 n, Ctrl+C 表现为无效。
- 单次迭代失败只跳过并继续; 循环直到凑齐 `TARGET_CANDIDATES` 个或达到
  `MAX_ITERATIONS`。

## 已知问题 / 经验

- `main.sh` 曾误写 `2>&1` 为 `"&2>1"`: `&` 会把循环放到后台导致 Ctrl+C 难终止,
  `2>1` 会把 stderr 写进名为 `1` 的文件。
- 依赖三个环境的绝对路径 (`/opt/envs/...`), 换机器需修改 `main_v8-0-0.sh` 顶部配置。
- 本流程的 candidate 命名 (candidate-8.n.0) 与 FIXED_RESIDUES 均为 v8-0-0
  项目专属, 套用到新项目时需整体替换。
