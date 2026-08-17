#!/usr/bin/env bash
# =====================================================================
# main_v8-0-0.sh - SolubleMPNN → colabfold-msa → OpenDDE → 一致性检查
#
# 按 AGENTS.md 的 5 步流程, 迭代生成 candidate-8.n.0:
#   1. SolubleMPNN 生成一条设计序列 (seed = SEED_BASE + n)
#   2. 按模板切分位置切成 3 片段 (1-79 : 80-85 : 86-158, 片段用设计序列自身残基)
#      提交 colabfold-msa (--pair-mode paired, 无 template)
#      paired.a3m 做预测 MSA, 提取 query 为 single.a3m
#   3. OpenDDE 预测 (opendde.json: 序列 + paired/single MSA 路径)
#   4. 一致性: ptm >= PTM_CUTOFF 且 backbone RMSD <= RMSD_CUTOFF (对齐 v8-0-0_mutate_conserved)
#   5. 通过则写入 report/summary.csv (marker, sequence)
#
# 目录约定: colabfold-msa 与 opendde 都是按记录名/作业名建子目录的批量工具,
# 统一输出到 BATCH_DIR (data/v8-0-0), 每个 candidate 的结果在
# CAND_DIR = BATCH_DIR/candidate-8.n.0 下 (paired/, unpaired/, seed_101/, ...)。
#
# 循环直到得到 TARGET_CANDIDATES 个可用 candidate。
# 每个步骤有断点检查, 中断后可重跑 (幂等); 单次迭代失败只跳过并继续。
# 用法: 在 batch/ 目录下运行  bash scripts/main.sh
# =====================================================================
set -euo pipefail

# Ctrl+C (SIGINT) / kill (SIGTERM): 立即退出循环。
# 不加 trap 时, 被中断的子进程只会让 run_iteration 返回 1,
# 主循环会当作"失败迭代"继续跑下一个 n, 表现为 Ctrl+C 无效。
trap 'log "收到 SIGINT, 退出循环"; exit 130' INT
trap 'log "收到 SIGTERM, 退出循环"; exit 143' TERM

# ---------------- 配置 ----------------
FIXED_RESIDUES="A1 A2 A3 A4 A5 A8 A11 A12 A14 A15 A16 A18 A19 A22 A23 A24 A25 A26 A27 A28 A29 A30 A33 A34 A36 A37 A38 A39 A40 A41 A42 A44 A45 A48 A49 A50 A51 A52 A53 A55 A56 A58 A59 A60 A61 A62 A63 A65 A66 A67 A69 A70 A71 A72 A73 A74 A76 A77 A79 A81 A82 A84 A85 A107 A108 A110 A111 A112 A115 A116 A118 A119 A122 A123 A125 A126 A127 A128 A129 A130 A132 A133 A134 A136 A138 A139 A141 A143 A144 A145 A148 A149 A150 A151 A152 A153 A155 A156 A158 A159 A160 A161 A162 A163 A165 A166 A167 A168 A169 A170 A171 A172 A173 A174 A176 A177 A179"
DESIGN_PDB="node/v8-0-0_mutate_conserved.pdb"
SOLUBLEMPNN_CKPT="/opt/envs/ProteinMPNN/models/solublempnn_v_48_020.pt"
SOLUBLEMPNN_PY="/opt/envs/ProteinMPNN/bin/python"
SOLUBLEMPNN_RUN="/opt/envs/ProteinMPNN/repo/LigandMPNN/run.py"
BIORAZER="/opt/envs/BioRazer/bin/biorazer"
BIORAZER_PY="/opt/envs/BioRazer/bin/python"
OPENDDE="/opt/envs/OpenDDE/bin/opendde"
# OpenDDE 环境变量 (与 /opt/envs/OpenDDE/etc/conda/activate.d/env_vars.sh 一致)
export OPENDDE_ROOT_DIR="${OPENDDE_ROOT_DIR:-/opt/envs/OpenDDE}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

SEED_BASE="${SEED_BASE:-101}"
TARGET_CANDIDATES="${TARGET_CANDIDATES:-3}"
MAX_ITERATIONS="${MAX_ITERATIONS:-50}"
PTM_CUTOFF="${PTM_CUTOFF:-0.8}"
RMSD_CUTOFF="${RMSD_CUTOFF:-1.5}"
EXCLUDE_LINKER="${EXCLUDE_LINKER:-80-85}"

FRAG1_LEN=79
MID_LEN=6
RENUMBER_SPEC="1_1,86_107"
SUMMARY_CSV="report/summary.csv"

# ---------------- 工具函数 ----------------
log() { echo "[$(date '+%F %T')] $*"; }

# 单次迭代; 返回 0=成功(可能已记录) 1=失败/未达标 2=设计序列重复
run_iteration() {
    local n=$1
    local CAND="candidate-8.${n}.0"
    local BATCH_DIR="data/v8-0-0"          # colabfold-msa / opendde 的批量输出根目录
    local CAND_DIR="${BATCH_DIR}/${CAND}"  # 该 candidate 的结果目录 (由批量工具按名字自动建)
    local SEED=$((SEED_BASE + n))

    log "==== 迭代 n=${n} candidate=${CAND} design_seed=${SEED} ===="
    mkdir -p "${CAND_DIR}"

    # ---------- Step 1: SolubleMPNN 生成设计序列 ----------
    local FA
    FA=$(ls "${CAND_DIR}"/seqs/*.fa 2>/dev/null | head -1 || true)
    if [ -z "${FA}" ]; then
        log "[${n}] SolubleMPNN 运行中 (seed=${SEED})..."
        ${SOLUBLEMPNN_PY} -u "${SOLUBLEMPNN_RUN}" \
            --pdb_path "${DESIGN_PDB}" \
            --out_folder "${CAND_DIR}" \
            --fixed_residues "${FIXED_RESIDUES}" \
            --omit_AA C \
            --model_type soluble_mpnn \
            --checkpoint_soluble_mpnn "${SOLUBLEMPNN_CKPT}" \
            --seed "${SEED}" \
            --batch 1 \
            --number_of_batches 1 || { log "[${n}] SolubleMPNN 失败"; return 1; }
        FA=$(ls "${CAND_DIR}"/seqs/*.fa | head -1)
    else
        log "[${n}] SolubleMPNN 结果已存在: ${FA}"
    fi

    # ---------- Step 1b: 取第 2 条记录(第一条设计), 查重, 切 3 片段 ----------
    local PREP DESIGNED_SEQ NORM_SEQ
    PREP=$("${BIORAZER_PY}" scripts/prepare_msa_input.py \
        -i "${FA}" -o "${CAND_DIR}/msa_input.fa" \
        --record-name "${CAND}" --frag1-len "${FRAG1_LEN}" --mid-len "${MID_LEN}") \
        || { log "[${n}] prepare_msa_input.py 失败"; return 1; }
    DESIGNED_SEQ=$(printf '%s\n' "${PREP}" | sed -n 's/^DESIGNED //p')
    NORM_SEQ=$(printf '%s\n' "${PREP}" | sed -n 's/^NORMALIZED //p')
    if [ -z "${DESIGNED_SEQ}" ] || [ -z "${NORM_SEQ}" ]; then
        log "[${n}] 错误: prepare_msa_input.py 输出异常, 跳过该迭代"
        return 1
    fi
    if [ -f "${SUMMARY_CSV}" ] && grep -q ",${DESIGNED_SEQ}$" "${SUMMARY_CSV}"; then
        log "[${n}] 设计序列已存在于 summary.csv, 跳过该迭代 (换 seed 继续)"
        return 2
    fi
    log "[${n}] msa_input.fa: $(sed -n '2p' "${CAND_DIR}/msa_input.fa")"

    # ---------- Step 2: colabfold-msa (输出到 BATCH_DIR, 生成 BATCH_DIR/<CAND>/ 即 CAND_DIR) ----------
    if [ ! -s "${CAND_DIR}/paired/paired.a3m" ]; then
        log "[${n}] colabfold-msa 提交中 (记录名 ${CAND})..."
        "${BIORAZER}" colabfold-msa \
            -i "${CAND_DIR}/msa_input.fa" \
            -o "${BATCH_DIR}" \
            --pair-mode paired \
            --no-fetch-templates || { log "[${n}] colabfold-msa 失败"; return 1; }
    else
        log "[${n}] MSA 已存在: ${CAND_DIR}/paired/paired.a3m"
    fi
    # 提取 paired.a3m 第一条(query)为 single.a3m (OpenDDE pairedMsaPath 需要)
    "${BIORAZER_PY}" scripts/extract_single_a3m.py \
        --paired-msa "${CAND_DIR}/paired/paired.a3m" \
        -o "${CAND_DIR}/single.a3m" || { log "[${n}] single.a3m 提取失败"; return 1; }
    log "[${n}] single.a3m 提取完成: ${CAND_DIR}/single.a3m"

    # ---------- Step 3: OpenDDE 预测 (输出到 BATCH_DIR, 生成 BATCH_DIR/<CAND>/seed_101/) ----------
    if [ ! -d "${CAND_DIR}/seed_101/predictions" ]; then
        log "[${n}] 生成 opendde.json 并提交 OpenDDE (seeds 101)..."
        "${BIORAZER_PY}" scripts/make_opendde_json.py \
            --template scripts/opendde_template.json \
            --sequence "${NORM_SEQ}" \
            --paired-msa "$(realpath "${CAND_DIR}/single.a3m")" \
            --unpaired-msa "$(realpath "${CAND_DIR}/paired/paired.a3m")" \
            --name "${CAND}" \
            -o "${CAND_DIR}/opendde.json" || { log "[${n}] make_opendde_json.py 失败"; return 1; }
        "${OPENDDE}" pred \
            -i "${CAND_DIR}/opendde.json" \
            -o "${BATCH_DIR}" \
            --need_atom_confidence True \
            --seeds 101 || { log "[${n}] OpenDDE 预测失败"; return 1; }
    else
        log "[${n}] OpenDDE 输出已存在: ${CAND_DIR}/seed_101/predictions"
    fi

    # ---------- Step 4/5: 一致性检查 + 写入 summary.csv ----------
    if [ -f "${SUMMARY_CSV}" ] && grep -q "^${CAND}," "${SUMMARY_CSV}"; then
        log "[${n}] ${CAND} 已在 summary.csv 中记录"
        return 0
    fi
    log "[${n}] 一致性检查 (ptm>=${PTM_CUTOFF}, rmsd<=${RMSD_CUTOFF}, excl linker ${EXCLUDE_LINKER})..."
    local CHECK
    CHECK=$("${BIORAZER_PY}" scripts/check_candidate.py \
        --cand-dir "${CAND_DIR}" --name "${CAND}" \
        --ref-pdb "${DESIGN_PDB}" \
        --renumber-script scripts/renumber_protein_structure_file.py \
        --renumber-spec "${RENUMBER_SPEC}" \
        --ptm-cutoff "${PTM_CUTOFF}" --rmsd-cutoff "${RMSD_CUTOFF}" \
        --exclude-linker "${EXCLUDE_LINKER}") || { log "[${n}] 一致性检查失败"; return 1; }
    printf '%s\n' "${CHECK}"

    if printf '%s\n' "${CHECK}" | grep -q '^VERDICT ACCEPT '; then
        echo "${CAND},${DESIGNED_SEQ}" >> "${SUMMARY_CSV}"
        log "[${n}] ACCEPTED -> ${SUMMARY_CSV} 追加 ${CAND}"
    else
        log "[${n}] REJECTED (ptm/rmsd 未达标)"
        return 1
    fi
    return 0
}

# ---------------- 主循环 ----------------
main() {
    local accepted=0 n=2 rc # n 从 2 开始，因为 8.0.0 是原始设计, 而 8.1.4 在之前被用来做质粒编号了，避免混淆
    mkdir -p "$(dirname "${SUMMARY_CSV}")"
    [ -f "${SUMMARY_CSV}" ] && accepted=$(grep -cE '^candidate-8\.[0-9]+\.0,' "${SUMMARY_CSV}" || true)
    log "启动: 已记录 ${accepted} 个 candidate, 目标 ${TARGET_CANDIDATES}, seed 基 ${SEED_BASE}, 最大迭代 ${MAX_ITERATIONS}"

    while [ "${accepted}" -lt "${TARGET_CANDIDATES}" ] && [ "${n}" -lt "${MAX_ITERATIONS}" ]; do
        if run_iteration "${n}"; then
            rc=0
        else
            rc=$?
        fi
        if [ "${rc}" -eq 0 ] && [ -f "${SUMMARY_CSV}" ] && grep -q "^candidate-8.${n}.0," "${SUMMARY_CSV}"; then
            accepted=$((accepted + 1))
            log "进度: ${accepted}/${TARGET_CANDIDATES} 个可用 candidate"
        fi
        n=$((n + 1))
    done

    log "== 循环结束: 共 ${accepted} 个可用 candidate (迭代至 n=$((n - 1))) =="
    if [ -f "${SUMMARY_CSV}" ]; then
        log "summary.csv 内容:"
        cat "${SUMMARY_CSV}"
    fi
    if [ "${accepted}" -lt "${TARGET_CANDIDATES}" ]; then
        log "警告: 未达到目标 ${TARGET_CANDIDATES} (达到 ${accepted}), 最后一次失败记录见上"
        return 1
    fi
}

main "$@"