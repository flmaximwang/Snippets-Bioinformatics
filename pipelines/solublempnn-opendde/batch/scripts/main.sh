#!/usr/bin/env bash
# 入口: 前台运行主循环, stderr 合并进 stdout 并 tee 到 report/v8-0-0.log
# (之前误写为 "&2>1": & 会把循环放后台导致 Ctrl+C 难终止, 2>1 会把 stderr 写进名为 "1" 的文件)
./scripts/main_v8-0-0.sh 2>&1 | tee report/v8-0-0.log
python ./scripts/summary_to_fasta.py report/summary.csv
MARK_RES_IDS=6,7,9,10,13,17,20,21,31,32,35,43,46,47,54,57,64,68,75,78,80,83,109,113,114,117,120,121,124,131,135,137,140,142,146,147,154,157,164,175,178
biorazer plot-seqlogo \
    -i report/summary.fa \
    -o report/summary_seqlogo_a.png \
    --renumber-res 0_1,85_107 \
    --res-id-range 1_85 \
    --first-tick-id 1 \
    --mark-res-ids $MARK_RES_IDS \
    --entropy-cutoff 0
biorazer plot-seqlogo \
    -i report/summary.fa \
    -o report/summary_seqlogo_b.png \
    --renumber-res 0_1,85_107 \
    --res-id-range 107_179 \
    --first-tick-id 1 \
    --mark-res-ids $MARK_RES_IDS \
    --entropy-cutoff 0