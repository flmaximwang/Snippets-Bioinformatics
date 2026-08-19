#!/bin/bash
# Daily rsync of entire PDB archive (PDB format)
# Logs to /mnt/data/public/RCSB/rsync_pdb.log

SRC="rsync://rsync.rcsb.org:33444/ftp_data/structures/divided/pdb/"
DST="/mnt/data/public/RCSB/"
LOG="/mnt/data/public/RCSB/rsync_pdb.log"

# 用 flock 防止重叠运行（万一前一次没跑完）
exec 200>"$DST/rsync_pdb.lock"
flock -n 200 || { echo "[$(date)] Previous rsync still running, skipped" | tee -a "$LOG"; exit 1; }

echo "=== $(date): Starting rsync ===" | tee -a "$LOG"
rsync -rlpt -v -z --delete --port=33444 "$SRC" "$DST" 2>&1 | tee -a "$LOG"
RSYNC_EXIT=${PIPESTATUS[0]}
echo "=== $(date): rsync exit code $RSYNC_EXIT ===" | tee -a "$LOG"
exit $RSYNC_EXIT
