#!/usr/bin/env bash
# 同步 bytedance/Protenix 的 infer_json_format.md 与 examples/ 到 docs/protenix2/
# 语义: 只覆盖同名文件, 不删除任何本地文件 (严格 cp, 无 rsync --delete)
# 浅克隆 + 稀疏检出, 不拉上游 git 历史
# 用法:
#   ./scripts/update_protenix2_docs.sh            # 默认 main (GitHub 当前内容)
#   ./scripts/update_protenix2_docs.sh v2.0.0     # 锁上游 tag
#   ./scripts/update_protenix2_docs.sh <ref>      # 任意上游 ref (branch/tag/commit)
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

REF="${1:-main}"
UPSTREAM="https://github.com/bytedance/Protenix.git"
SRC_DOCS="docs/infer_json_format.md"
SRC_EXAMPLES="examples"
DST="docs/protenix2"
STATE="scripts/protenix2_upstream.ref"

echo "==> 上游 ref: $REF"

TMP="$(mktemp -d "${TMPDIR:-/tmp}/protenix-docs.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

# 只下载 commits+trees + docs/examples 的 blob (~26 MB), 无历史
git clone --depth 1 --filter=blob:none --sparse --no-checkout --branch "$REF" "$UPSTREAM" "$TMP"
git -C "$TMP" sparse-checkout set docs examples
git -C "$TMP" checkout

mkdir -p "$DST/examples"
cp -f "$TMP/$SRC_DOCS" "$DST/infer_json_format.md"
cp -rf "$TMP/$SRC_EXAMPLES/." "$DST/examples/"

SHA="$(git -C "$TMP" rev-parse HEAD)"
printf '%s %s\n' "$REF" "$SHA" > "$STATE"

echo "==> 完成: docs/protenix2 已对齐 $REF ($SHA)"
echo "    记录: $STATE (未提交, 请自行 review 后 git add/commit)"
