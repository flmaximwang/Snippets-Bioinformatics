#!/usr/bin/env bash
# 同步 aurekaresearch/OpenDDE 的 docs/inference_instructions.md 与 examples/ 到 docs/opendde/
# 语义: 只覆盖同名文件, 不删除任何本地文件 (严格 cp, 无 rsync --delete)
# 浅克隆 + 稀疏检出, 不拉上游 git 历史
# 用法:
#   ./scripts/update_opendde_docs.sh            # 默认 main (GitHub 当前内容)
#   ./scripts/update_opendde_docs.sh v0.1.0     # 锁上游 tag
#   ./scripts/update_opendde_docs.sh <ref>      # 任意上游 ref (branch/tag/commit)
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

REF="${1:-main}"
UPSTREAM="https://github.com/aurekaresearch/OpenDDE.git"
SRC_DOCS="docs/inference_instructions.md"
SRC_EXAMPLES="examples"
DST="docs/opendde"
STATE="scripts/opendde_upstream.ref"

echo "==> 上游 ref: $REF"

TMP="$(mktemp -d "${TMPDIR:-/tmp}/opendde-docs.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

# 只下载 commits+trees + docs/examples 的 blob (~26 MB), 无历史
git clone --depth 1 --filter=blob:none --sparse --no-checkout --branch "$REF" "$UPSTREAM" "$TMP"
git -C "$TMP" sparse-checkout set docs examples
git -C "$TMP" checkout

mkdir -p "$DST/examples"
cp -f "$TMP/$SRC_DOCS" "$DST/inference_instructions.md"
cp -rf "$TMP/$SRC_EXAMPLES/." "$DST/examples/"

# 截断所有 *.a3m 到前 100 行, 避免仓库过大
find "$DST/examples" -type f -name '*.a3m' -print0 | while IFS= read -r -d '' f; do
  head -n 100 "$f" > "$f.head" && mv "$f.head" "$f"
done

SHA="$(git -C "$TMP" rev-parse HEAD)"
printf '%s %s\n' "$REF" "$SHA" > "$STATE"

echo "==> 完成: docs/opendde 已对齐 $REF ($SHA)"
echo "    记录: $STATE (未提交, 请自行 review 后 git add/commit)"
