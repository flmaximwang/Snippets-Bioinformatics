#!/usr/bin/env python3
"""生成指定数量的 UUID → 写入剪贴板

用法:
  ./gen_uuid.py         # 生成 1 个
  ./gen_uuid.py 5       # 生成 5 个
  ./gen_uuid.py 10      # 生成 10 个
"""

import subprocess
import sys
import uuid


def main() -> None:
    count = int(sys.argv[1]) if len(sys.argv) > 1 else 1

    uuids = [str(uuid.uuid4()) for _ in range(count)]
    output = "\n".join(uuids)

    # macOS 剪贴板
    subprocess.run("pbcopy", text=True, input=output)

    print(f"UUID v4 × {count}")
    print(output)
    print("已复制到剪贴板 ✓")


if __name__ == "__main__":
    main()
