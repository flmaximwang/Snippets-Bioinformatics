#!/usr/bin/env python3
"""Print the current timestamp and copy it to the clipboard automatically.

Usage:
    timestamp.py                # default format, auto-copies to clipboard
    timestamp.py "%Y-%m-%d"     # custom strftime format, auto-copies
    timestamp.py -n             # print only, do not touch the clipboard
    timestamp.py -h             # help
"""

import argparse
import platform
import subprocess
import sys
from datetime import datetime

DEFAULT_FMT = "%Y-%m-%dT%H:%M:%S"


def _write_stdin(cmd, text, shell=False, stderr=None):
    """Pipe text into cmd's stdin; return True if it exits 0."""
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=shell, stderr=stderr)
    proc.communicate(input=text.encode("utf-8"))
    return proc.returncode == 0


def copy_to_clipboard(text):
    """Copy text to the system clipboard. Returns True on success."""
    system = platform.system()
    if system == "Darwin":
        return _write_stdin(["pbcopy"], text)
    if system == "Windows":
        # `clip` is a cmd.exe builtin, so the shell is required.
        return _write_stdin("clip", text, shell=True)
    if system == "Linux":
        # Try X11 tools first, then Wayland; fall through on missing tool.
        candidates = [
            ["xclip", "-selection", "clipboard"],
            ["xsel", "--clipboard", "--input"],
            ["wl-copy"],
        ]
        for cand in candidates:
            try:
                if _write_stdin(cand, text, stderr=subprocess.DEVNULL):
                    return True
            except FileNotFoundError:
                continue
        return False
    return False


def main():
    parser = argparse.ArgumentParser(
        prog="timestamp.py",
        description="Print the current timestamp and copy it to the clipboard.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--clipboard",
        action="store_true",
        default=argparse.SUPPRESS,
        help="copy to clipboard (on by default; kept for compatibility)",
    )
    parser.add_argument(
        "-n",
        "--no-copy",
        action="store_true",
        default=argparse.SUPPRESS,
        help="print only; do not copy to clipboard",
    )
    parser.add_argument(
        "format",
        nargs="?",
        default=DEFAULT_FMT,
        help="strftime format string",
    )
    args = parser.parse_args()

    try:
        timestamp = datetime.now().strftime(args.format)
    except ValueError as e:
        parser.error(f"invalid format string: {e}")

    print(timestamp)

    if not getattr(args, "no_copy", False):
        if copy_to_clipboard(timestamp):
            print("(copied to clipboard)", file=sys.stderr)
        else:
            print(
                f"(failed to copy to clipboard on {platform.system()})",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
