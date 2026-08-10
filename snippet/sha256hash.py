#!/usr/bin/env python3
"""Usage: sha256hash <string>

Outputs the SHA256 hex digest of the input string (converted to UPPERCASE
first) and copies it to clipboard. Supports macOS (pbcopy) and Windows (clip)."""
import hashlib, platform, subprocess, sys

if len(sys.argv) < 2:
    print("Usage: sha256hash <string>", file=sys.stderr)
    sys.exit(1)

input_str = ' '.join(sys.argv[1:]).upper()
digest = hashlib.sha256(input_str.encode('utf-8')).hexdigest()

print(digest)

# Copy to clipboard — cross-platform
system = platform.system()
if system == 'Darwin':
    subprocess.run(['pbcopy'], input=digest.encode('utf-8'), check=True)
elif system == 'Windows':
    subprocess.run(['clip'], input=digest.encode('utf-8'), check=True)
else:
    print("⚠ Clipboard copy not supported on this OS", file=sys.stderr)
    sys.exit(0)

print("✓ Copied to clipboard", file=sys.stderr)
