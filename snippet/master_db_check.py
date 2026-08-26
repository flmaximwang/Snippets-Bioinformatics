#!/usr/bin/env python3
"""Detect truncated/corrupt MASTER binary target .pds files.

Run with multiprocessing (default: all cores). Usage:

    python3 master_db_check.py /mnt/data/public/MASTER/pds [--njobs N]

Binary layout (little-endian; ProteinStruct/TargetStruct, MASTER v1.6):
  header ints (bytes 0..27):
    0  BB coords offset, 1 CA coords offset, 2 numres,
    3  pdbinfo offset,   4 seq offset,       5 dihed distr offset,
    6  dist distr offset (_os_distdistr)
  distance distribution section @ _os_distdistr:
    +0 double dcut, +8 double dstep,
    +16 int[numres*numbin] offset table (numbin=ceil(dcut/dstep)),
    then per non-empty bin: int numelem + numelem ints at _os_distdistr+off
Expected total size = _os_distdistr + 16 + numres*numbin*4
                      + sum over non-empty bins of (1+numelem)*4

Root cause this checks for: an interrupted `createPDS` write (kill / OOM /
power loss) truncates the .pds mid-write. The distance-distribution section is
written LAST, so a partial write leaves the offset table pointing past EOF and
`master` search fails with "Error: could not read distance distribution"
(TargetStruct.cpp readDistDistr). The DB build's incremental skip tests only
file EXISTENCE (`[ ! -e $pout ]`), so an existing-but-truncated .pds is never
rebuilt.
"""
import os
import sys
import struct
import math
from multiprocessing import Pool

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


def scan(path):
    """Return (path, reason) if the file is truncated/corrupt, else (path, None)."""
    sz = os.path.getsize(path)
    try:
        with open(path, 'rb') as f:
            if sz < 28:
                return path, "too-small(%d)" % sz
            data = f.read()
    except OSError as e:
        return path, "io-error(%s)" % e

    try:
        hdr = struct.unpack_from('<7i', data, 0)
    except struct.error:
        return path, "short-header"
    _bb, _ca, numres, _pdb, _seq, _dihed, dist = hdr
    if dist <= 0 or dist + 16 > sz:
        return path, "bad-dist-offset(%d)" % dist
    if dist + 16 + 16 > sz:
        return path, "truncated(dist-section-head beyond eof)"
    dcut, dstep = struct.unpack_from('<2d', data, dist)
    if not (dstep > 0):
        return path, "bad-dstep(%r)" % dstep
    numbin = int(math.ceil(dcut / dstep))
    tbl = dist + 16
    ntab = numres * numbin
    if tbl + ntab * 4 > sz:
        return path, "truncated(offset-table beyond eof)"
    offs = struct.unpack_from('<%di' % ntab, data, tbl)
    expected = tbl + ntab * 4
    for o in offs:
        if o == 0:
            continue
        p = dist + o
        if p + 4 > sz:
            return path, "truncated(data@%d beyond eof %d)" % (p, sz)
        numelem = struct.unpack_from('<i', data, p)[0]
        if numelem < 0 or p + 4 + numelem * 4 > sz:
            return path, "truncated(numelem=%d @%d beyond eof %d)" % (numelem, p, sz)
        expected += (1 + numelem) * 4
    if sz != expected:
        return path, "size-mismatch(actual=%d expected=%d diff=%+d)" % (sz, expected, sz - expected)
    return path, None


def iter_pds(root):
    for dirpath, _dirs, files in os.walk(root):
        for fn in files:
            if fn.endswith('.pds'):
                yield os.path.join(dirpath, fn)


def main(argv):
    if len(argv) < 2:
        sys.stderr.write("usage: %s <pds_root> [--njobs N]\n" % argv[0])
        return 1
    root = argv[1]
    njobs = os.cpu_count() or 1
    if '--njobs' in argv:
        njobs = int(argv[argv.index('--njobs') + 1])

    paths = list(iter_pds(root))
    n = len(paths)
    print("scanning %d .pds files with %d workers..." % (n, njobs), file=sys.stderr)

    bad = []
    with Pool(njobs) as pool:
        it = pool.imap_unordered(scan, paths, chunksize=256)
        if tqdm is not None:
            it = tqdm(it, total=n, desc="scanning pds", unit="file", leave=False)
        for path, reason in it:
            if reason:
                bad.append((path, reason))

    bad.sort()
    print("scanned %d .pds files" % n)
    print("truncated/corrupt: %d" % len(bad))
    for p, r in bad:
        print("%s\t%s" % (p, r))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
