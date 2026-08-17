# PyMOL extra_fit: `*, ref and c. A` → invalid selections (root cause, verified)

PyMOL 3.1.8 (/opt/envs/pymol). Tested on RPXDock ASYM v1 poses of 3e2c 1-79
(chain A + chain B, 474 atoms each, single state).

## Symptom

```
PyMOL>extra_fit *, ASYM__body_job0_comp0_3e2c_1-79__body_job0_comp1_3e2c_1-79__top0_0 and c. A, \
PyMOL>    method=align, cycles=5, cutoff=2.0, mobile_state=-1, target_state=-1
ASYM__..._job0_top0_0 RMSD =    0.000 (474 atoms)
 ExecutiveAlign: invalid selections for alignment.
```

Prints ONE RMSD line (first object aligned fine), then dies. The selection syntax
is NOT wrong — the crash happens on the reference object's own iteration.

## Root cause — two layers

1. **extra_fit is a wrapper** (`modules/pymol/fitting.py`): it loops over every
   object in the mobile selection and calls `method(mobile='?sele & ?obj',
   target='?sele & ?reference', ...)`. The `reference` argument is handled as a
   **bare object name**:

   ```python
   elif reference in models: models.remove(reference)   # 'X and c. A' != 'X' → not removed
   else: _self.select(sele_name, reference, merge=1)     # reference object STAYS in the loop
   ```

   Passing `obj and c. A` fails the string equality, so the reference object is
   never removed from the loop. (The target expression `?sele & ?ref and c. A`
   IS parsed correctly — chain A is honored. The `and c. A` is not the problem.)

2. **align cannot align a selection to another selection of the same object**
   (C side, `layer3/Executive.cpp:7836`):

   ```c
   vla2 = SelectorGetResidueVLA(G, sele2, use_structure, mobile_obj);
   ```

   and `SelectorGetResidueVLA` (`layer3/Selector.cpp:3561`):

   ```c
   for(SeleAtomIterator iter(G, sele); iter.next();) {
       if(iter.obj == exclude) continue;
   ```

   → the target residue list can NEVER contain atoms of the mobile object. When
   the loop reaches the reference object itself: `align(ref, ref and c. A)` →
   every target atom belongs to mobile_obj → empty target → `nb == 0` →
   `ExecutiveAlign: invalid selections for alignment` → CmdException aborts
   extra_fit.

## Verified behavior (PyMOL 3.1.8)

- `align(obj, obj)`, `align(obj, obj and c. A)`, `align(obj and c. B, obj and c. A)`
  → all `invalid selections for alignment`. Cross-object only.
- `mobile_state=-1 / target_state=-1` is NOT the cause: identical failure with
  states=0; cross-object align with states=-1 works fine. The failure happens
  while building residue VLAs, before any state-dependent coordinate access.
  (Python state 0 = "all states" → C -1; Python -1 → C -2, undocumented — don't pass.)

## Corrected command

```pymol
extra_fit * and not ASYM__body_job0_comp0_3e2c_1-79__body_job0_comp1_3e2c_1-79__top0_0, \
    ASYM__body_job0_comp0_3e2c_1-79__body_job0_comp1_3e2c_1-79__top0_0 and c. A, \
    method=align, cycles=5, cutoff=2.0
```

Verified: all 10 poses align onto ref chain A, RMSD ≈ 0.001 (474 atoms). The
reference object stays fixed; only the other poses move.

| extra_fit(selection, reference, ...) | Result |
|---|---|
| `* and not <ref>`, `<ref> and c. A` | ✅ each pose onto ref chain A (474 atoms) |
| `*`, `<ref>` (bare name) | ✅ runs, but target = WHOLE ref object (both chains) — fit dragged by the other chain |
| `*`, `<ref> and c. A` | ❌ invalid selections (self-alignment on ref object) |
| `* and c. A`, `<ref> and c. A` | ❌ same |

## Same-object chain superposition

To superimpose chain B onto chain A of ONE object:

- `pair_fit obj and c. B, obj and c. A` — WORKS (matches by residue ID + atom
  name; verified RMSD 0.0006).
- `fit obj and c. B, obj and c. A` — does NOT work cleanly: default matchmaker=0
  pre-intersects `((mobile) in (target))` which is empty for different chains →
  "No atoms selected" + garbage return value (6.1e11 observed); matchmaker=1
  hits an internal "no atoms left after refinement" error.
- `align`/`super`/`extra_fit` — invalid selections (exclude rule above).

## Version / source note

Behavior verified on the installed PyMOL 3.1.8 and against
schrodinger/pymol-open-source master C sources (Executive.cpp, Selector.cpp).
Earlier note (2026) claiming `extra_fit *, ref and c. A` works is WRONG: both it
and `* and c. A` fail at the same self-alignment iteration — the `*` form just
prints RMSD lines for the other objects first, so it looks like it ran.
