# Retired examples

Nothing in this directory is a result.

These scripts were proof-of-concept chains, written to show that the
independently validated components compose into a working link. They did
that. They were never checked against a published experiment, so their
numbers could only be compared against themselves, and the protocol
replications replaced them:

- `analysis/val_gobby/` — Gobby et al. 2004, time-bin, QBER vs distance
- `analysis/validation/validate_duplinskiy_*.py` — Duplinskiy et al. 2017,
  polarisation

They are kept because they run and because the reasoning written into
them is worth reading. Nothing in `src/` or `run_all.py` imports from
here, and **no entry in the harness roster runs or checks any of it**.
Treat every number here as an illustration.

## What is committed here

| Artifact | Written by | Checked by |
|---|---|---|
| `val_system/val_system--seed42.csv` | `val_system.py` | nothing |
| `val_system/val_system--seed42.png` | `val_system.py` | nothing |

Both carry a RETIRED banner in their header. They were regenerated on
2026-08-20, and that regeneration moved 14 of 46 rows — the largest
change was 0.27 pp in QBER and 2.5 % in sifted count. The previous
version had been written before the memory fixes in `8220c5b`, so it had
been out of date for two commits. Nothing noticed, because nothing checks
these files.

A third file, `val_system--seed42_table.csv`, used to sit here. No script
in the tree wrote it: it came from a table-extraction pass (`49de71a`)
whose list later dropped `val_system`. It could not be regenerated or
checked, so it was deleted on 2026-08-20. It is in the git history if it
is ever needed.

## One claim these scripts used to make that has moved

`val_system.py` argued that a time-bin chain is blind to birefringence.
That is now measured by
`analysis/validation/validate_gobby_impairments.py`, on the right
topology. This chain builds a **balanced** interferometer, where both
time bins travel on one polarisation, so a rotation affects both equally
and cancels. The Gobby replication runs the **polarisation-multiplexed**
topology, where the two arms leave on orthogonal polarisations and a
rotation does reach them. The blindness is a property of the topology,
not of the encoding, and only the newer script makes that distinction.
