# Retired examples

Nothing in this directory is a result.

These scripts were proof-of-concept chains, built to show that the
independently validated components compose into a working link. They did
that. They were never literature-validated, so their numbers could only
ever be checked against themselves, and the protocol replications
replaced them:

- `analysis/val_gobby/` — Gobby et al. 2004, time-bin, QBER vs distance
- `analysis/validation/validate_duplinskiy_*.py` — Duplinskiy et al. 2017,
  polarisation

They are kept because they run and because the reasoning written into
them is worth reading. Nothing in `src/` or `run_all.py` imports from
here, and **no entry in the harness roster regenerates or checks any of
it**. Treat every number below as an illustration.

## What is committed here, and what generates it

| Artifact | Written by | Checked by |
|---|---|---|
| `val_system/val_system--seed42.csv` | `val_system.py` | nothing |
| `val_system/val_system--seed42.png` | `val_system.py` | nothing |
| `val_system/val_system--seed42_table.csv` | **nothing** | nothing |

That third row is not a typo. `val_system--seed42_table.csv` is an
orphan: it was produced by a table-extraction pass (`49de71a`) whose
roster later dropped `val_system`, and no script in the tree writes it
now. It cannot be regenerated, cannot be checked, and cannot carry a
provenance header from a generator that no longer exists. It survives
only because it was committed before its producer went away.

It is left in place rather than deleted, because deleting committed
output is the repository owner's call and the file is harmless once
labelled. If it is ever removed, this row is the reason.

## The one claim these scripts used to carry that has since moved

`val_system.py` argued that a time-bin chain is blind to birefringence.
That is now measured properly by
`analysis/validation/validate_gobby_impairments.py`, and measured on the
right topology — this chain builds a **balanced** interferometer, where
both time bins ride one polarisation so a rotation cancels. The Gobby
replication runs the **polarisation-multiplexed** topology, where a
rotation does reach the arms. The blindness is a property of a topology,
not of an encoding, and only the newer script makes that distinction.
