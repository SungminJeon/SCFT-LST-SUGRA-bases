round10plus_cap9_excluded — evidence for the 9-round phase2 chain cap
=====================================================================
2026-07-01.

The production phase2 chain stops at 9 rounds (run_pipeline*.sh -> run_chain.sh
MAX=9). These are the gravity blocks the chain would produce PAST that cap,
captured by a full T_H<=10 chain re-run in the clean-final lineage
(--use-lst-T --no-2mix --no-canonical, chain MAX=13, --save-nonsugra; the chain
converges by itself at round 13).

  400 raw / 399 distinct (by IF eigenvalue signature), ALL at T_H<=9
  (T_H=1: 383, T_H 2..5,9: 16).  object counts: 10->365, 11->33, 12->1.
  Round profile (T_H<=10):  r9=2590 -> r10=365 -> r11=33 -> r12=1 -> r13=0.

ALL 399 ARE INVALID, hence excluded from the catalog
----------------------------------------------------
They are deep su2/su3-stacking on the (-1) curves of a minimal base (mostly
T_H=1, which has the most (-1) room), filled to the c.c. budget (8 per -1;
c(su2)=1, c(su3)=2). A genuine SUGRA base attaching this many su2/su3 needs a
matching number of su2n3 / su2n3mix partners, and those attach at most 6:

  - su2 >= 7  ->  needs > 6 su2n3 partners  ->  impossible        (378 blocks)
  - su2 <= 6  ->  but 0 su2n3/su2n3mix partners present           ( 21 blocks)

  -> 0 / 399 survive the su2n3 partner requirement.

So the 9-round cap loses NO valid blocks: it is a justified practical bound, not
a truncation of real data. (This is a practical post-filter the pipeline does not
itself enforce — the chain passes them on gravitational-anomaly + NHC + signature
alone — but they are trivially excluded by inspection, as here.)

Files
-----
  by_combo/<combo>.cat   the 400 raw blocks, grouped by external combo
  summary.csv            per block: T_H, combo, n_objects, n_su2, n_su3,
                         n_su2n3type, verdict (all INVALID)

NOT part of the catalog. Kept only as evidence for the cap-9 rationale.
