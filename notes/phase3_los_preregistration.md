# Phase 3 Preregistration: Locking the Line-of-Sight Pair Cut

Status: DRAFT — awaiting advisor approval. Per the preregistration
discipline, the scientific sweep will not be produced or examined until the
definitions and targets below are approved. (Tooling has been built and
functionally tested on throwaway subsamples; a timing pilot measures compute
cost only.)

## Purpose

Choose a single fixed `rpar` cut (maximum line-of-sight separation for pair
selection, in redshift space) to be used identically in the simulation grid
and the BOSS reruns. The choice is made from predefined completeness and
contamination criteria measured in the simulation, where both redshift-space
and true LOS coordinates are known per pair — NOT from the strength or
stability of the bridge signal, which would risk tuning the selection to the
result. Bridge stacks are produced afterward only to document the
consequence of the choice.

## Proposed definitions (for approval)

- **Primary association criterion**: a pair is physically associated if its
  true (real-space) LOS separation satisfies `|r_par,real| <= 20 h^-1 Mpc`.
  Note this is a real-space *cylinder* membership definition, not a
  boundness claim; it still admits real-space chance alignments within the
  cylinder.
- **Robustness criterion**: true 3D separation `<= 2.0 x r_perp`.
  Classifications are repeated under this definition to check that the
  chosen cut is not an artifact of the primary criterion.
- **Candidate cuts**: `rpar,rsd <= 5, 10, 20, 30 h^-1 Mpc`.
- **Completeness target**: `>= 90%` of truth pairs retained.
- **Interloper ceiling**: `<= 10%` of selected pairs failing the primary
  criterion.
- **Decision rule**: the smallest candidate cut satisfying both. If no
  candidate satisfies both, do not choose; report the table and revisit the
  targets together before any further analysis.

## Measurement design

- **Truth catalog**: pair search in *real space* (`--rpar-space real`,
  cut = the primary criterion itself), same transverse bin
  (`r_perp = 8-12 h^-1 Mpc`, mass13 sample). Completeness at cut C is the
  fraction of truth pairs whose *redshift-space* separation is `<= C` —
  measured before any RSD cut, so there is no completeness-denominator
  truncation.
- **Deep RSD catalog**: production-like pair search in redshift space with
  `rpar <= 50 h^-1 Mpc`; nested subsets at the candidate cuts. Interloper
  fraction at C is the fraction of selected pairs failing the association
  criterion. A boundary check reports truth pairs near the 50 edge.
- **Consequence stacks** (after the choice): fixed-separation pair stacks of
  each nested subset on the unsmoothed and pre-smoothed 8-arcmin
  0.1 h^-1 Mpc maps, each with a single-halo template from the *same* map
  (matched-template rule); residual bridge statistics reported vs cut.

## Tooling (built, functionally tested)

- `analysis/sim/find_pairs_sim.py --rpar-space {redshift,real}`; catalogs
  now record `r_parallel_real` and `r_parallel_rsd` explicitly.
- `analysis/sim/los_contamination_table.py` — the decision table.
- `analysis/sim/run_los_sweep.py` + `summarize_los_sweep.py` — consequence
  stacks with matched templates.

## Execution order after approval

1. Truth + deep pair searches (two runs).
2. `los_contamination_table.py` -> decision table -> apply the rule.
3. Record the locked cut in `publication_forward_plan.md`.
4. Consequence stacks and summary.
5. Only then: BOSS reruns at the locked cut (separate sign-off on scope).
