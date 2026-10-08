# Draft: Limitations — what "explained by depth" can and cannot establish

*Draft for the manuscript's limitations section, decided by the project owner on 2026-10-08. Per
the 2026-10-07 plan, the manuscript itself stays untouched until the confirmatory validation
(plan step 3) is done. This text waits here until then.*

> **"Explained by depth" is established only where the design can establish it.** The
> diagnosis needs two things at the level of biological replicates: a raw difference that is
> detected, and a difference at equal depth that is not. With few replicates the first is out of
> reach. In the conditions of probe p11b — a pure depth artifact (identical biology, one group at
> 30% capture), 100 cells per mouse and the reference `mi_3bin` metric — the raw difference was
> detected and explained by depth in 5 of 40 datasets with 4 mice per group, in 25 of 40 with 6 or
> 8, 26 of 40 with 10, 30 of 40 with 12, 34 of 40 with 16, and 40 of 40 (95% Clopper–Pearson interval 91–100%) with
> 20 (`validation/probes/p11b_design.log`, engine `d96d747`). With 4 vs 4 mice the exact
> permutation over replicates reaches p < 0.05 only at complete separation of the mice
> (p = 0.029; the next attainable values are 0.057 and 0.086). On typical designs of three to eight
> replicates per group the tool therefore protects against a false SUPPORTED — no number of mice
> gave one in these 360 datasets — but it returns INCONCLUSIVE instead of proving the artifact. An
> INCONCLUSIVE verdict at such a design means "the claim is not established", not "the
> difference is technical".

## Notes for the methods section

- The number of mice for p11b was chosen by a rule committed before it was run (smallest N whose
  lower 95% bound of P(explained by depth) over 40 datasets is at least 0.80; commits `cda7dd9`,
  `0112bb1`), and the design run itself exposed a bias of the earlier per-cell quantile thinning
  that was fixed before the result above (dev-set journal, D1–D1b).
- The required number of replicates depends on the between-replicate spread of the metric
  relative to the artifact: per-mouse MI from 100 cells is noisy (between-mouse SD 0.04–0.08
  against a raw difference of 0.035–0.075). A metric that is less noisy per replicate, or more
  cells per replicate, needs fewer.
