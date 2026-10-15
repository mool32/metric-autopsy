# Red flags

*Loaded on demand by the `metric-autopsy` skill. If you see **any** of these, stop and
investigate before making a claim — each has a real example attached.*

| # | Red flag | Gate | Our example |
|---|----------|------|-------------|
| 1 | QC differs > 1.5× between compared groups | 1 | TMS: old males detect fewer genes than young males ([audit](../validation/flagship_audit/REPORT.md)) |
| 2 | Metric uses zero / non-zero expression as a feature | 0 | MI 3-bin: bin 0 = "not detected" |
| 3 | Effect vanishes on a different platform | 6 | ACP: an entropy anticorrelation on SmartSeq2 vanished on 10x (internal checklist, no script) |
| 4 | Effect reverses sign at a different parameter | 0 | ACP: the anticorrelation reversed at low depth (internal checklist, no script) |
| 5 | The "negative control" also shows the effect | 5 | MI: HK coupling also declines in males |
| 6 | The effect vanishes at equal depth (or after n_genes matching) | 2 | MI human skin: all coupling vanishes |
| 7 | Effect only in one subset but you didn't stratify | 1 | TMS MI: GATE 1 flags the male stratum and passes the female one ([audit](../validation/flagship_audit/REPORT.md)) |
| 8 | Metric correlates with a physical/geometric constant | 4 | β tracks conduction geometry |
| 9 | Two metrics correlate but share input data | 0 | ACP: E_intra and E_inter from one matrix |
| 10 | You can't state the mechanism in one sentence | 4 | "coupling decreases" — but *why*? |
| 11 | QC distributions don't overlap between groups | 1–2 | TMS males: every old mouse detects fewer genes than every young one ([audit](../validation/flagship_audit/REPORT.md)) |
| 12 | QC differs in an interaction, not the main effect | 1 | TMS: the drop is in the old males; whether it is a sex-by-age interaction the data cannot establish ([audit](../validation/flagship_audit/REPORT.md)) |

## How to read this table

- **Rows 1, 7, 11, 12** are all the *factorial-QC* failure and are the cheapest to check and the
  most likely to fire. Run GATE 1 first, every time.
- **Rows 2, 4, 9** are *mathematical dependency* (GATE 0): the metric would move even on data with
  no biology in it.
- **Rows 3, 6** are *platform artifacts* (GATES 2/6): the effect is real in the numbers but not in
  the biology; it does not survive a change of technology or equal depth. (Equalize depth by
  thinning, not by selecting cells on n_genes, which biology also moves — and only for a
  composition estimand; a content estimand needs spike-ins.)
- **Rows 8, 10** are *category errors* (GATE 4): you are measuring a physical constant, or you
  cannot name the process — the metric is not about the biology you claim.

A single red flag is not a death sentence; an *unexamined* red flag is. The job of the autopsy is
to convert each flag into an explicit pass/fail with a number attached.
