# Deviations from pre-registration v1

After the tag `v0.3.0-prereg` the engine and the protocol change only through an entry here
(D-series, one entry per change: what, why, when, and which results it touches).

## D1 — the loading of anchor R1's background (2026-10-10; after the results, outside v1)

Written before R1's run, as the owner decided on 2026-10-10.

**What broke.** In v1's run (workflow run 37990113304, job `anchors`) R1 raised
`ValueError('setting an array element with a sequence.')`, and the engine never ran.
- B2.npz stores its counts sparse, and `anchors._load` returned them as a scipy CSR matrix.
- `anchors.r1` passed that matrix to `SimpleData`, whose constructor converts X with
  `np.asarray(X, dtype=float)` (`src/metric_autopsy/core.py:122`). That accepts only a dense array.
- None of R1's three claims has a verdict: Xist, the Y-gene score, and the female-vs-female sham.

**Why it was not caught.** The dry runs' backgrounds are simulated files stored dense, and R2's
B3 is dense too. No test loaded a sparse background through `anchors._load`.
`validation/prereg/test_anchors.py` now does, on synthetic sparse files only. On the tag's
`anchors.py` it fails with this ValueError at `core.py:122`.

**What changes: only the loading in `validation/prereg/anchors.py`.**
- `_load` takes the genes an anchor names, and `r1` names Xist and the four Y genes.
- A dense file (B3) is returned as before.
- A sparse file is cut and made dense by the rule that cut B1 and B2 for the panel (`panel.py`,
  `N_GENES = 2000`): the named genes, then the most expressed by mean count over the file's
  cells (ties by column order), 2,000 genes in all, in column order.
  - Every B2 donor has at least 231 cells, above the panel's 200, so the panel's rule counts
    all of B2's cells.
  - The reason for the cut is memory. B2 dense in float64 is 19,475 × 53,384 × 8 B, about
    8.3 GB, and the engine copies its data for its perturbations, on a 16 GB runner. Cut, it is
    about 0.3 GB. No new rule is introduced.
  - A consequence of the cut, recorded before the run: R1's metric, `log_cp10k_mean`
    (unchanged), divides by the sum of a cell's counts over the genes it is given. R1's library
    size is therefore the sum over the 2,000 genes kept, as for the panel's backgrounds.
- Nothing else in `anchors.py` changes: not the metric, not the claims, not the SESOI, not the
  injection. The allowed outcomes stay exactly R1's row of v1.md section 3.2:
  - SUPPORTED where the engine's replicate rule gives an effect verdict (at least 3 mice of
    each sex), else INCONCLUSIVE;
  - the sham: NO DETECTABLE EFFECT or INCONCLUSIVE;
  - GATE 0's refusal: allowed, never definite.

**What it touches: only R1.**
- v1's results stay as published on `results/panel-v1`: the scores, the criteria and the
  anchors, with R1 "did not run".
- R1 runs once, after the results and outside v1, in `.github/workflows/anchors-r1.yml`. The run
  commit is the tag plus three files of this change: this file, `anchors.py` and
  `test_anchors.py`. The workflow itself cannot be in that commit, because GITHUB_TOKEN cannot push
  a commit that adds a workflow; its blob is recorded with the result.
- The result goes, whatever it is, to the new branch `results/panel-v1-r1`. A second start finds
  the branch and stops. If R1 fails again, it is not fixed without the owner's decision.

**The engine is the tag's, byte for byte.** The workflow checks this before it installs anything,
as `frozen.py` checks the frozen paths: `git diff` against the tag is empty for `src/`.
- `src/` git tree: `97cbab3c437c9ddb118a0690e6e0d24f8a16069f`.
- sha256 of the list of its files' sha256 (`find src -type f | LC_ALL=C sort | xargs sha256sum`):
  `17e9a04f6fdd82e98606da89ad98cdce9a7f22fc35c1b31235efcb0637b132e8`.
- The files:

```
05246af6fa1e8b2213a738ab33bbfa07aa96a76274fed406ae23df2d8ad04221  src/metric_autopsy/__init__.py
50b7734f3e4322a7540397836c74021975a6d295c8891043e7cae210fe4193d9  src/metric_autopsy/cli.py
0ce473e79f1a17107cadc5f85bdd8c39ac8739d5144fd86147d6862898c4efd2  src/metric_autopsy/core.py
0633adf21c93072a9c3d2bce76d5f1147b452fad242b19faeee103731ae2366a  src/metric_autopsy/effect.py
f990991966f23571079c5f155591808173418d98e523c71305898b7a0270253b  src/metric_autopsy/equalize.py
a20a0e611ea5faddf6f9c09db56bbb3e4b259d3688b3510d7b58eefa3ffe05e1  src/metric_autopsy/gates.py
b80a0cba8fa3418fa5f97e1d8e32eb7e4bd9d38eabaccd5631cb8b897f17251a  src/metric_autopsy/injected_signal.py
db0078c17fdeef04efea74357e4f39568371a1d75361c987de1e05cfe239df52  src/metric_autopsy/mcp_server.py
e30aef6a3c7f3f1feeeee6b44bf6f20af44000d5003bb79d915112b47a8e573a  src/metric_autopsy/metrics.py
5569fb1a200e26b384590b66cf47f33c4f4ff325a4fe97a311f82e77de794fe8  src/metric_autopsy/provenance.py
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  src/metric_autopsy/py.typed
383d2b33a419ef063c6d4eed2cf4da8997ddf31ca80d85d23b72ed59d56b7c1e  src/metric_autopsy/qc.py
c4ab62e88ffadc3eb15fdc32d2f1ebc6a87d88e3de1f84b962e28e7ee346ecff  src/metric_autopsy/report.py
8e7b086124ecf435e223c4e65593c7ecfb4c4685449f7ef835ec88d268773a4b  src/metric_autopsy/stats.py
```

**The outcome (written after the run, 2026-10-10).** R1 ran once: workflow run 38033821978,
run commit b0ce92d, result on `results/panel-v1-r1` at 36a04bf (`r1/anchors_r1.json`). It did not
crash. Two of its three claims are outside their allowed sets:

| Claim | Verdict (cause) | Allowed | In the allowed set |
|---|---|---|---|
| R1 Xist | NOT SUPPORTED — metric invalid (GATE 4): the injected 2-fold Xist signal moves the metric by +0.1627 (95% CI +0.1626 to +0.1628), below delta_min 0.25 | SUPPORTED | no |
| R1 Y genes | UNIDENTIFIABLE: B2's donor ids are partially crossed with sex (some hold cells of both sexes), so neither a nested nor a paired analysis is valid | SUPPORTED | no |
| R1 sham (female vs female) | NO DETECTABLE EFFECT (TOST within ±0.5) | NO DETECTABLE EFFECT, INCONCLUSIVE | yes |

- The design is UNIDENTIFIABLE for Xist too (its fields), but metric validity decides first.
- The mice per age and sex that R1 counted add up to 56 over B2's 50 donor ids.
- Nothing was changed after the run, and nothing will be without the owner's decision.
- `test_anchors.py` gained a skip without scipy after the run, for CI's core-only jobs (a6d3e15).
  The run used the version at 965cfc4.

**The causes (an analysis after the outcome, 2026-10-10).** It re-runs no claim and changes no
verdict (`validation/exploratory/b2_sex_structure.py` and `.log`, workflow run 38045140802). The
headline stands: two of R1's three claims are outside their allowed sets. Their causes differ.
- **The Y genes, and the design of the Xist claim: UNIDENTIFIABLE, which reads B2 correctly.**
  B2's 50 donor ids are samples. By the share of their cells with Xist, 34 are single males
  (Fltp_adult, STZ, VSG, spikein_drug), 10 single females (NOD, NOD_elimination) and 6 pools of
  both sexes (Fltp_P16 at 2 weeks, Fltp_2y at 20 months and over). No development stage holds a
  single male and a single female, so within a stage the sexes can be compared only inside the
  pools. The protocol erred: `r1_allowed` counted the ids of each sex over all stages (16 female,
  40 male: a pool counts for both), and checked neither that an id is one animal of one sex nor
  that both sexes share a stratum. By the engine's own design rules R1's allowed set was
  UNIDENTIFIABLE: with the pools the ids are partially crossed with sex
  (`effect._replicate_design`), and without them no stage holds both sexes (GATE 1 stops). The
  source atlas lists the six pools as "mixed", and in them each cell's sex was assigned from a
  score of the Y-chromosome genes (theislab/mouse_cross-condition_pancreatic_islet_atlas at
  3af65e4: `2_annotate_Fltp_P16.py`, `2_annotate_Fltp_2y.py`), so a Y-gene comparison inside the
  pools would be circular.
- **Xist: "metric invalid" by GATE 4, which is false: an engine defect.** GATE 4's module probe
  raises Xist 2-fold in every cell (the other genes thinned to half, against a sham that thins
  every gene) and measures the metric over all cells. 74.2% of B2's cells have no Xist, where
  the probe changes nothing, so the response is +0.1627 over all cells, below delta_min 0.25, and
  +0.6310 over the 5,027 cells with Xist. v1 never tested this path: its panel used the module probe
  only on the useless metric N6c. In the anchors the probe ran on four valid metrics and passed on
  three: R2a's G2M score, R1's Y-gene score, and R1's sham, which is the same Xist metric on the
  female cells only (+0.6149 in the log).
- **The sham, female against female: NO DETECTABLE EFFECT, which is correct.**

v0.1.1 did not confirm Xist because of GATE 1, on TMS (`validation/flagship_audit/REPORT.md`);
v0.3 did not either, because of GATE 4 and the design, on B2. v0.3 did not fix v0.1.1's Xist
error.
