# Draft: "Validation, not benchmarking"

*Draft paragraph for the manuscript's positioning section. Per the 2026-10-07 plan, the
manuscript itself stays untouched until the confirmatory validation (plan step 3) is done.
This text waits here until then.*

> **Validation, not benchmarking.** Benchmarking studies (Weber et al., 2019) and suites
> such as scIB (Luecken et al., 2022) ask a comparative, method-level question: across
> datasets with ground truth, which method performs best? metric-autopsy asks an
> inference-level question: for this number, computed on these data, is the reported
> difference produced by the biological attribute it is claimed to measure, or by the
> measurement process? This is construct validity (Cronbach & Meehl, 1955; Messick, 1995)
> in its causal form: a measure is valid for an attribute if variation in the attribute
> causally produces variation in the measure (Borsboom, Mellenbergh & van Heerden, 2004).
> Because ground truth is rarely available for the claim under test, validity cannot be
> scored; it has to be argued (Kane, 2013). The pre-registration states the interpretive
> argument, and each gate tests one of its inferences against a named rival explanation
> (Platt, 1964; Shadish, Cook & Campbell, 2002). A benchmark can tell you which method is
> best on average; only a validity argument can tell you whether the number in your figure
> means what your sentence says.

## Gate → validity concept (notes for the methods section)

| Gate | Validity concept | Design consequence |
|---|---|---|
| 0 | discriminant validity against method variance (Campbell & Fiske, 1959); robustness (ICH Q2(R2)) | test response to the construct as well as invariance to nuisance; express shifts in affine-invariant units (Stevens, 1946); separate attenuation (reliability; Spearman, 1904) from bias |
| 1–2 | internal validity; balance and matching (Shadish, Cook & Campbell, 2002; Stuart, 2010) | n_genes is downstream of the biology, so matching on it is a bad control (Cinelli, Forney & Pearl, 2022); the correction must follow the pre-registered estimand |
| 4 | construct validity proper; causal theory of validity (Borsboom et al., 2004; Platt, 1964) | operationalize through planted signal (binomial thinning; Gerard, 2020) and a dose–response curve |
| 5 | convergent/discriminant validity; negative controls (Lipsitch, Tchetgen Tchetgen & Cohen, 2010) | many expression-matched negative pairs give an empirical null (Schuemie et al., 2014; Gagnon-Bartsch & Speed, 2012) |
| 6 | external validity / generalizability (Messick, 1995) | require a different platform; the biological replicate is the unit (Squair et al., 2021; Zimmerman et al., 2021) |
| 7 | statistical-conclusion validity (Shadish, Cook & Campbell, 2002) | "no effect" requires an equivalence test against a pre-registered SESOI (Lakens, 2017) |
| pre-reg | argument-based validation (Kane, 2013) | the pre-registration is the interpretation/use argument; each gate tests one inference |

## References to add

- Borsboom D, Mellenbergh GJ, van Heerden J (2004). The concept of validity. *Psychological Review* 111:1061–1071.
- Campbell DT, Fiske DW (1959). Convergent and discriminant validation by the multitrait-multimethod matrix. *Psychological Bulletin* 56:81–105.
- Cinelli C, Forney A, Pearl J (2022). A crash course in good and bad controls. *Sociological Methods & Research*.
- Cronbach LJ, Meehl PE (1955). Construct validity in psychological tests. *Psychological Bulletin* 52:281–302.
- Gagnon-Bartsch JA, Speed TP (2012). Using control genes to correct for unwanted variation in microarray data. *Biostatistics* 13:539–552.
- Gerard D (2020). Data-based RNA-seq simulations by binomial thinning. *BMC Bioinformatics* 21:206.
- ICH (2023). Q2(R2) Validation of analytical procedures.
- Kane MT (2013). Validating the interpretations and uses of test scores. *Journal of Educational Measurement* 50:1–73.
- Lakens D (2017). Equivalence tests: a practical primer for t tests, correlations, and meta-analyses. *Social Psychological and Personality Science* 8:355–362.
- Lipsitch M, Tchetgen Tchetgen E, Cohen T (2010). Negative controls: a tool for detecting confounding and bias in observational studies. *Epidemiology* 21:383–388.
- Luecken MD, et al. (2022). Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods* 19:41–50. *(already cited)*
- Messick S (1995). Validity of psychological assessment. *American Psychologist* 50:741–749.
- Platt JR (1964). Strong inference. *Science* 146:347–353.
- Schuemie MJ, et al. (2014). Interpreting observational studies: why empirical calibration is needed to correct p-values. *Statistics in Medicine* 33:209–218.
- Shadish WR, Cook TD, Campbell DT (2002). *Experimental and Quasi-Experimental Designs for Generalized Causal Inference.* Houghton Mifflin.
- Spearman C (1904). The proof and measurement of association between two things. *American Journal of Psychology* 15:72–101.
- Squair JW, et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications* 12:5692. *(already cited)*
- Stevens SS (1946). On the theory of scales of measurement. *Science* 103:677–680.
- Stuart EA (2010). Matching methods for causal inference: a review and a look forward. *Statistical Science* 25:1–21.
- Weber LM, et al. (2019). Essential guidelines for computational method benchmarking. *Genome Biology* 20:125.
- Zimmerman KD, Espeland MA, Langefeld CD (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications* 12:738.

*Check volume and page numbers against the publisher records before the text goes into the
manuscript.*
