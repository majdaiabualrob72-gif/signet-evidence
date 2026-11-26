# Use of generative AI in this work

This record supports the declaration in the manuscript (Section 6.4 and Declarations) and follows the publisher's request for the tool, version, dates, access route and a full description of use.

## Tools

- Main tool: Claude Opus 5.5 (Anthropic; model identifier claude-opus-5-5), a large language model, accessed without modification through claude.ai in an agentic workspace with a Linux shell, Python, web search and page fetching, a scholarly search connector and a document editor.
- Dates: 17 to 19 November 2026.
- Fourth tool, for the revision of 22 November 2026: Claude Fable 5.2 (Anthropic; model identifier claude-fable-5-2), in the same claude.ai conversation and workspace. At about 16:50 Doha time that day the model of the conversation was switched to Claude Opus 5.5, which completed the final checks of the built files, the records of the revision and the package.
- Second tool: ChatGPT (OpenAI; model version not recorded). It helped the author to run the reproduction commands, executed part of the checks in its own environment and compiled `REPRODUCTION_LOG.md`. The author also used it to produce the audit of version 8, to write the R reconstruction reported in version 8, and to review the code and the mathematics, including the proofs; no record of that review was kept, so it is disclosed but not cited as evidence.
- Third tool, used only to guide the re-run of version 1.1.0 on 19 November 2026: Claude Sonnet 5.5 (Anthropic; model identifier claude-sonnet-5-5), in the same claude.ai conversation. It gave the set-up and run instructions and read the logs the author uploaded; it executed nothing on the author's computer.
- No model is an author. The author is accountable for all content.

## What each party did

| Task | Done by the model | Done by the author |
|---|---|---|
| Research problem and programme | Proposed the methods in response to the prompts | Wrote the prompts (S8a to S8f), set aims and constraints, chose to pursue the work |
| Methodology | Invented DECODE and then SIGNET; chose the rung structure, certificates, test and tolerance construction | None beyond the prompts |
| Mathematics | Stated and proved Theorems 1 to 3 and Propositions S1 to S4; re-derived them in the final audit | Not recorded in this archive |
| Literature | Searched (SciSpace and web), read abstracts and pages, built the novelty boundary and reference list | Not recorded in this archive |
| Real-data extraction | Wrote and froze the protocol, extracted estimates from opened abstracts with locations, wrote amendments | Not recorded in this archive |
| Code | Wrote the package, tests, simulation and analysis scripts, audit scripts and figure code | Re-ran the tests, the real-data analysis, the simulation tables and the audit checks on a second computer (macOS) on 18 November 2026, with assistance from ChatGPT; all checks passed (`REPRODUCTION_LOG.md`). Ran `scripts/reproduce_v110.sh` for version 1.1.0 on 19 November 2026; all checks passed (`reproduction_logs/REPRODUCTION_v1.1.0_20261119_0940.log`) |
| Analyses | Ran every simulation and analysis in its workspace | Re-ran the subset listed above |
| Figures | Generated all figures with Python code | None |
| Writing | Drafted the manuscript, supplement, cover letter, declarations, simulated peer review and this package | Not recorded in this archive |

## Known risks of AI-generated research and how they were addressed

| Risk | Mitigation in this work | Remaining action |
|---|---|---|
| Fabricated or inaccurate references | Only references located during the work were cited; on 22 November 2026 every reference was checked against Crossref, Europe PMC, DataCite or the publisher's page, and 18 entries were corrected or completed (`REFERENCE_VERIFICATION.md`) | Reading of the cited passages against the full texts by the author |
| Errors in proofs | Proofs written in full (Supplementary S1) and re-derived in the final audit; results checked numerically by unit tests | Independent re-derivation |
| Implementation errors | 50 unit checks, independent-solver verification, independent audit checks, fixed seeds, archived outputs, re-runs of versions 1.0.0 and 1.1.0 on a second computer, and a second implementation in R written by the same model | Code review or re-implementation by an independent person |
| Extraction errors | Each value recorded with its access URL and location; protocol frozen and hashed before extraction; hormone-therapy and vitamin E values matched the published abstracts in a spot check | Check of every value against full texts |
| Overstated claims | Executed and proposed work are labelled; negative results reported; post hoc analyses labelled | Reading by independent investigators |
| Self-evaluation bias | Stress tests built to break the method; every failure and redesign reported | Independent re-analysis |

## Files

- `prompts/`: the prompts, verbatim, with hashes and the accompanying chat instructions.
- `experiments/realdata/PROTOCOL_realdata_v1.md` and `.sha256`: the protocol frozen before extraction.
- `REPRODUCTION_LOG.md`: the re-run of version 1.0.0 on a second computer.
- `REPRODUCTION_RECORD_v1.1.0.md` and `reproduction_logs/`: the runs of version 1.1.0 in the development environment and on the author's computer.
- The conversation transcripts are not part of this archive.

## Pre-submission audit (GENESIS V, 17 to 18 November 2026)

The same model audited the manuscript against a structured editorial prompt (S8d): it revised overclaims ("exact inference" became finite-sample inference with a conservative reported set), added the hypotheses missing from Corollary 1 and Theorems 2 and 3, integrated five newly found related works (references 17, 18, 22, 28 and 39), renumbered the references, moved the Highlights before the References as the journal requires, removed a tint from Figure 4, rebuilt the reference verification list, compared three journals and assembled the package.

## Final scientific gate (GENESIS VI, 18 November 2026)

A sixth prompt (S8e) asked the same model to act as a final gate. It re-derived the theorems, re-ran every test, added independent checks (`audit/`), regenerated every simulation table, re-ran the real-data analysis with a finer scan and spot-checked extracted estimates against the published abstracts. It found and corrected three errors: an overstated description of the cross-solver check, a reporting default that the package did not yet implement, and an edge case that the simulation code handled incorrectly. None changed a reported simulation number. Because the auditor is the model that produced the work, this audit is not independent.

## Re-run on a second computer (18 November 2026)

The author re-ran the four test scripts, the real-data analysis, the default sets, the simulation table generators, the Table 3 audit, the independent checks, the default-equivalence check and two short simulation runs on macOS (Python 3.13.7 with the pinned NumPy, SciPy, CVXPY and Clarabel). All checks passed; the qualifications are listed in `REPRODUCTION_LOG.md`. This establishes computational reproducibility, not the correctness of the logic the code implements.

## Final package (18 November 2026)

Figure 1 was redrawn as a step-by-step guide with the sigmoidoscopy question worked through each step, the example outputs were regenerated so that they show the default, and this README gained a guide to running a new review. The analysis code (`signet/`, `tests/`, `audit/`, `experiments/`) is byte-identical to the version that was re-run on the second computer.

## Version 8 (18 November 2026)

At the author's request, the same model shortened the main text (from about 8,400 to about 6,100 words), opened it with the hormone-therapy case, added Box 1 (SIGNET with ROBINS-I direction judgements), wrote the R interface (`r/`), drafted a pre-submission enquiry, a list of suggested reviewers and a prospective-test protocol, and moved the longer version 7 text of condensed sections to Supplementary S10. A separate implementation of the worked example in R, reported in `audit/independent_R_check_report.md`, reproduced the reported limits; who wrote it, and whether an AI tool assisted, is stated in the manuscript by the author. The analysis code (`signet/`, `tests/`, `audit/` scripts, `experiments/`) is unchanged.

## Version 9 and software version 1.1.0 (18 November 2026)

The author supplied a structured audit of version 8, produced with ChatGPT (OpenAI; model version not recorded). In response, the same model (Claude Opus 5.5):

- changed the inference code (`signet/core.py`, `signet/spec.py`; see `CHANGELOG.md`): explicit bracketing of outer limits with warnings, an outer set with exact end points, a dense check of the hull by default, a validated upper bound for directional p-values, unit-length constraint rows, a certificate-uniqueness check, certificates frozen at a chosen covariance, a dual-feasibility check that returns the non-homogeneous offset, correlations between cells and expected standard errors for a frozen certificate choice (`certificate_se`) in the specification;
- wrote `tests/test_v11.py` (14 checks) and re-ran the 23 version 1.0.0 checks;
- wrote and ran `experiments/realdata/v11_checks.py` (real data under version 1.1.0, inversion and certificate checks, correlation sensitivity), `experiments/validation/validation_v11.py` (4,000 base literatures under five settings), `experiments/validation/validation_default_v11.py` (the reporting default under archived and hardened settings), `experiments/validation/paired_reruns.py` (per-literature re-runs of the parts that had kept only summaries) and `experiments/validation/paired_contrasts.py` (paired re-analysis of Table 4);
- added `outer_scan=0` to the archived benchmark scripts, so that they reproduce version 1.0.0 under version 1.1.0;
- found, during the re-runs, that the missing-controls stress scenario was not reproducible across runs (a loop over a Python set), fixed the loop, regenerated the scenario and the tables and figures built from it, and kept the version 1.0.0 outputs;
- verified the updated reference 26 and five new references on opened pages, and revised the manuscript, supplement (new section S11), title page, cover letter, README and R README.

All of this ran in the model's development environment. The author's re-run of part of it is recorded below. The wording about the R reconstruction was corrected: its report is archived, its source code is not, and the `r/` folder is an interface that calls the Python engine.

## Re-run of version 1.1.0 by the author (19 November 2026)

The author ran `scripts/reproduce_v110.sh` on a Mac (macOS, Apple M1; Python 3.13.7, NumPy 2.5.5, SciPy 1.18.1, Clarabel 0.11.1), following set-up instructions from Claude Sonnet 5.5. The 37 unit checks and the 7 R-interface checks passed, the real-data checks were regenerated, and 400 of 400 sampled base-simulation analyses agreed with the archive. The full validation, the paired re-runs, the tables and the figures were not re-run. Details and checksums: `REPRODUCTION_RECORD_v1.1.0.md`; log: `reproduction_logs/REPRODUCTION_v1.1.0_20261119_0940.log`.

## Version 9.1 (19 November 2026)

Claude Opus 5.5 updated the manuscript, supplement (S8, S9, S11), title page, cover letter, response to the audit and this archive's documentation to report the re-run above, qualified two statements of finite-sample inference in the manuscript, and the matching ones in `README.md` and `CITATION.cff`, as exact for known variances and approximate when variances are estimated. The code was not changed.

## Final manuscript version (19 November 2026)

At the author's request, Claude Opus 5.5 shortened the manuscript of version 9.1 from about 7,300 to about 6,400 words by tightening sentences; no section, theorem, table, figure, reference or reported number was removed, and the detail moved out of Section 7.7 remains in Supplementary S11. The code was not changed.

## Forensic audit and software version 1.2.0 (22 November 2026)

The author supplied a forensic audit of the final submission package with 29 findings (three major, thirteen detailed, thirteen minor) and a structured prompt instructing the model to resolve each one (`prompts/S8f_forensic_revision_prompt.txt`). The audit was produced by Claude Fable 5.2 (Anthropic) in a separate conversation, the same model that then revised the work. It is therefore neither independent of the work nor a human review.

In response, Claude Fable 5.2 (Anthropic; model identifier claude-fable-5-2), in the same claude.ai conversation:

- reconstructed from the code what the random-effects comparator had done (modified Hartung-Knapp interval in the simulation; Wald interval for the vitamin E pool; single published estimates elsewhere) and recomputed it with three intervals;
- added BA-S with shared bias terms, an ablation with the conditional test alone at the full level, three readings of an assumption group, the proportion of bounded sets and the coverage among them, the share of empty sets, a held-out check of the default on new seeds and three new stress scenarios, real-data sensitivity analyses (inflated variances; trial validity as a named assumption), the weight of the most precise cell under the certificate rule, and a sensitivity analysis of the absolute thresholds (`experiments/audit_v12/`; Supplementary S12);
- changed the package (`CHANGELOG.md`): no directional p-value for an empty set, `default_capacity`, a report that prints the default with its capacity and the set without tolerance, named thresholds and a scale warning; wrote `tests/test_v12.py` (13 checks) and re-ran the 37 earlier checks;
- wrote a second implementation of the method in R (`r/reconstruction/`), which agreed with the package in 80 of 80 quantities. The same model wrote both, so this is a consistency check and not independent verification. The R reconstruction cited in versions 8 to 10 is no longer cited as evidence, because its source was never archived; the author states that ChatGPT (model version not recorded) wrote it;
- checked every reference against Crossref, Europe PMC, DataCite or the publisher's page, corrected or completed 18 entries (among them the author list of reference 17 and the first author of reference 32) and added four references, and added the DOI of the Crossref record to further entries (`REFERENCE_VERIFICATION.md`);
- revised the manuscript, the supplement and Figures 1 and 5, and prepared a correction register, a verification package for an independent statistician and a submission-readiness report.
- worked in four passes: revision; recomputation of every reported number from the per-literature records by a separate script; an adversarial reading of the revised text; and an inspection of every page of the built Word and PDF files, which led to layout corrections (tables on landscape pages, column widths, code lines) and to the repair of PDF copies that had lost their equations.
- reported the unfavourable result of its own held-out check: two designs biased against their stated directions reduced the coverage of the default to 0.710 (Supplementary S12.4; Abstract; Sections 7.4, 8.1, 9 and 11).

The model of the conversation was switched to Claude Opus 5.5 during the fourth pass; Claude Opus 5.5 completed that pass, updated the records (testing and readiness reports, correction register, this file) and assembled the package.

All of this ran in the models' environment on 22 November 2026.

## Author's re-run of version 1.2.0 and statements of the author (22 November 2026)

The author ran `scripts/reproduce_v120.sh` on a second computer (MacBook Air, Apple M1, macOS; Python 3.13.7, NumPy 2.5.3, SciPy 1.18.1, Clarabel 0.11.1; R 4.5.1 with lpSolve 5.6.23, quadprog 1.5.8 and jsonlite 2.0.0). A first attempt at 18:12 stopped because NumPy was not installed (`reproduction_logs/REPRODUCTION_v1.2.0_20261122_1812.log`). In the second run (`REPRODUCTION_v1.2.0_20261122_1813.log`), the 50 unit checks, the four example specifications (0 differences), the samples of every analysis of the forensic-audit revision, the 400 legacy simulation analyses, the 7 checks of the R interface and the 80 comparisons of the R implementation passed. One comparison reported two differences: the lower limit of the group-tolerant rung-3 set for the control outcome of the hormone-therapy question, at correlations 0 and 0.75, lay within 10⁻⁶ of zero on the log scale and differed from the archived values by 4.6 × 10⁻⁸ and 6.1 × 10⁻⁸. This exceeds the strict tolerance of that comparison (10⁻⁹ plus 10⁻⁶ relative) but lies far below the bisection tolerance of the method (10⁻⁴); no reported number changes, and the tolerance was not changed after the run.

The author states that he checked the references and the extracted values against the full texts, contributed ideas to the design of the method, and confirms the statements added on 22 November 2026. After this archive was prepared, the author obtained feedback on the manuscript of version 11 from two AI review tools, ChatGPT (OpenAI) and paperreview.ai; revisions of the manuscript made in response are recorded in the submission records, not in this archive.

No human statistician has independently verified the mathematics, the code or the analyses. The author has left that verification to peer review; the verification package lists the claims, hand calculations and commands for the reviewers, and nothing in this archive should be read as an independent sign-off.
