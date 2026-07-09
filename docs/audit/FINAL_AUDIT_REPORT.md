# Final Audit Report

Date: 2026-07-08
Auditor: independent technical review
Project: `tumoronly`
Repository root: `/Users/marianaboroni/Documents/Variants_analysis`

## 1. Executive Summary

This audit round focused on converting the previous partial package-readiness
assessment into a reproducible internal release-readiness artifact.

The main blocking operational issue from the previous round was resolved:

- missing `Suggests` (`httr2`, `Biostrings`, `rtracklayer`) were installed into
  the project-local library `r_libs`;
- `R CMD build .` completed successfully;
- `R CMD check --no-manual tumoronly_0.1.0.tar.gz` completed with `Status: OK`;
- the package was installed locally with `--install-tests`;
- installed-package tests were confirmed to be present under the installed
  package tree.

This materially improves package reproducibility and release hygiene.

It does **not** by itself upgrade the scientific maturity of the tool. The
pipeline is substantially more coherent than the earlier script-based version,
but several capabilities remain experimental or only partially validated for
Brazilian tumor-only use.

## 2. Verdict

**B. pronta para commit em branch de desenvolvimento**

### Justification

Approved for development-branch consolidation because:

- package structure is real and loadable;
- CLI is functional and materially more coherent;
- maftools-based plotting is implemented natively;
- full package build/check now passes in the audited environment;
- installed-package testing is reproducible.

Not upgraded to cohot validation / routine research use because:

- several scientific layers remain explicitly experimental;
- no matched-normal or orthogonal-validation evidence was introduced in this
  round;
- large-scale performance, chunk equivalence, and external/live integrations are
  still not fully audited;
- Brazilian population logic and ancestry remain interpretation aids, not fully
  validated production classifiers.

## 3. Critical Findings

None confirmed in this round.

## 4. High Findings

1. Scientific readiness still lags behind packaging readiness.
   - Evidence: adaptive filtering, ancestry, ML, and Brazilian population logic
     are described and implemented as experimental layers.
   - Severity: HIGH
   - Recommended action: validate against curated Brazilian ground truth with
     matched normal and orthogonal confirmation before any promotion beyond
     development use.

2. Installed-package tests require explicit `--install-tests`.
   - Evidence: `test_check("tumoronly")` on the first local install failed with
     `No test files found`; reinstall with `R CMD INSTALL --install-tests -l r_libs tumoronly_0.1.0.tar.gz`
     materialized the test tree under the installed package.
   - Severity: HIGH
   - Recommended action: document this in contributor/release instructions or
     provide a dedicated audit helper script.

## 5. Medium Findings

1. CLI help text was incomplete relative to implemented commands.
   - Evidence: ancestry-related commands existed in `exec/tumoronly` but were
     not listed in the zero-argument usage text before this round.
   - Severity: MEDIUM
   - Action taken: usage text updated.

2. Test count claims remain imprecise if expressed as a single scalar.
   - Evidence: the repository contains 19 `tests/testthat/test-*.R` files and
     102 `test_that()` blocks by direct count in this environment. Previous
     reports cited larger totals without an auditable counting method.
   - Severity: MEDIUM
   - Recommended action: report test file count, `test_that()` block count, and
     `R CMD check` result explicitly instead of a single ambiguous total.

3. Dual MAF/plotting lineage remains visible in the codebase.
   - Evidence: the newer native maftools flow lives in `R/plots_maftools.R`,
     while legacy/overlapping MAF-oncoplot logic remains in `R/maf_reporting.R`.
   - Severity: MEDIUM
   - Recommended action: finish convergence on one plotting/reporting API to
     reduce maintenance ambiguity.

## 6. Low Findings

1. Local project library now contains additional installed dependencies.
   - Severity: LOW
   - Recommended action: document library bootstrap expectations for fresh clones.

## 7. Architecture

### Status
Partial approval.

### Evidence

- Formal package files exist and load: `DESCRIPTION`, `NAMESPACE`, `exec/tumoronly`.
- Public workflow is orchestrated through `run_tumor_only()`.
- Native maftools plotting is used in `R/plots_maftools.R`.
- Prior architecture observations remain relevant from `docs/REFACTOR_AUDIT.md`.

### Residual concerns

- overlapping reporting/plotting code paths;
- experimental modules integrated into the same package surface;
- some prior audit artifacts still describe older repo states and should be
  versioned or marked historical.

## 8. Inputs

### Status
Partial approval.

### Evidence

- tests cover VCF parsing/build handling and multiple input paths;
- fixture inventory includes VCF and multiple COSMIC-related tables;
- multi-format ingestion is wired through `read_variant_input()`.

### Not fully tested in this round

- large custom TSV/ANNOVAR/Funcotator input matrix;
- malformed/truncated gzip edge cases;
- chunked vs integral equivalence on large files.

## 9. Filters

### Status
Partial approval.

### Evidence

- existing tests cover adaptive filtering, authenticity, labels, predictors,
  Brazilian escalation, COSMIC context, and reporting behavior;
- review layers escalate borderline findings to `REVIEW` rather than silently
  promoting them.

### Not fully tested in this round

- comprehensive per-filter stress matrix across all technical fields and real
  Brazilian cohorts.

## 10. Predictors

### Status
Partial approval.

### Evidence

- predictor tests are present and passing;
- predictor inventory is exported by the workflow;
- computational evidence is represented separately from filter status.

## 11. Authenticity

### Status
Partial approval.

### Evidence

- authenticity tests are present and passing;
- technical authenticity and biological support are separated conceptually and in
  implementation/documentation;
- the prior audit concern about hotspot/driver evidence masking technical
  weakness appears to have been addressed by the current code path.

## 12. Prioritization

### Status
Partial approval.

### Evidence

- prioritization tests and native maftools plot tests pass;
- workflow keeps prioritization downstream of filtering and authenticity.

## 13. COSMIC

### Status
Partial approval.

### Evidence

- same-build COSMIC preparation and matching are covered by tests;
- cache reuse and context-specific annotation behavior are exercised;
- `R CMD check` now passes with COSMIC-related tests active.

### Not tested

- cross-build liftover with real chain-file workflow in this round.

## 14. Brazilian Population

### Status
Partial approval.

### Evidence

- Brazilian population evidence tests pass;
- review escalation exists for possible Brazilian germline patterns.

### Not tested

- validation against curated Brazilian matched-normal cohorts.

## 15. Ancestry

### Status
Partial approval.

### Evidence

- ancestry tests pass locally;
- ancestry commands exist in the CLI and are now listed in usage text.

### Not tested

- external ancestry backends and large real-world cohort validation in this
  round.

## 16. OncoKB

### Status
Partial approval.

### Evidence

- `httr2` is now installed locally, removing a package-level environment blocker;
- OncoKB tests are present and passing.

### Not tested

- live API behavior with real token, rate limiting, and production error modes.

## 17. Plots

### Status
Approved for development use.

### Evidence

- `maftools::read.maf()` validation path is exercised;
- native `maftools` plotting path is implemented in `R/plots_maftools.R`;
- plot tests pass;
- no custom `ggplot2` imitation is used in the primary modern plotting path.

## 18. Report

### Status
Partial approval.

### Evidence

- report generation tests pass;
- fallback HTML builder writes output successfully in tests;
- workflow exports audit-oriented tables and plot manifests.

### Not fully audited in this round

- end-user usability review of the dashboard/report as a human triage interface.

## 19. ML

### Status
Not approved for promotion beyond experimental use.

### Evidence

- review/ML tests pass;
- architecture supports training registry and explicit activation flows.

### Limitation

- no new real-world validation evidence was generated in this round.

## 20. Performance

### Status
Not tested in this round.

### Evidence

- no fresh benchmark over 150k / 1M / 4M variant inputs was executed in this
  session.

## 21. Reproducibility

### Status
Improved materially.

### Evidence

- package tarball builds;
- full `R CMD check` is now green in the audited environment;
- package can be installed locally with tests included;
- installed package contains test files under:
  `/Users/marianaboroni/Documents/Variants_analysis/r_libs/tumoronly/tests/testthat`

## 22. Privacy

### Status
Not fully audited in this round.

### Evidence

- no credential leakage was observed during package install/check;
- deeper LGPD/privacy review was not rerun end-to-end here.

## 23. Tests

### Environment

- R: `4.1.2`
- OS: `Darwin 25.3.0 x86_64`
- local library: `/Users/marianaboroni/Documents/Variants_analysis/r_libs`

### Results confirmed in this round

- `R CMD build .` -> OK
- `R CMD check --no-manual tumoronly_0.1.0.tar.gz` -> **Status: OK**
- package installed locally with:
  `R CMD INSTALL --install-tests -l r_libs tumoronly_0.1.0.tar.gz`
- installed-package test tree confirmed present
- installed-package test execution completed successfully via:
  `testthat::test_dir(system.file("tests/testthat", package="tumoronly"), package="tumoronly", load_package="installed")`

### Test inventory

- `tests/testthat/test-*.R` files: **19**
- direct `test_that()` block count: **102**
- installed-package test runtime observed in this session: **~70 s**
- installed-package test summary observed in this session:
  - skips: **1** (cross-build error-path test became inapplicable once `rtracklayer` was installed)
  - warnings: **6** (all related to absent FASTA for REF validation in COSMIC prep tests)

### Notes

- a prior install without `--install-tests` caused `test_check("tumoronly")` to
  fail with `No test files found`; this is an installation-mode issue, not a
  code failure.

## 24. Limitations

This report should not be interpreted as:

- clinical validation;
- Brazilian matched-normal validation;
- final performance certification;
- complete external-integration certification for OncoKB or cross-build COSMIC.

## 25. Correction Plan

1. Add a reproducible audit helper or Makefile target that:
   - installs `Suggests`,
   - installs the package with `--install-tests`,
   - runs `R CMD build`,
   - runs `R CMD check`,
   - runs installed-package tests.
2. Converge legacy vs modern MAF/reporting code paths.
3. Execute a dedicated performance audit on large VCF/TSV inputs.
4. Audit chunked vs integral equivalence.
5. Validate Brazilian population and ancestry logic against curated real cohorts.
6. Run live OncoKB integration tests with a non-leaking audit token workflow.

## 26. Approval Matrix

| Area                 | Aprovada | Parcial | Reprovada | Evidência |
| -------------------- | -------: | ------: | --------: | --------- |
| Arquitetura          |          | X       |           | package real + auditorias prévias + sobreposição residual |
| CLI                  | X        |         |           | `doctor`, uso, comandos e ajuste de help |
| VCF                  |          | X       |           | testes locais presentes e passing |
| MAF                  | X        |         |           | `maftools::read.maf()` + testes |
| TSV                  |          | X       |           | ingestão multi-formato existe; matriz completa não reexecutada |
| Chunking             |          |         | X         | não testado nesta rodada |
| Filtros              |          | X       |           | testes passing; validação científica ainda parcial |
| Adaptativo           |          | X       |           | testes passing; ainda experimental |
| Preditores           |          | X       |           | testes passing; validação externa pendente |
| Autenticidade        |          | X       |           | testes passing; separação técnico/biológico implementada |
| Priorização          |          | X       |           | testes passing; refinamento futuro |
| COSMIC               |          | X       |           | same-build coberto; cross-build não reexecutado |
| População brasileira |          | X       |           | testes passing; coorte real pendente |
| Ancestralidade       |          | X       |           | testes passing; validação real pendente |
| OncoKB               |          | X       |           | `httr2` instalado; API live não auditada |
| MAFtools             | X        |         |           | pipeline nativo confirmado |
| Relatório            |          | X       |           | geração OK; usabilidade não reavaliada por completo |
| Revisão humana       |          | X       |           | fluxo existe; governança real ainda parcial |
| ML                   |          | X       |           | arquitetura/testes OK; pronto só como experimental |
| Performance          |          |         | X         | não testado nesta rodada |
| Reprodutibilidade    | X        |         |           | build/check/install-tests OK |
| Privacidade          |          | X       |           | sem vazamento observado; auditoria profunda não refeita |
