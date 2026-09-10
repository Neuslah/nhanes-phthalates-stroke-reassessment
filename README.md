# NHANES urinary phthalates and prevalent stroke: local reproducibility package

## Version history

**v1.0.0** — Initial public reproducibility release, published 21 August 2026.

**v1.0.1** — Minimal reproducibility correction: adds a portable aggregate reproduction route for the four corrected Supplementary Table S3 variables and documents their variable-specific non-missing denominators. It does not change the analytic sample, primary models, FDR family, secondary/exploratory model definitions, or primary statistical results.

The dataset pipeline emits `dataset/s3_four_variable_summary.csv` under the user-specified reproduction output directory. It compares the analytic group with the urinary-creatinine-unavailable group for hypertension, diabetes, hyperlipidemia, and ever smoking. `prevalence` is a proportion (0–1); SMD is the absolute binary standardized mean difference using each variable's non-missing denominator. `config/s3-four-variable-reference.csv` supplies the final aggregate validation anchors. This route does **not** claim to reproduce every cell of Supplementary Table S3. No participant-level S3 file is saved or released.

## Study scope

This package reproduces the analysis of NHANES 2003–2018 data for cross-sectional associations between ten urinary phthalate metabolites and prevalent self-reported stroke among adults aged 20 years or older. It contains code, configuration, aggregate verification anchors, and environment records. It does not contain participant-level data.

## Data source

All source components are official CDC/NCHS NHANES public-use files. `config/nhanes-components.csv` provides the exact component URLs, cycle-specific cache filenames, and SHA-256 identities used by the reconstruction. `R/01_download_or_import_nhanes.R` downloads a missing source file or verifies an existing cache file before translating it. For isolated environments without online codebook access, a separately held translated official-source cache may be supplied and is verified against `config/translated-cache-manifest.csv`; no translated data are redistributed. The package does not redistribute raw XPT files or any derived analytic RDS.

## Analysis hierarchy

- Primary: sex-stratified survey-weighted single-metabolite Models 1–3. The primary multiplicity family is the 20 fully adjusted Model 3 tests (10 metabolites × 2 sexes), controlled using Benjamini–Hochberg FDR.
- Secondary: sex-interaction tests and the 2005–2018 common-cycle sensitivity analysis.
- Exploratory: restricted cubic spline analyses and grouped molar-sum analyses. Era-stratified analyses are descriptive/exploratory. Male MCNP/MCOP subgroup global-interaction tests are non-confirmatory robustness analyses outside the primary FDR family.
- Secondary-exploratory mixture: WQS, qgcomp, and repeated-holdout WQS are unweighted. They do not support complex-survey national inference.

## Survey design

The deterministic survey analyses use cycle-appropriate phthalate subsample weights: `WTSA2YR` for 2011–2012 and `WTSB2YR` for the other cycles. The eight-cycle weight is the selected two-year subsample weight divided by 8. Survey designs use `SDMVSTRA`, `SDMVPSU`, nesting, and Taylor linearization through the R `survey` package.

## Restricted cubic splines

The current RCS implementation uses exactly three explicit knots at the ordinary empirical 10th, 50th, and 90th percentiles (`quantile(..., type = 7)`) within the corresponding sex-specific model-complete sample. It does not use package-default knot selection.

## Execution order

Use R 4.5.2 or a compatible R 4.5.x environment with the package versions recorded under `environment/`. Supply explicit directories; no script depends on a particular drive, username, or working directory.

1. Create empty directories for an official-source cache, a reproduction workspace, and aggregate outputs.
2. Run the complete deterministic pipeline:

   ```text
   Rscript --vanilla R/00_run_deterministic_pipeline.R <package_root> <official_source_cache> <translated_cache_or_dash> <workspace_root> <output_root>
   ```

   This performs official-source identity checks/import, dataset reconstruction, primary models, primary FDR, deterministic secondary/exploratory analyses, the current grouped and RCS routes, and lightweight mixture validation.
3. Run the package tests:

   ```text
   Rscript --vanilla tests/run_tests.R <package_root> <workspace_root> <output_root>
   ```

   When full mixture outputs are available, `tests/03_inherited_mixture_output_tests.R` checks the registered WQS, qgcomp, and repeated-holdout aggregate anchors and schemas. Stage 5-C used this test with the previously verified full-run outputs instead of rerunning the hours-long models.

4. To fit the full mixture models, run scripts 11–13 with their final argument set to `run`, in numerical order. Full WQS must precede repeated-holdout WQS. The fixed settings are embedded and fail-fast; they are not user-selectable scientific options.

Repeated-holdout partition identity is verified from the serialized R object content. Regenerated compressed RDS file bytes may differ because compression metadata are not a scientific partition identifier; the canonical file SHA-256 is retained separately as provenance.

Intermediate translated and analytic RDS files are written only to the user-declared reproduction workspace. Do not place the workspace inside a public repository.

## Reproducibility anchors

The deterministic reconstruction stops if any of these anchors differs:

| Anchor | Expected |
|---|---:|
| Final analytic sample | 9,129 |
| Stroke events | 375 |
| Female | 4,611 |
| Male | 4,518 |
| Common-cycle sample / events | 8,311 / 339 |
| Mixture Female sample / events | 4,179 / 175 |
| Mixture Male sample / events | 4,132 / 164 |

The participant-flow anchor file also verifies all preceding selection counts.

## Package contents

- `R/`: portable analysis scripts and one deterministic runner.
- `config/nhanes-components.csv`: official-source component manifest.
- `config/translated-cache-manifest.csv`: identity and dimension gate for an optional separately held translated official-source cache.
- `config/variable-dictionary.csv`: source-to-derived variable dictionary.
- `environment/`: reproduction-environment records; these describe the successful reproduction environment and are not presented as the unknown historical-original environment.
- `tests/`: static, anchor, schema, survey-weight, current grouped/RCS, mixture-sample, and repeated-holdout partition checks.
- `outputs/`: aggregate example/import audit only; runtime outputs belong in the user-declared output directory.

Historical grouped and RCS implementations were superseded during reproducibility QC and are not active execution routes in this package.

## Limitations

This code reproduces the analysis; it does not convert cross-sectional associations into causal evidence. The outcome is self-reported prevalent stroke. The mixture analyses are deliberately unweighted and exploratory and do not support design-based national estimates.

## License

This reproducibility package is released under the MIT License. See `LICENSE`.
