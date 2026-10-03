# NHANES Phthalates–Stroke Reassessment — v1.1.0

This is a local release candidate. VERSION = 1.1.0; TAG_DRAFT = v1.1.0.
PUBLIC_RELEASE_READY = YES; PUBLIC_RELEASE_EXECUTED = NO.
No official v1.1.0 tag, release date or version DOI has been assigned.

v1.1.0 supersedes v1.0.0 / v1.0.1 for reproduction of the current manuscript results. Earlier public source files remain byte-identical under `history/public-v1.0.1`; registered execution scripts and the isolated second mixed driver remain under `history`. Historical code and tests are not entrypoints. Do not run them as the current analysis.

## Inputs and runtime

Use Python 3.10+ (standard library only), R 4.5.2 and the recorded package versions in `environment/package-versions.csv`. Install/configure these libraries before execution using normal R library configuration, such as R_LIBS_USER. No script installs packages or embeds a user library path.

The reconstruction requires 66 translated official NHANES components, 16 official phthalate/alcohol XPTs and 24 official codebooks. Exact component URLs and official XPT identities are in `config/nhanes-components.csv`; codebook URLs and hashes are in `config/official-codebook-urls.csv`. `config/source-inputs.csv` records registered input identities. The independently acquired RDS objects matched the registered objects exactly; alternate serialization hashes are recorded separately in `config/validated-acquisition-sha256.csv`. Both serialization identities are verified input routes and do not change values, factor levels or participant membership. All inputs and reconstructed participant-level data remain outside this package.

From any working directory, with explicit paths:

```text
python <package_root>/scripts/acquire_inputs.py --package-root <package_root> --input-root <empty_input_root> --source-cache <official_source_cache> --rscript <Rscript>
python <package_root>/scripts/run_pipeline.py --package-root <package_root> --input-root <input_root> --work-root <empty_external_workspace> --rscript <Rscript>
```

An optional `--translated-cache <cache_root>` uses a separately held registered cache whose hashes/dimensions are verified before import. Without it, official XPTs are translated with nhanesA. A source cache may already hold the registered official XPTs; otherwise they are downloaded and verified. The current runner refuses a nonempty analysis workspace and never imports historical models, partial estimates or completed imputations.

## Current execution and expected outputs

The runner verifies sources/alcohol coding, reconstructs all component joins, applies the registered cholesterol proxy, closes participant selection, and fits primary Models 1–3. It then runs common-cycle and categorical cycle sensitivities, sex interactions, RCS, grouped molar sums, male subgroup global interactions, era analyses, Table 1, fresh MI FCS and MI Model 3, ordinary WQS and qgcomp. Each step writes explicit external outputs and logs. Expected canonical CSVs under `expected` are comparison references, never computed results.

Core counts are 10,528/385 stroke events (female 5,222/198; male 5,306/187); common-cycle counts are 9,366/348 (female 4,642/181; male 4,724/167). The first eight metabolites use the 2003–2018 framework; MCNP/MCOP use 2005–2018. The primary Model 3 BH family contains 20 tests with zero q < .05 survivors. Definitions, covariates, survey weights and inference families remain those of the registered corrected analysis.

MI uses four separate sex/window frameworks, m = 50, maxit = 20, seed = 20260925, the registered predictor matrices/methods and nnet.maxit = 1000. Each framework restores its recorded initialization RNG state and stops on automatic predictor changes, logged events or convergence failures. Survey MI pooling retains the registered finite complete-data df rule. Ordinary WQS uses seed = 2026, q = 4, validation = 0, b = 1000 and a sequential gWQS plan; qgcomp bootstrap uses seed = 2026, q = 4 and B = 1000. RH-WQS is withdrawn and no active entrypoint executes it.

Comparison policy was fixed before R execution: absolute difference <= 1e-10 + 1e-8 * abs(expected). Integer sample counts/identities and strings are exact; missing values must agree. Failures stop execution; the tolerance is never relaxed. The `config/lod_exact_comparison.csv` input preserves source-audit annotations; the analysis uses official released exposure values without applying a second LOD substitution.

CLEAN_RUN_VALIDATION = PASS. Fresh aggregate outputs are included only following the full validation gate. Participant data, SEQN lists, completed MI/model objects and credentials are not redistributed. MIT license is retained. Publication remains a separately authorized action.
