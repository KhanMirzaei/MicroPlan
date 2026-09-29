# MicroPlan

Plan microbiome sample processing before laboratory work begins. MicroPlan produces randomized, batch-balanced sample sheets, checks the declared study design for confounding and insufficient independent replication, and documents capacity, costs and sample-loss scenarios.

**Version 0.1.0 is a research prototype.** Name availability has not been checked. It supports 16S, ITS and shotgun study-planning contexts. It does not analyse abundance tables or estimate statistical power.

## Install

Python 3.10 or newer is required. There are no runtime dependencies beyond the Python standard library.

From the extracted MicroPlan directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
microplan --help
```

On Windows, activate `.venv\Scripts\activate`. Alternatively install the supplied wheel with `python -m pip install /path/to/microplan-0.1.0-py3-none-any.whl`. The package is not published to PyPI. Source installation can fetch its build dependency; the supplied wheel requires no runtime dependency downloads. The tool itself makes no network calls.

## Quick start

```bash
microplan plan --config examples/study.json --out results
```

The included example is entirely synthetic: 12 subjects in two groups, each with baseline and follow-up specimens. It illustrates allocation mechanics, not an adequate sample size recommendation. Edit the sample table and configuration for the real study. Runs are noninteractive; invalid input exits with an error rather than asking questions.

## What to provide

1. A CSV or TSV sample manifest with one row per planned processing sample.
2. A JSON configuration describing the experimental unit, assay, endpoint, processing capacities and balancing factors.

The experiment's biological groups and time points must already be defined. MicroPlan randomizes **processing allocation**, not patient treatment or collection times.

Required sample columns:

| Column | Meaning |
|---|---|
| sample_id | Unique processing-sample identifier, including separate IDs for technical replicates |
| group | Planned biological group; must be constant within the declared experimental unit |
| Your experimental-unit column | For example subject_id, animal_id, cage_id or independently sampled plot_id |
| Every balance column | Categorical factors to distribute across batches, e.g. group, timepoint, site |

Optional sample columns:

- `specimen_id`: shared by technical replicates from one collected specimen. If omitted, each processing row is assumed to represent a different specimen.
- `fixed_extraction_batch`, `fixed_sequencing_batch`: a required batch ID for that sample.
- `allowed_extraction_batches`, `allowed_sequencing_batches`: semicolon-separated permitted batch IDs, e.g. `E1;E2`. Blank means unrestricted.
- Other metadata are preserved. Blank required values, duplicate IDs and inconsistent specimen metadata are rejected. Identifiers remain text, including leading zeros.

All input rows represent biological processing samples. Stage controls are generated from configuration. Generated IDs start with `CONTROL__`; this prefix is reserved. Collection/field blanks and other custom control workflows need separate laboratory planning in this release.

The **experimental unit is not automatically inferred**. If treatment is assigned to a cage or plot, declaring each animal or aliquot independent can misrepresent replication. Repeated samples and technical replicates do not increase the number of independent units. Crossover/time-varying treatment groups within a unit are not supported.

## Configure the design

See [the complete example](examples/study.json). Required settings:

- `samples`: CSV/TSV path, relative to the configuration file unless absolute.
- `experimental_unit`: column defining independent units.
- `assay`: `16s`, `its`, or `shotgun`.
- `endpoint`: `community`, `diversity`, `taxon`, or `exploratory`—used for analysis guidance, not power calculation.
- `stages`: both extraction and sequencing, each with named batches and total capacities.

Optional settings:

- `balance_columns`: categorical columns; defaults to `["group"]` and must include group. Continuous covariates should be represented by scientifically justified strata if balancing is needed. The tool does not choose cutoffs.
- `seed`: default 42; `search_restarts`: default 100, at most 10,000.
- `keep_unit_together`: per-stage flag keeping an experimental unit's biological samples in one batch. This is a design choice, not a universal recommendation. Batch constraints may make it impossible.
- `controls_per_batch`: explicit counts. Extraction supports `extraction_blank` and `mock_community`; sequencing supports `library_blank` and `library_positive`. No universal count is assumed. Omitted controls produce reminders.
- `budget`: currency, total limit and nonnegative costs per independent unit, specimen, extraction and library. All values must be supplied when assessing cost. Controls count toward processing cost.
- `loss_scenarios`: whole-unit loss rates, independent processing-sample failure rate and simulation count. Defaults are 0%, 10%, 20% unit loss, 0% sample failure and 500 simulations; these are illustrative scenarios, not estimated rates.

Extraction controls are allocated to sequencing too. Additional sequencing controls consume only sequencing capacity. Sequencing allocation also attempts to distribute incoming extraction batches. Unused positions are reported; positions are generic integers, not instrument-specific well labels or barcode assignments. No physical collection schedule or storage constraint is inferred—encode batch eligibility explicitly.

## Outputs

- **REPORT.md:** design findings, limitations and next steps.
- **sample_sheet.tsv:** all samples and controls with extraction/sequencing batch and position.
- **extraction_sheet.tsv / sequencing_sheet.tsv:** stage-specific ordered work lists.
- **independent_units.tsv:** independent units, distinct specimens and processing samples per group.
- **balance.tsv / capacity.tsv:** recorded-factor balance, overlap and utilization.
- **loss_scenarios.tsv:** simulated remaining samples, units and complete specimen schedules.
- **budget.json:** estimated cost including controls, or a not-assessed status.
- **ANALYSIS_PLAN.md:** endpoint-specific draft considerations.
- **README.md:** what the run did and how to read outputs.
- **manifest.json / resolved_config.json:** settings, seed, source hashes and warnings.

No browser, server or graphs. The output directory must be new or empty. `COMPLETE` means generation finished; it does not certify study adequacy. A plan with fixed confounding or excess budget is returned with explicit warnings for review. Infeasible capacities or contradictory restrictions fail before results are written.

## What it does not do

This first version does not determine the best number of participants, optimize sampling times or sequencing depth, estimate microbiome power from pilot tables, assign treatments, generate a blinded codebook, fit microbiome models, handle every laboratory restriction, or guarantee a globally optimal allocation. Its attrition simulation is operational planning, **not a statistical power simulation**.

A future power module must be validated for a specified endpoint and experimental-unit structure. Community differences, diversity and individual taxa cannot be assigned a single universal microbiome sample-size calculation. See [methods and references](docs/METHODS.md) and [validation status](docs/VALIDATION.md).

## HPC and development

The CLI can run in an ordinary batch job without a web port, GPU or external account. [Slurm example](examples/slurm.sh). It has not been tested on an actual cluster. The allocation search is serial; repeated restarts, large manifests and many constraints increase runtime. Input is limited to 10,000 biological processing rows; this is a guard, not a scaling benchmark.

```bash
python -m unittest discover -s tests -v
```

Keep generated results and real sample metadata outside the source repository. MIT license. Author: Mohammadali Khan Mirzaei.
