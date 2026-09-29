# Methods, assumptions and limitations

## Processing allocation

For each stage, subtract configured stage-control counts from every batch's total capacity. Form indivisible bundles by experimental unit when keep_unit_together is true; otherwise each processing sample is an allocation unit. Generated controls are individual allocation units. Fixed/permitted batch restrictions are intersected within a bundle.

Each random restart shuffles bundles, then prioritizes the most constrained and largest bundles. A greedy step chooses a feasible batch minimizing the incremental sum of squared deviations from capacity-proportional categorical target counts, with a small occupancy penalty. Randomized ties provide alternative allocations. The lowest categorical-imbalance score across successful restarts is retained. Rows are randomized within each selected batch to assign processing positions. Different configurations may consume the random stream differently even with the same seed.

The objective balances sample counts, not equal numbers of independent units in each batch or all factor interactions. Repeated samples can dominate the objective. The separate independent-unit table exposes total replication. Balancing is marginal; users must review interactions and continuous covariates. This heuristic does not prove optimality or feasibility when a search fails. Hard constraints take precedence over balance.

Extraction is planned first, including generated extraction blanks and mock communities. Its complete output becomes the sequencing input; sequencing-stage blanks/positives are then added. Extraction controls are never silently omitted from sequencing capacity or cost. Sequencing balances source extraction batches alongside requested categorical factors. Library-positive controls are conceptual planned material and are not generated physically by the software.

## Confounding checks

Categorical tables report sample counts, descriptive Cramer's V and the number of connected components in the observed factor–batch incidence graph. More than one component means some additive factor contrasts cannot be separated from batch. V>0.5 prompts a descriptive imbalance warning, not a significance claim or validated threshold. A connected graph does not establish good precision, adequate replication or lack of higher-order confounding.

Checks are repeated between group and other requested balancing factors, to expose pre-existing group/site/time imbalance that processing randomization cannot repair. Only recorded factors are checked. The software does not infer randomization, exchangeability or causal identification from balanced tables.

Independent replication is counted using the declared experimental unit. Group membership must be constant within it; fewer than two units in a group is flagged. More than two units does not establish adequate power. Distinct specimen counts require correct specimen_id annotations, and multiple processing samples from one specimen are not independent biological samples.

## Costs and attrition

Budget totals multiply declared unit/specimen/extraction/library counts by explicitly supplied rates. Extraction controls count in both laboratory stages; library controls count only in sequencing. Fixed fees, storage, repeats and collection-control costs are not separately modeled. A provided rate must include any relevant costs if users intend a complete budget.

Attrition simulations independently remove entire experimental units at the specified rate, then independently fail processing samples among surviving units. Removing a unit removes all its samples. Remaining units have at least one surviving sample. A complete schedule requires at least one surviving processing sample for every original specimen of that unit; it does not require every technical replicate. These are simplified scenarios, not estimators of real dropout behavior. Time-varying, batch-wide, informative or group-specific failure are not modeled. Percentiles summarize Monte Carlo retention, not confidence intervals or hypothesis-test power.

## Scope and next validation

No pilot abundance tables, power estimates, selected sample-size recommendations, biological effect simulations or tests are produced. Endpoint-specific text suggests issues for a future analysis plan and must be adapted to the actual design. Validate any future power module against appropriate established implementations and simulation scenarios before making scientific planning claims.

## Background sources

These sources motivate design considerations; they do not validate this implementation or prescribe its control counts, balancing score or thresholds.

- [A guide to human microbiome research: study design, sample collection, and bioinformatics analysis](https://pubmed.ncbi.nlm.nih.gov/32604176/).
- [Conducting a Microbiome Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5074386/).
- [Guidelines for preventing and reporting contamination in low-biomass microbiome studies](https://www.nature.com/articles/s41564-025-02035-2).
- [The rise to power of the microbiome: power and sample size calculation for microbiome studies](https://www.nature.com/articles/s41385-022-00548-1).

MicroPlan complements existing design and power software; it makes no claim to be the first microbiome planning tool.
