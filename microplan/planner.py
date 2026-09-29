import json
import platform
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from . import __version__
from .core import validate,allocate,audit,loss_scenarios,budget_summary,write_tsv,sha


def run(config_path,out):
    config_path=Path(config_path).resolve();out=Path(out).resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):raise ValueError('Output directory must be new or empty')
    cfg,rows=validate(json.loads(config_path.read_text()),config_path.parent)
    source=(config_path.parent/cfg['samples']).resolve();unit=cfg['experimental_unit']
    warnings=[];unit_counts=[]
    for group in sorted({r['group'] for r in rows}):
        selected=[r for r in rows if r['group']==group]
        n=len({r[unit] for r in selected});specimens=len({r['specimen_id'] for r in selected})
        unit_counts.append(dict(group=group,n_independent_units=n,n_specimens=specimens,n_processing_samples=len(selected)))
        if n<2:warnings.append(f'{group}: fewer than two declared independent experimental units. Technical replicates or repeated samples do not supply independent replication.')
    if len(unit_counts)<2:warnings.append('Only one study group: a between-group contrast cannot be assessed.')
    if len({r[unit] for r in rows})<len(rows):warnings.append(f'Repeated observations share {unit}; downstream tests must account for that dependence. This tool trusts your declared experimental unit.')
    if len({r['specimen_id'] for r in rows})<len(rows):warnings.append('Multiple processing samples represent the same specimen. They are technical replicates, not extra biological specimens.')
    for stage in ('extraction','sequencing'):
        controls=cfg['stages'][stage].get('controls_per_batch',{})
        if stage=='extraction':
            if not controls.get('extraction_blank'):warnings.append('No extraction blanks reserved. Include appropriate blanks, especially for low-biomass samples.')
            if not controls.get('mock_community'):warnings.append('No mock-community extraction controls reserved; consider suitable positive controls for the assay.')
        if stage=='sequencing' and not controls.get('library_blank'):warnings.append('No library/preparation blanks reserved. For amplicon workflows, consider appropriate no-template PCR controls.')
    if cfg['assay']=='shotgun':warnings.append('Shotgun planning must also consider host DNA, biomass and usable microbial reads. This version does not optimize sequencing depth.')
    if cfg['assay'] in ('16s','its'):warnings.append('Record primer region, amplification protocol and reference database. This tool does not assess primer coverage or amplification bias.')
    extraction,score1=allocate(rows,'extraction',cfg['stages']['extraction'],cfg['balance_columns'],unit,cfg['seed'],cfg['search_restarts'])
    sheet,score2=allocate(extraction,'sequencing',cfg['stages']['sequencing'],cfg['balance_columns'],unit,cfg['seed']+1,cfg['search_restarts'])
    balances=[]
    for stage in ('extraction','sequencing'):
        table,issues=audit(sheet,stage,cfg['balance_columns'],unit);balances+=table;warnings+=issues
    other=[c for c in cfg['balance_columns'] if c!='group']
    if other:
        study=[dict(r,study_group_batch=r['group']) for r in rows]
        _,issues=audit(study,'study_group',other,unit)
        warnings.extend('Study group comparison: '+issue.replace('study_group: ','').replace('batch','group') for issue in issues if 'only one occupied' not in issue)
    budget=budget_summary(rows,sheet,cfg)
    if budget['status']=='over_budget':warnings.append(f"Estimated cost {budget['estimated_total']:.2f} exceeds the supplied budget {budget['limit']:.2f} {budget['currency']}.")
    losses=loss_scenarios(rows,cfg)
    warnings.append('Allocation balances recorded categorical margins, not every interaction or unmeasured confounder. A feasible plan is not evidence of adequate power.')
    warnings.append('Collection schedule, plate geometry, forbidden wells, indexing and reagent compatibility must be reviewed before laboratory use; numbered positions are generic processing positions.')
    out.mkdir(parents=True,exist_ok=True)
    fields=['sample_id','sample_type','source_stage','specimen_id',unit,'group']
    fields=list(dict.fromkeys(fields+[k for row in rows for k in row]+['extraction_batch','extraction_position','sequencing_batch','sequencing_position']))
    sheet.sort(key=lambda r:(r['sequencing_batch'],r['sequencing_position']))
    write_tsv(out/'sample_sheet.tsv',sheet,fields)
    write_tsv(out/'extraction_sheet.tsv',sorted((r for r in sheet if r.get('extraction_batch')),key=lambda r:(r['extraction_batch'],r['extraction_position'])),fields)
    write_tsv(out/'sequencing_sheet.tsv',sheet,fields)
    write_tsv(out/'balance.tsv',balances,['stage','factor','batch','level','n_samples','cramers_v','overlap_components'])
    write_tsv(out/'independent_units.tsv',unit_counts,['group','n_independent_units','n_specimens','n_processing_samples'])
    write_tsv(out/'loss_scenarios.tsv',losses,list(losses[0]))
    utilization=[]
    for stage in ('extraction','sequencing'):
        for b in cfg['stages'][stage]['batches']:
            occupied=[r for r in sheet if r.get(stage+'_batch')==b['id']]
            utilization.append(dict(stage=stage,batch=b['id'],capacity=b['capacity'],used=len(occupied),biological=sum(r['sample_type']=='biological' for r in occupied),controls=sum(r['sample_type']!='biological' for r in occupied),unused=b['capacity']-len(occupied)))
    write_tsv(out/'capacity.tsv',utilization,list(utilization[0]))
    (out/'budget.json').write_text(json.dumps(budget,indent=2)+'\n')
    (out/'resolved_config.json').write_text(json.dumps(dict(cfg,samples=str(source)),indent=2)+'\n')
    manifest=dict(tool='MicroPlan',version=__version__,created_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        config=cfg,inputs=[dict(path=str(p),sha256=sha(p)) for p in (config_path,source)],
        search=dict(algorithm='Random-restart constrained greedy allocation; randomized ties; categorical marginal imbalance objective',restarts=cfg['search_restarts'],extraction_score=score1,sequencing_score=score2),
        warnings=warnings,power_analysis='not_performed',treatment_randomization='not_performed',status='completed_with_warnings' if warnings else 'completed')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    report=['# Microbiome study-design report','',cfg.get('study_name','Untitled microbiome study'),'',
        f"Assay: {cfg['assay']}. Intended endpoint: {cfg['endpoint']}. Experimental unit: {unit}.",'',
        f"Allocated {len(rows)} biological processing samples and {len(sheet)-len(rows)} controls across extraction and sequencing batches.",
        'All extraction controls are carried through to sequencing. Sequencing-stage controls do not consume extraction positions.','',
        '## Independent replication','']
    for row in unit_counts:report.append(f"- {row['group']}: {row['n_independent_units']} independent units; {row['n_specimens']} specimens; {row['n_processing_samples']} processing samples.")
    report+=['','## Allocation and budget','',
        'Treatments and collection times were not assigned or changed. Processing batches and positions were assigned under the supplied capacity/restriction rules. '
        'The heuristic attempts to reduce categorical imbalance; it does not prove optimality.','',f"Budget status: {budget['status']}."]
    if 'estimated_total' in budget:report.append(f"Estimated total: {budget['estimated_total']:.2f} {budget['currency']}; supplied limit: {budget['limit']:.2f}.")
    report+=['','## Attrition scenarios','',
        'loss_scenarios.tsv simulates whole-unit loss plus independent processing-sample failure using the declared rates. '
        'Its percentiles describe simulated retention, not confidence intervals or power. It does not model time-dependent dropout or group-dependent loss.','',
        '## Checks and limitations','']+['- '+w for w in warnings]
    report+=['','## Next steps','',
        'Review every flagged confounding condition and the lab-specific constraints. Define the primary contrast, meaningful effect and analysis model before estimating power. '
        'Do not interpret an allocation file or a within-budget plan as approval of sample size. See ANALYSIS_PLAN.md for endpoint-specific considerations.','']
    (out/'REPORT.md').write_text('\n'.join(report))
    endpoint={
        'community':'Prespecify a community-distance metric and group/time contrasts. Consider dispersion as well as centroid differences. Permutations must respect the declared experimental units; simple row permutations are inappropriate for repeated or clustered observations.',
        'diversity':'Prespecify the diversity measure and its interpretation under sampling depth. Use a model or contrast appropriate to independent, paired or repeated measurements; avoid treating technical replicates as independent.',
        'taxon':'Prespecify taxon/feature filtering, compositional assumptions, abundance model, covariates and multiplicity correction. An assay-appropriate differential-abundance method is required; relative-abundance differences are not automatically absolute changes.',
        'exploratory':'Define exploratory questions and reserve independent validation where possible. Do not select a primary endpoint after observing which analysis is significant.'}[cfg['endpoint']]
    (out/'ANALYSIS_PLAN.md').write_text('# Draft analysis considerations\n\n'+endpoint+'\n\nInclude relevant batch, collection-site and host/environmental covariates where estimable. Complete confounding cannot be repaired by adding a batch covariate. Review storage duration, biomass, extraction kit/lot, collection timing and assay-specific controls.\n\nThis is guidance, not a fitted model or a completed statistical analysis plan. Power, effect sizes, sequencing depth and optimal sampling times were not calculated.\n')
    (out/'README.md').write_text('''# What this run did

Validated input identities and capacities; counted declared independent units and specimens; assigned extraction and sequencing batches with randomized constrained greedy search; reserved controls; carried extraction controls forward; checked categorical batch overlap; estimated supplied costs; simulated operational sample losses. No abundance data were analysed, no treatment allocation was changed, no power calculation or graphs were generated, and no network calls were made.

- REPORT.md: findings and limitations.
- sample_sheet.tsv: every biological/control row with both processing assignments.
- extraction_sheet.tsv and sequencing_sheet.tsv: stage-specific work lists in processing order. Position numbers are not instrument well labels or indexes.
- independent_units.tsv: independent-unit, specimen and processing counts by group.
- balance.tsv: categorical sample counts, descriptive Cramer's V and overlap components; not p-values.
- capacity.tsv: used, unused and control positions per batch.
- loss_scenarios.tsv: simulated remaining samples, units and complete specimen schedules; not power.
- budget.json: explicit-cost estimate, including controls, or a not-assessed status.
- ANALYSIS_PLAN.md: endpoint-specific considerations requiring study-specific review.
- resolved_config.json and manifest.json: configuration, seed, search settings, source hashes, software version and warnings.

Read identifiers as text. Empty group/unit/specimen fields on generated control rows mean not applicable. Existing subject/group/timepoint labels were preserved. Technical replicates require a shared specimen_id; otherwise each row is assumed to represent a distinct specimen. The experimental unit is supplied by the researcher and cannot be inferred reliably from this table.

The same version, configuration, row order and seed reproduce the allocations. Processing randomization helps balance known factors but cannot remove biological confounding. The output contains sensitive metadata if supplied; review it before sharing. COMPLETE means file generation finished, not that the design is statistically adequate.
''')
    (out/'COMPLETE').write_text('Design files generated. Review REPORT.md before use. No power calculation performed.\n')
    return dict(biological_samples=len(rows),controls=len(sheet)-len(rows),independent_units=sum(r['n_independent_units'] for r in unit_counts),budget_status=budget['status'],output=str(out))
