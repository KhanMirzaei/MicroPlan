"""Auditable processing allocation; never assigns treatments or estimates power."""
import csv
import hashlib
import json
import math
import random
import statistics
from collections import Counter,defaultdict
from pathlib import Path

STAGES=('extraction','sequencing')
RESERVED={'sample_type','source_stage','extraction_batch','sequencing_batch','extraction_position','sequencing_position'}


def integer(value,label,minimum=0):
    if type(value) is not int or value<minimum:raise ValueError(f'{label} must be an integer >= {minimum}')
    return value


def number(value,label):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
        raise ValueError(f'{label} must be a finite nonnegative number')
    return value


def read_samples(path,unit,columns):
    path=Path(path)
    with path.open(encoding='utf-8-sig',newline='') as handle:
        reader=csv.DictReader(handle,delimiter=',' if path.suffix.lower()=='.csv' else '\t')
        header=reader.fieldnames or []
        if len(set(header))!=len(header) or any(not x for x in header):raise ValueError('Duplicate or blank sample-table columns')
        if RESERVED&set(header):raise ValueError('Input uses reserved output columns: '+', '.join(sorted(RESERVED&set(header))))
        required={'sample_id','group',unit,*columns}
        if not required<=set(header):raise ValueError('Missing columns: '+', '.join(sorted(required-set(header))))
        rows=[]
        for row in reader:
            if None in row or any(v is None for v in row.values()):raise ValueError('Ragged sample table')
            if any(not row[c].strip() for c in required):raise ValueError('Blank identifier/group/balancing value')
            if any(row[c]!=row[c].strip() for c in required):raise ValueError('Identifier and balancing values must not have surrounding whitespace')
            row=dict(row,sample_type='biological',source_stage='input')
            row.setdefault('specimen_id',row['sample_id'])
            if not row['specimen_id'].strip():raise ValueError('Blank specimen_id')
            if row['sample_id'].startswith('CONTROL__'):raise ValueError('CONTROL__ sample prefix is reserved')
            rows.append(row)
    if not rows:raise ValueError('No biological samples supplied')
    if len({r['sample_id'] for r in rows})!=len(rows):raise ValueError('Duplicate sample_id')
    membership=defaultdict(set);specimens=defaultdict(set)
    for row in rows:
        membership[row[unit]].add(row['group'])
        specimens[row['specimen_id']].add(tuple(row.get(c,'') for c in [unit,'group','timepoint']))
    if any(len(x)>1 for x in membership.values()):raise ValueError('An experimental unit belongs to multiple groups; crossover/time-varying treatment designs are not supported')
    if any(len(x)>1 for x in specimens.values()):raise ValueError('Replicates of a specimen have inconsistent unit/group/timepoint metadata')
    return rows


def validate(cfg,base):
    allowed={'samples','experimental_unit','assay','endpoint','balance_columns','stages','seed','search_restarts','loss_scenarios','budget','study_name'}
    if set(cfg)-allowed:raise ValueError('Unknown configuration settings: '+', '.join(sorted(set(cfg)-allowed)))
    for key in ['samples','experimental_unit','assay','endpoint','stages']:
        if key not in cfg:raise ValueError(f'Missing configuration setting: {key}')
    if cfg['assay'] not in ('16s','its','shotgun'):raise ValueError('assay must be 16s, its, or shotgun')
    if cfg['endpoint'] not in ('community','diversity','taxon','exploratory'):raise ValueError('Unsupported endpoint')
    cols=cfg.get('balance_columns',['group'])
    if not isinstance(cols,list) or not all(isinstance(c,str) and c for c in cols):raise ValueError('balance_columns must be a list of column names')
    if 'group' not in cols:raise ValueError('balance_columns must include group')
    if len(set(cols))!=len(cols):raise ValueError('Duplicate balance column')
    cfg=dict(cfg,balance_columns=cols,seed=integer(cfg.get('seed',42),'seed'),search_restarts=integer(cfg.get('search_restarts',100),'search_restarts',1))
    if cfg['search_restarts']>10000:raise ValueError('search_restarts must be <=10000')
    if set(cfg['stages'])!=set(STAGES):raise ValueError('Configure exactly extraction and sequencing stages')
    for stage in STAGES:
        spec=cfg['stages'][stage]
        if set(spec)-{'batches','controls_per_batch','keep_unit_together'}:raise ValueError(f'{stage}: unknown stage setting')
        if not isinstance(spec.get('batches'),list) or not spec['batches']:raise ValueError(f'{stage}: specify batches')
        names=[]
        for b in spec['batches']:
            if set(b)!={'id','capacity'}:raise ValueError(f'{stage}: each batch needs id and capacity')
            if not isinstance(b['id'],str) or not b['id'] or any(not(c.isalnum() or c in '_-') for c in b['id']):raise ValueError('Batch IDs must be simple identifiers')
            integer(b['capacity'],'capacity',1);names.append(b['id'])
        if len(names)!=len(set(names)):raise ValueError('Duplicate batch IDs')
        controls=spec.get('controls_per_batch',{})
        valid={'extraction_blank','mock_community'} if stage=='extraction' else {'library_blank','library_positive'}
        if set(controls)-valid:raise ValueError(f'{stage}: unsupported control type')
        for kind,count in controls.items():integer(count,kind)
        if not isinstance(spec.get('keep_unit_together',False),bool):raise ValueError('keep_unit_together must be true/false')
        if any(sum(controls.values())>=b['capacity'] for b in spec['batches']):raise ValueError(f'{stage}: controls consume batch capacity')
    rows=read_samples(base/cfg['samples'],cfg['experimental_unit'],cols)
    if len(rows)>10000:raise ValueError('Prototype supports <=10000 biological rows')
    for r in rows:
        for stage in STAGES:
            names={b['id'] for b in cfg['stages'][stage]['batches']}
            fixed=r.get(f'fixed_{stage}_batch','')
            permitted=r.get(f'allowed_{stage}_batches','')
            if fixed and fixed not in names:raise ValueError(f'{r["sample_id"]}: unknown fixed {stage} batch')
            if permitted and not set(permitted.split(';'))<=names:raise ValueError(f'{r["sample_id"]}: unknown allowed {stage} batch')
            if fixed and permitted and fixed not in permitted.split(';'):raise ValueError('Fixed batch conflicts with allowed batches')
    loss=cfg.get('loss_scenarios',{})
    if set(loss)-{'unit_loss_rates','sample_failure_rate','simulations'}:raise ValueError('Unknown loss_scenarios setting')
    rates=loss.get('unit_loss_rates',[0,.1,.2])
    if not isinstance(rates,list) or not rates:raise ValueError('unit_loss_rates needs a nonempty list')
    for value in rates+[loss.get('sample_failure_rate',0)]:
        if number(value,'loss rate')>1:raise ValueError('Loss rates must be between 0 and 1')
    simulations=integer(loss.get('simulations',500),'simulations',10)
    if simulations>10000:raise ValueError('simulations must be <=10000')
    if 'budget' in cfg:
        budget=cfg['budget'];fields={'currency','limit','per_unit','per_specimen','per_extraction','per_library'}
        if set(budget)!=fields:raise ValueError('budget must specify currency, limit, per_unit, per_specimen, per_extraction, per_library')
        if not isinstance(budget['currency'],str) or not budget['currency']:raise ValueError('Specify budget currency')
        for k in fields-{'currency'}:number(budget[k],k)
    return cfg,rows


def allowed(row,stage,names):
    fixed=row.get(f'fixed_{stage}_batch','')
    if fixed:return {fixed}
    values=row.get(f'allowed_{stage}_batches','')
    return set(values.split(';')) if values else set(names)


def balance_score(counts,totals,capacities):
    denom=sum(capacities.values())
    return sum((counts[b][key]-n*capacities[b]/denom)**2/max(1,n*capacities[b]/denom)
               for b in capacities for key,n in totals.items())


def allocate(rows,stage,spec,columns,unit,seed,restarts):
    rng=random.Random(seed);control_count=sum(spec.get('controls_per_batch',{}).values())
    capacities={b['id']:b['capacity']-control_count for b in spec['batches']};names=list(capacities)
    if len(rows)>sum(capacities.values()):raise ValueError(f'{stage}: {len(rows)} incoming samples exceed {sum(capacities.values())} available positions after controls')
    groups=defaultdict(list)
    for row in rows:
        key=('unit',row[unit]) if spec.get('keep_unit_together') and row['sample_type']=='biological' else ('sample',row['sample_id'])
        groups[key].append(row)
    bundles=[]
    for group in groups.values():
        choices=set(names)
        for row in group:choices &= allowed(row,stage,names)
        if not choices:raise ValueError(f'{stage}: contradictory restrictions within allocation unit')
        if not any(capacities[b]>=len(group) for b in choices):raise ValueError(f'{stage}: allocation unit cannot fit any permitted batch')
        bundles.append((group,choices))
    def features(row):
        if row['sample_type']!='biological':return [('control_type',row['sample_type'])]
        result=[(c,row[c]) for c in columns]
        if stage=='sequencing':result.append(('source_extraction_batch',row['extraction_batch']))
        return result
    totals=Counter(key for row in rows for key in features(row))
    best=None
    for attempt in range(restarts):
        shuffled=list(bundles);rng.shuffle(shuffled);shuffled.sort(key=lambda x:(len(x[1]),-len(x[0])))
        counts={b:Counter() for b in names};used=Counter();assignment={};failed=False
        for group,choices in shuffled:
            feasible=[b for b in names if b in choices and used[b]+len(group)<=capacities[b]]
            if not feasible:failed=True;break
            increments=Counter(key for row in group for key in features(row));rng.shuffle(feasible)
            def cost(b):
                delta=0
                for key,inc in increments.items():
                    target=totals[key]*capacities[b]/sum(capacities.values())
                    delta+=((counts[b][key]+inc-target)**2-(counts[b][key]-target)**2)/max(1,target)
                return delta + .1*((used[b]+len(group))**2-used[b]**2)/capacities[b]
            selected=min(feasible,key=cost)
            counts[selected].update(increments);used[selected]+=len(group)
            for row in group:assignment[row['sample_id']]=selected
        if failed:continue
        score=balance_score(counts,totals,capacities)
        if best is None or score<best[0]:best=(score,assignment)
    if best is None:raise ValueError(f'{stage}: randomized search found no feasible allocation; this is not proof that none exists. Increase restarts or revise constraints.')
    assignment=best[1];output=[dict(row,**{stage+'_batch':assignment[row['sample_id']]}) for row in rows]
    for b in names:
        for kind,n in sorted(spec.get('controls_per_batch',{}).items()):
            for i in range(n):
                output.append(dict(sample_id=f'CONTROL__{stage}__{b}__{kind}__{i+1}',sample_type=kind,source_stage=stage,**{stage+'_batch':b}))
        selected=[r for r in output if r[stage+'_batch']==b];rng.shuffle(selected)
        for pos,row in enumerate(selected,1):row[stage+'_position']=pos
    return output,best[0]


def audit(rows,stage,columns,unit):
    biological=[r for r in rows if r['sample_type']=='biological'];details=[];issues=[]
    batch=stage+'_batch'
    for column in columns:
        counts=Counter((r[batch],r[column]) for r in biological)
        batches=sorted({r[batch] for r in biological});levels=sorted({r[column] for r in biological});n=len(biological)
        bt=Counter(r[batch] for r in biological);lt=Counter(r[column] for r in biological)
        chi=sum((counts[b,l]-bt[b]*lt[l]/n)**2/(bt[b]*lt[l]/n) for b in batches for l in levels)
        v=math.sqrt(chi/(n*min(len(batches)-1,len(levels)-1))) if min(len(batches),len(levels))>1 else 0.
        # Disconnected incidence graph implies nonestimable additive level contrasts with batch.
        graph=defaultdict(set)
        for (b,l),count in counts.items():
            if count:graph[('batch',b)].add(('level',l));graph[('level',l)].add(('batch',b))
        pending=set(graph);components=0
        while pending:
            components+=1;stack=[pending.pop()]
            while stack:
                for neighbour in graph[stack.pop()]:
                    if neighbour in pending:pending.remove(neighbour);stack.append(neighbour)
        if len(levels)>1 and components>1:issues.append(f'{stage}: {column} and batch have disconnected overlap; some level contrasts are inseparable from batch in an additive model.')
        elif v>.5:issues.append(f'{stage}: substantial {column}/batch association (descriptive Cramer V={v:.3f}); inspect balance.')
        for b in batches:
            for l in levels:
                details.append(dict(stage=stage,factor=column,batch=b,level=l,n_samples=counts[b,l],cramers_v=v,overlap_components=components))
        if len(levels)>1 and len(batches)==1:issues.append(f'{stage}: only one occupied batch; between-batch robustness cannot be assessed.')
    return details,issues


def loss_scenarios(rows,cfg):
    setting=cfg.get('loss_scenarios',{});unit=cfg['experimental_unit'];rng=random.Random(cfg['seed']+1009)
    units=sorted({r[unit] for r in rows});groups=sorted({r['group'] for r in rows});output=[]
    bygroup={g:[r for r in rows if r['group']==g] for g in groups}
    simulations=setting.get('simulations',500);failure=setting.get('sample_failure_rate',0)
    for rate in setting.get('unit_loss_rates',[0,.1,.2]):
        counts={g:[] for g in groups};remaining={g:[] for g in groups};complete={g:[] for g in groups}
        for _ in range(simulations):
            kept={u for u in units if rng.random()>=rate}
            survived={r['sample_id'] for r in rows if r[unit] in kept and rng.random()>=failure}
            for g,group in bygroup.items():
                observed=[r for r in group if r['sample_id'] in survived]
                counts[g].append(len(observed));remaining[g].append(len({r[unit] for r in observed}))
                # A complete schedule requires every distinct specimen, not every technical replicate.
                target=defaultdict(set);seen=defaultdict(set)
                for r in group:target[r[unit]].add(r['specimen_id'])
                for r in observed:seen[r[unit]].add(r['specimen_id'])
                complete[g].append(sum(seen[u]==samples for u,samples in target.items()))
        for g in groups:
            def quantile(x,q):return sorted(x)[min(len(x)-1,int((len(x)-1)*q))]
            output.append(dict(group=g,unit_loss_rate=rate,sample_failure_rate=failure,simulations=simulations,
                mean_remaining_samples=statistics.mean(counts[g]),mean_remaining_units=statistics.mean(remaining[g]),
                units_p05=quantile(remaining[g],.05),units_p95=quantile(remaining[g],.95),
                mean_complete_schedule_units=statistics.mean(complete[g]),interpretation='Operational attrition scenario, not statistical power'))
    return output


def budget_summary(rows,sheet,cfg):
    if 'budget' not in cfg:return dict(status='not_assessed',reason='No costs/budget supplied')
    b=cfg['budget'];unit=cfg['experimental_unit']
    quantities=dict(per_unit=len({r[unit] for r in rows}),per_specimen=len({r['specimen_id'] for r in rows}),
        per_extraction=sum(bool(r.get('extraction_batch')) for r in sheet),per_library=len(sheet))
    total=sum(quantities[k]*b[k] for k in quantities)
    return dict(status='within_budget' if total<=b['limit'] else 'over_budget',currency=b['currency'],estimated_total=total,limit=b['limit'],quantities=quantities,
        assumptions='All allocated samples and controls incur the supplied stage costs. No fixed fees, failed-run repeats or storage costs unless included in those rates.')


def write_tsv(path,rows,fields):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(rows)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
