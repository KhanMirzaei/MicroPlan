import copy
import csv
import json
import tempfile
import unittest
from pathlib import Path
from microplan.core import allocate,audit,validate,loss_scenarios,budget_summary
from microplan.planner import run
from microplan.cli import main


def samples():
    return [dict(sample_id=f's{i}_{t}',specimen_id=f'p{i}_{t}',subject_id=f'{i:03}',group='a' if i<4 else 'b',timepoint=str(t),sample_type='biological',source_stage='input') for i in range(8) for t in range(2)]


def config():
    return dict(samples='samples.csv',experimental_unit='subject_id',assay='16s',endpoint='community',balance_columns=['group','timepoint'],seed=12,search_restarts=20,
        stages=dict(extraction=dict(batches=[dict(id='E1',capacity=12),dict(id='E2',capacity=12)],controls_per_batch={'extraction_blank':1,'mock_community':1},keep_unit_together=True),
                    sequencing=dict(batches=[dict(id='S1',capacity=16),dict(id='S2',capacity=16)],controls_per_batch={'library_blank':1},keep_unit_together=True)),
        loss_scenarios=dict(unit_loss_rates=[0,.2],sample_failure_rate=0,simulations=30))


class PlannerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.cfg=config();self.rows=samples();self.save()
    def tearDown(self):self.temp.cleanup()
    def save(self):
        with (self.root/'samples.csv').open('w',newline='') as f:
            fields=[k for k in self.rows[0] if k not in ('sample_type','source_stage')]
            w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(self.rows)
        (self.root/'study.json').write_text(json.dumps(self.cfg))
    def plan_stage(self,rows=None):
        return allocate(rows or self.rows,'extraction',self.cfg['stages']['extraction'],self.cfg['balance_columns'],'subject_id',12,20)
    def test_reproducible_allocation(self):
        a=self.plan_stage();b=self.plan_stage();self.assertEqual(a,b)
    def test_identifiers_preserved(self):
        cfg,rows=validate(self.cfg,self.root)
        self.assertEqual(rows[0]['subject_id'],'000')
    def test_capacity_and_unique_positions(self):
        rows,_=self.plan_stage()
        for b in ['E1','E2']:
            selected=[r for r in rows if r['extraction_batch']==b]
            self.assertLessEqual(len(selected),12)
            self.assertEqual(sorted(r['extraction_position'] for r in selected),list(range(1,len(selected)+1)))
    def test_unit_kept_together(self):
        rows,_=self.plan_stage()
        for subject in {r['subject_id'] for r in self.rows}:
            self.assertEqual(len({r['extraction_batch'] for r in rows if r.get('subject_id')==subject}),1)
    def test_fixed_batch_respected(self):
        for r in self.rows:
            if r['subject_id']=='000':r['fixed_extraction_batch']='E2'
        rows,_=self.plan_stage()
        self.assertTrue(all(r['extraction_batch']=='E2' for r in rows if r.get('subject_id')=='000'))
    def test_allowed_batches_respected(self):
        for r in self.rows:
            if r['subject_id']=='000':r['allowed_extraction_batches']='E1'
        rows,_=self.plan_stage()
        self.assertTrue(all(r['extraction_batch']=='E1' for r in rows if r.get('subject_id')=='000'))
    def test_conflicting_unit_restrictions(self):
        self.rows[0]['fixed_extraction_batch']='E1';self.rows[1]['fixed_extraction_batch']='E2'
        with self.assertRaisesRegex(ValueError,'contradictory'):self.plan_stage()
    def test_overcapacity(self):
        for b in self.cfg['stages']['extraction']['batches']:b['capacity']=8
        with self.assertRaisesRegex(ValueError,'exceed'):self.plan_stage()
    def test_carry_controls_forward(self):
        first,_=self.plan_stage()
        second,_=allocate(first,'sequencing',self.cfg['stages']['sequencing'],['group','timepoint'],'subject_id',13,20)
        ids={r['sample_id'] for r in first if r['sample_type']!='biological'}
        self.assertEqual(len(ids),4);self.assertTrue(ids<={r['sample_id'] for r in second})
        self.assertEqual(len(second),22)
        self.assertTrue(all(not r.get('extraction_batch') for r in second if r['source_stage']=='sequencing'))
    def test_complete_confounding_detected(self):
        rows=[dict(r,extraction_batch='E1' if r['group']=='a' else 'E2') for r in self.rows]
        details,issues=audit(rows,'extraction',['group'],'subject_id')
        self.assertTrue(any('inseparable' in x for x in issues));self.assertAlmostEqual(details[0]['cramers_v'],1)
    def test_balanced_groups_not_flagged_confounding(self):
        rows,_=self.plan_stage();details,issues=audit(rows,'extraction',['group'],'subject_id')
        self.assertFalse(any('inseparable' in x for x in issues))
        self.assertEqual(details[0]['overlap_components'],1)
    def test_no_loss_exact_counts(self):
        result=loss_scenarios(self.rows,self.cfg)
        zero=[r for r in result if r['unit_loss_rate']==0]
        self.assertTrue(all(r['mean_remaining_units']==4 and r['mean_remaining_samples']==8 for r in zero))
    def test_total_unit_loss(self):
        self.cfg['loss_scenarios']['unit_loss_rates']=[1]
        self.assertTrue(all(r['mean_remaining_samples']==0 for r in loss_scenarios(self.rows,self.cfg)))
    def test_all_sample_failure(self):
        self.cfg['loss_scenarios']['sample_failure_rate']=1
        self.assertTrue(all(r['mean_remaining_units']==0 for r in loss_scenarios(self.rows,self.cfg)))
    def test_technical_replicates_not_extra_units_or_specimens(self):
        self.rows.append(dict(self.rows[0],sample_id='technical_rep'))
        self.save();summary=run(self.root/'study.json',self.root/'out')
        self.assertEqual(summary['independent_units'],8)
        report=(self.root/'out/REPORT.md').read_text();self.assertIn('technical replicates',report.lower())
    def test_budget_includes_controls(self):
        first,_=self.plan_stage();sheet,_=allocate(first,'sequencing',self.cfg['stages']['sequencing'],['group'],'subject_id',13,20)
        self.cfg['budget']=dict(currency='EUR',limit=1000,per_unit=1,per_specimen=2,per_extraction=3,per_library=4)
        budget=budget_summary(self.rows,sheet,self.cfg)
        self.assertEqual(budget['estimated_total'],8+32+60+88)
    def test_reject_duplicate_id(self):
        self.rows[1]['sample_id']=self.rows[0]['sample_id'];self.save()
        with self.assertRaisesRegex(ValueError,'Duplicate sample_id'):validate(self.cfg,self.root)
    def test_reject_unit_in_multiple_groups(self):
        self.rows[1]['group']='b';self.save()
        with self.assertRaisesRegex(ValueError,'multiple groups'):validate(self.cfg,self.root)
    def test_reject_inconsistent_specimen(self):
        self.rows[1]['specimen_id']=self.rows[0]['specimen_id'];self.save()
        with self.assertRaisesRegex(ValueError,'inconsistent'):validate(self.cfg,self.root)
    def test_reject_unknown_setting(self):
        self.cfg['power']=.8
        with self.assertRaisesRegex(ValueError,'Unknown'):validate(self.cfg,self.root)
    def test_reject_unknown_batch(self):
        for r in self.rows:r['fixed_extraction_batch']='BAD'
        self.save()
        with self.assertRaisesRegex(ValueError,'unknown fixed'):validate(self.cfg,self.root)
    def test_reject_invalid_loss(self):
        self.cfg['loss_scenarios']['unit_loss_rates']=[1.1]
        with self.assertRaisesRegex(ValueError,'between'):validate(self.cfg,self.root)
    def test_full_report_and_provenance(self):
        summary=run(self.root/'study.json',self.root/'out')
        self.assertEqual(summary['biological_samples'],16)
        for file in ['README.md','REPORT.md','sample_sheet.tsv','ANALYSIS_PLAN.md','loss_scenarios.tsv','capacity.tsv','manifest.json','COMPLETE']:
            self.assertTrue((self.root/'out'/file).exists())
        manifest=json.loads((self.root/'out/manifest.json').read_text())
        self.assertEqual(manifest['power_analysis'],'not_performed')
        self.assertEqual(len(manifest['inputs']),2)
        self.assertFalse(list((self.root/'out').glob('*.html'))+list((self.root/'out').glob('*.png')))
    def test_no_overwrite(self):
        out=self.root/'out';out.mkdir();(out/'keep.txt').write_text('keep')
        with self.assertRaisesRegex(ValueError,'empty'):run(self.root/'study.json',out)
        self.assertEqual((out/'keep.txt').read_text(),'keep')
    def test_infeasible_does_not_leave_output(self):
        for b in self.cfg['stages']['sequencing']['batches']:b['capacity']=3
        self.save()
        with self.assertRaises(ValueError):run(self.root/'study.json',self.root/'out')
        self.assertFalse((self.root/'out').exists())
    def test_cli_invalid_returns_error(self):
        with self.assertRaises(SystemExit) as e:main(['plan','--config',str(self.root/'absent.json'),'--out',str(self.root/'out')])
        self.assertEqual(e.exception.code,2)
