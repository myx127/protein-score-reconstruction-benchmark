"""Scientific contracts using synthetic values only; run with unittest."""
from __future__ import annotations
import ast,hashlib,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pandas as pd
from psbench import frozen as f
from psbench.io import ROOT,check_design,outside_repository,read_matrix
from psbench.recovery import METHODS,aggregate,conditional_bootstrap

def fixture():
    rng=np.random.default_rng(1701)
    units=pd.read_csv(ROOT/'metadata/paper_direction_units.csv')
    seeds=pd.read_csv(ROOT/'metadata/paper_random_seeds.csv')
    batches={}
    for name,n in [('HMACS_Batch1',60),('HMACS_Batch2',70)]:
        ids=np.array([name+'_TOY_'+str(i) for i in range(n)])
        batches[name]=f.Bundle(ids,ids.copy(),[],np.empty((n,0)),{})
    rows=[]
    for u in units.itertuples():
        b=batches[u.test_batch];n=len(b.participant)
        truth=rng.normal(size=n);err=rng.normal(scale=rng.uniform(.2,.8),size=n)
        rows.append(pd.DataFrame(dict(direction_unit_id=u.direction_unit_id,test_batch=u.test_batch,participant_id=b.participant,record_id=b.record,reference_full=truth,direct=truth+err,within_cf_ridge=truth+.4*err,cross_batch_ridge=truth+.75*err,**{'CRNP_v0.2.0':truth+.85*err})))
    return units,seeds,batches,pd.concat(rows,ignore_index=True)

class FrozenCodeTests(unittest.TestCase):
    def test_extracted_bodies_and_vendored_R_hashes(self):
        source=(ROOT/'src/psbench/frozen.py').read_text('utf-8');nodes={n.name:n for n in ast.parse(source).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        for item in json.loads((ROOT/'metadata/code_provenance.json').read_text('utf-8')):
            if item['source_family']=='frozen Stage24':
                n=nodes[item['component']];body=ast.get_source_segment(source,n)
                if isinstance(n,ast.ClassDef):body='@dataclass\n'+body
                self.assertEqual(hashlib.sha256(body.encode()).hexdigest(),item['source_function_SHA256'],item['component'])
            else:
                p=ROOT/'vendor/crnp/R'/item['source_file']
                self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),item['SHA256'],str(p.name))

    def test_score_normalization_and_column_identity(self):
        x=np.array([[1.,2.,4.],[3.,-1.,2.]])
        b=f.Bundle(np.array(['a','b']),np.array(['1','2']),['p','q','r'],x,{'p':0,'q':1,'r':2})
        w=pd.Series([2.,-1.,3.],index=['p','q','r'])
        np.testing.assert_allclose(f.weighted_score(b,w,['p','q','r']),x@w.to_numpy()/6)
        observed=f.weighted_score(b,w,['q','p']);np.testing.assert_allclose(observed,(2*x[:,0]-x[:,1])/3)
        # The residual includes the changed denominator, not only omitted mass.
        self.assertFalse(np.allclose(f.weighted_score(b,w,['p','q','r'])-observed,3*x[:,2]/6))
        np.testing.assert_array_equal(f.feature_matrix(b,['r','p']),x[:,[2,0]])

    def test_nmse_uses_sample_variance_and_gain_identity(self):
        rng=np.random.default_rng(5);truth=rng.normal(size=100);direct=truth+rng.normal(size=100);pred=.8*truth+.2*direct
        m=f.metrics(truth,direct,pred)
        self.assertAlmostEqual(m['risk'],np.mean((truth-pred)**2)/np.var(truth,ddof=1),places=14)
        gain=f.metrics(truth,direct,direct)['risk']-m['risk']
        self.assertAlmostEqual(gain,m['alignment_benefit']-m['correction_cost'],places=14)

    def test_cluster_bootstrap_retains_every_record_per_draw(self):
        ids=np.repeat(['a','b','c'],2);b=f.Bundle(ids,np.arange(6).astype(str),[],np.empty((6,0)),{})
        idx=f.participant_bootstrap_indices(b,981)
        self.assertEqual(len(idx),6)
        for pair in idx.reshape(-1,2):self.assertEqual(pair[0]//2,pair[1]//2);self.assertEqual(pair[1],pair[0]+1)

    def test_nested_folds_are_disjoint(self):
        rng=np.random.default_rng(122);x=rng.normal(size=(130,4));y=x[:,0]+rng.normal(size=130)
        p,a,audit=f.within_crossfit(x,y,np.ones(130,bool),919)
        self.assertTrue(np.isfinite(p).all());self.assertEqual(len(a),5)
        self.assertTrue(all(v['train_test_overlap']==0 for v in audit))
        self.assertEqual(sum(v['test_n'] for v in audit),130)

    def test_one_se_selects_largest_eligible_alpha(self):
        rng=np.random.default_rng(11);x=rng.normal(size=(100,3));y=rng.normal(size=100)
        alpha,audit=f.select_alpha_one_se(x,y,11);best=audit.mean_cv_mse.idxmin()
        eligible=audit[audit.mean_cv_mse<=audit.loc[best,'mean_cv_mse']+audit.loc[best,'se_cv_mse']+1e-15]
        self.assertEqual(alpha,eligible.alpha.max())

class InputBoundaryTests(unittest.TestCase):
    def test_source_target_participant_overlap_is_rejected(self):
        names=['a','b','c'];ids=np.array(['SYN_1','SYN_2'])
        b=f.Bundle(ids,ids.copy(),names,np.ones((2,3)),dict(zip(names,range(3))))
        weights=pd.DataFrame({'module':['toy']*3,'protein':names,'weight':[1.,1.,1.]})
        masks=pd.DataFrame({'mask_id':['m','m'],'protein':['a','b']})
        units=pd.DataFrame([dict(direction_unit_id='u',module='toy',mask_unit_id='m',train_batch='source',test_batch='target')])
        with self.assertRaisesRegex(ValueError,'overlap'):check_design(weights,masks,units,{'source':b,'target':b})

    def test_complete_coverage_is_rejected_for_partial_mask_benchmark(self):
        names=['a','b','c'];ids=np.array(['SYN_1','SYN_2']);other=np.array(['SYN_3','SYN_4'])
        b=f.Bundle(ids,ids.copy(),names,np.ones((2,3)),dict(zip(names,range(3))))
        c=f.Bundle(other,other.copy(),names,np.ones((2,3)),dict(zip(names,range(3))))
        w=pd.DataFrame({'module':['toy']*3,'protein':names,'weight':[1.,1.,1.]});m=pd.DataFrame({'mask_id':['m']*3,'protein':names})
        u=pd.DataFrame([dict(direction_unit_id='u',module='toy',mask_unit_id='m',train_batch='source',test_batch='target')])
        with self.assertRaisesRegex(ValueError,'observed members'):check_design(w,m,u,{'source':b,'target':c})

    def test_governed_paths_cannot_enter_repository(self):
        for p in [ROOT,ROOT/'inputs/raw.csv',ROOT/'outputs/private']:
            with self.assertRaises(ValueError):outside_repository(p)
        self.assertEqual(outside_repository(ROOT.parent/'approved_work'),ROOT.parent/'approved_work')

    def test_matrix_preserves_string_identifiers(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'matrix.csv';n=70
            pd.DataFrame({'participant_id':[f'{i:03d}' for i in range(n)],'record_id':[f'{i:04d}' for i in range(n)],'p':np.arange(n)}).to_csv(p,index=False)
            b=read_matrix(p);self.assertEqual(b.participant[0],'000');self.assertEqual(b.record[0],'0000')

    def test_nonfinite_and_duplicate_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'bad.csv';d=pd.DataFrame({'participant_id':[str(i) for i in range(70)],'record_id':[str(i) for i in range(70)],'p':np.arange(70,dtype=float)})
            d.loc[2,'p']=np.nan;d.to_csv(p,index=False)
            with self.assertRaises(ValueError):read_matrix(p)
            d.loc[2,'p']=2;d.loc[2,['participant_id','record_id']]=d.loc[1,['participant_id','record_id']];d.to_csv(p,index=False)
            with self.assertRaises(ValueError):read_matrix(p)

class AggregateTests(unittest.TestCase):
    def test_ratios_follow_equal_cell_aggregation(self):
        units,seeds,batches,p=fixture();cells=f.cell_metrics(p,units);a=aggregate(cells,units)
        for r in a[a.summary_level=='overall'].itertuples():
            raw=cells[cells.method==r.method].groupby('direction').risk.mean().mean()
            self.assertAlmostEqual(r.risk,raw,places=14)
            self.assertAlmostEqual(r.relative_gain,(r.direct_risk-r.risk)/r.direct_risk,places=14)

    def test_failed_setting_is_retained_in_aggregate(self):
        units,seeds,batches,p=fixture();cells=f.cell_metrics(p,units);cells.loc[cells.index[0],'status']='FAILED'
        a=aggregate(cells,units);r=a[(a.summary_id=='equal_direction')&(a.method=='direct')].iloc[0]
        self.assertEqual(r.status,'FAILED_RETAINED');self.assertTrue(np.isnan(r.risk))

    def test_fast_bootstrap_matches_literal_paired_resampling(self):
        units,seeds,batches,p=fixture();fast=conditional_bootstrap(p,units,batches,seeds,2)
        lookup={(r.stream,int(r.replicate)):int(r.seed) for r in seeds.itertuples()};slow=[]
        for rep in [1,2]:
            sampled={k:f.participant_bootstrap_indices(b,lookup[(k+'_target',rep)]) for k,b in batches.items()}
            frames=[g.reset_index(drop=True).iloc[sampled[g.test_batch.iloc[0]]] for _,g in p.groupby('direction_unit_id',sort=False)]
            slow.append(aggregate(f.cell_metrics(pd.concat(frames,ignore_index=True),units),units,rep))
        slow=pd.concat(slow,ignore_index=True);cols=['risk','direct_risk','method_gain','relative_gain','empirical_headroom','transport_loss','captured_headroom']
        np.testing.assert_allclose(fast[cols],slow[cols],rtol=1e-12,atol=1e-14,equal_nan=True)

    def test_published_aggregates_reconcile_with_72_cells(self):
        d=pd.read_csv(ROOT/'results/published/HMACS_72_settings.csv');o=pd.read_csv(ROOT/'results/published/HMACS_overall.csv')
        self.assertEqual(len(d),72);self.assertEqual(d.groupby('direction').size().tolist(),[36,36])
        for r in o.itertuples():
            value=d.groupby('direction')[r.method+'__risk'].mean().mean()
            self.assertAlmostEqual(r.risk,value,places=12)
        for method in ['cross_batch_ridge','CRNP_v0.2.0']:
            self.assertEqual(int((d[method+'__relative_gain']>0).sum()),60)

if __name__=='__main__':unittest.main()
