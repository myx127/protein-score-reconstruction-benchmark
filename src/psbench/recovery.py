"""Portable recovery adapter and paired fixed-prediction inference."""
from __future__ import annotations
import json,math
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from . import frozen as f
from .io import ROOT,check_design,enforce_versions,load_config,outside_repository,read_matrix,resolve_path,run_r_job,save_json,sha256,work_directory

METHODS=['direct','within_cf_ridge','cross_batch_ridge','CRNP_v0.2.0']

def aggregate(cells,units,replicate=0):
    """Equal module-mask means within direction, then equal direction means.

    Unlike the original fixed-size driver this accepts small synthetic designs.
    The manuscript design still requires exactly 36 settings in each direction.
    """
    directions=list(dict.fromkeys(units.direction))
    expected=units.groupby('direction').size().to_dict()
    by_direction={}
    for d in directions:
        block=cells[cells.direction==d];values={}
        for m in METHODS:
            part=block[(block.method==m)&(block.status=='PASS')]
            values['risk_'+m]=float(part.risk.mean()) if len(part)==expected[d] else math.nan
        base=block[(block.method=='direct')&(block.status=='PASS')]
        for k in ['empirical_headroom','transport_loss']:
            values[k]=float(base[k].mean()) if len(base)==expected[d] else math.nan
        by_direction[d]=values
    overall={k:float(np.mean([by_direction[d][k] for d in directions])) for k in next(iter(by_direction.values()))}
    rows=[]
    for level,label,values in [('direction',d,by_direction[d]) for d in directions]+[('overall','equal_direction',overall)]:
        rd,h=values['risk_direct'],values['empirical_headroom']
        for m in METHODS:
            risk=values['risk_'+m];gain=rd-risk
            rows.append(dict(replicate=replicate,summary_level=level,summary_id=label,method=m,risk=risk,direct_risk=rd,method_gain=gain,relative_gain=gain/rd if rd>0 else math.nan,empirical_headroom=h,transport_loss=values['transport_loss'],captured_headroom=gain/h if h>0 and m!='direct' else math.nan,status='PASS' if np.isfinite(risk) else 'FAILED_RETAINED'))
    return pd.DataFrame(rows)

def conditional_bootstrap(predictions,units,batches,seeds,n_boot):
    """Resample target participants once per batch/draw and reuse across cells.

    Only risk-related quantities are recalculated. The original score vectors,
    fits and preprocessing remain fixed; no bootstrap refitting is implied.
    """
    lookup={(str(r.stream),int(r.replicate)):int(r.seed) for r in seeds.itertuples()}
    cached=[]
    for u in units.itertuples():
        q=predictions[predictions.direction_unit_id==u.direction_unit_id].set_index(['participant_id','record_id'])
        b=batches[u.test_batch];q=q.loc[pd.MultiIndex.from_arrays([b.participant,b.record])]
        truth=q.reference_full.to_numpy(float);pred=q[METHODS].to_numpy(float)
        if not np.isfinite(truth).all() or not np.isfinite(pred).all():raise ValueError('Conditional bootstrap requires finite locked prediction vectors')
        cached.append((u,truth,(truth[:,None]-pred)**2))
    out=[]
    for rep in range(1,n_boot+1):
        sampled={batch:f.participant_bootstrap_indices(b,lookup[(batch+'_target',rep)]) for batch,b in batches.items()}
        rows=[]
        for u,truth,squares in cached:
            idx=sampled[u.test_batch];variance=np.var(truth[idx],ddof=1)
            risks=np.mean(squares[idx,:],axis=0)/variance if variance>0 else np.full(4,np.nan)
            headroom=risks[0]-risks[1];transport=risks[2]-risks[1]
            for j,m in enumerate(METHODS):
                rows.append(dict(direction=u.direction,method=m,risk=float(risks[j]),empirical_headroom=float(headroom),transport_loss=float(transport),status='PASS' if np.isfinite(risks[j]) else 'FAILED_RETAINED'))
        out.append(aggregate(pd.DataFrame(rows),units,rep))
    return pd.concat(out,ignore_index=True) if out else pd.DataFrame()

def run(config_path,workdir,n_boot=None,allow_version_drift=False):
    config,base=load_config(config_path);synthetic=config.get('data_kind')=='synthetic'
    if config.get('data_kind') not in {'synthetic','governed'}:raise ValueError('data_kind must be synthetic or governed')
    out=work_directory(workdir,synthetic)
    if (out/'RUN_STATUS.json').exists():raise FileExistsError('Run directory already contains a completed run; choose a new directory')
    environment=enforce_versions(allow_version_drift or synthetic)
    paths={k:resolve_path(v,base) for k,v in config['batches'].items()}
    if not synthetic:
        for p in paths.values():outside_repository(p)
    batches={k:read_matrix(p) for k,p in paths.items()}
    weights_path=resolve_path(config.get('weights',str(ROOT/'metadata/hmacs_module_weights.csv')),base)
    masks_path=resolve_path(config.get('masks',str(ROOT/'metadata/paper_coverage_masks.csv')),base)
    units_path=resolve_path(config.get('units',str(ROOT/'metadata/paper_direction_units.csv')),base)
    seeds_path=resolve_path(config.get('seeds',str(ROOT/'metadata/paper_random_seeds.csv')),base)
    weights,masks,units,seeds=[pd.read_csv(p) for p in [weights_path,masks_path,units_path,seeds_path]]
    check_design(weights,masks,units,batches)
    if not synthetic and config.get('enforce_paper_fingerprints',True):
        expected=json.loads((ROOT/'metadata/expected_governed_inputs.json').read_text('utf-8'))
        for name,p in paths.items():
            if name not in expected or sha256(p)!=expected[name]['SHA256']:raise ValueError('Input differs from the frozen manuscript matrix; supply the approved frozen version or explicitly disable fingerprint enforcement for a new-data analysis')
    count=int(n_boot if n_boot is not None else config.get('n_bootstrap',2000))
    available=seeds[seeds.role=='primary_target_bootstrap'].replicate.max()
    if count<1 or count>int(available):raise ValueError('Requested bootstrap count exceeds the registered seed table')
    r_out=out/'r_scores';r_out.mkdir()
    payload=dict(weights=str(weights_path),masks=str(masks_path),units=str(units_path),seeds=str(seeds_path),batches={k:str(p) for k,p in paths.items()},output=str(r_out),strict_environment=not(synthetic or allow_version_drift))
    run_r_job('recovery',payload,out,config.get('rscript'))
    crnp=pd.read_csv(r_out/'crnp_scores.csv',dtype={'participant_id':str,'record_id':str})
    with threadpool_limits(limits=1):
        predictions,failures,audit,folds=f.compute_predictions('synthetic' if synthetic else 'governed',out,batches,weights,masks,units,seeds,crnp)
        f.atomic_csv(failures,out/'failures.csv')
        if len(failures) or len(predictions.direction_unit_id.unique())!=len(units):raise RuntimeError('Some settings failed; failure records retained in selected work directory')
        f.atomic_csv(predictions,out/'restricted_predictions.csv.gz',compress=True)
        f.atomic_csv(audit,out/'fit_audit.csv');f.atomic_csv(folds,out/'nested_fold_audit.csv')
        cells=f.cell_metrics(predictions,units);point=aggregate(cells,units)
        boot=conditional_bootstrap(predictions,units,batches,seeds,count)
        inference=f.percentile_inference(point,boot,out)
    for name,d in [('cell_metrics.csv',cells),('point_estimates.csv',point),('bootstrap_aggregates.csv',boot),('confidence_intervals.csv',inference)]:f.atomic_csv(d,out/'aggregate'/name)
    status=dict(status='PASS',data_kind='synthetic' if synthetic else 'governed',settings=len(units),directions=units.direction.nunique(),bootstrap_draws=count,methods=METHODS,crnp_version='0.2.0',reference='normalized fixed-weight score',interval_scope='fixed preprocessing/models/predictions; paired target participant bootstrap',environment=environment,manuscript_estimates_reproduced=False if synthetic else None)
    save_json(status,out/'RUN_STATUS.json')
    return status
