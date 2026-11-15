"""Single-thread timing with a persistent R worker and interleaved methods."""
from __future__ import annotations
import gc,json,math,os,queue,shutil,subprocess,tempfile,threading,time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_info,threadpool_limits
from . import frozen as f
from .io import ROOT,check_design,enforce_versions,load_config,outside_repository,read_matrix,resolve_path,save_json,sha256,work_directory

METHODS=['direct','cross_batch_ridge','CRNP_v0.2.0']

class RSession:
    """Keep R resident while Python and R measurement commands interleave."""
    def __init__(self,payload,out,rscript=None):
        executable=rscript or os.environ.get('RSCRIPT') or shutil.which('Rscript')
        if not executable:raise RuntimeError('Rscript is required; see docs/INSTALL.md')
        self.parent=Path(tempfile.gettempdir()).resolve()
        self.temp=Path(tempfile.mkdtemp(prefix='psbench_r_',dir=self.parent)).resolve()
        self.proc=None;self.log=None
        try:
            shutil.copyfile(ROOT/'R/worker.R',self.temp/'worker.R')
            shutil.copytree(ROOT/'vendor/crnp/R',self.temp/'crnpR')
            save_json(dict(payload,task='timing',vendor_r='crnpR'),self.temp/'job.json')
            env=os.environ.copy();env.update(LC_ALL='C',LANG='C')
            self.log=(out/'r_timing.log').open('w',encoding='utf-8')
            self.proc=subprocess.Popen([str(executable),'--vanilla','worker.R','job.json'],cwd=self.temp,env=env,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,text=True,encoding='utf-8',errors='replace',bufsize=1)
            self.messages=queue.Queue()
            def reader():
                for line in self.proc.stdout:
                    if line.startswith('RESULT:'):self.messages.put(json.loads(line[7:]))
                    else:self.log.write(line);self.log.flush()
                self.messages.put({'action':'error','reason':'R worker exited; inspect r_timing.log'})
            self.reader=threading.Thread(target=reader,daemon=True);self.reader.start()
            self.ready=self.receive()
            if self.ready.get('action')!='ready':raise RuntimeError(str(self.ready))
        except BaseException:
            self.close();raise
    def receive(self):
        try:message=self.messages.get(timeout=3600)
        except queue.Empty:raise RuntimeError('Timed out waiting for the R timing worker')
        if message.get('action')=='error':raise RuntimeError(str(message))
        return message
    def request(self,**command):
        self.proc.stdin.write(json.dumps(command)+'\n');self.proc.stdin.flush()
        return self.receive()
    def close(self):
        if self.proc is not None:
            if self.proc.poll() is None:self.proc.terminate()
            self.proc.wait(timeout=30)
            for stream in [self.proc.stdin,self.proc.stdout]:
                if stream is not None:stream.close()
        if self.log is not None:self.log.close()
        if self.temp.parent!=self.parent or not self.temp.name.startswith('psbench_r_'):raise RuntimeError('Invalid cleanup boundary')
        shutil.rmtree(self.temp)

def _fit(c):
    # Score and residual construction belong to fit time, as in the manuscript.
    full=c['source_full']@c['w_full']/np.abs(c['w_full']).sum()
    direct=c['source_obs']@c['w_obs']/np.abs(c['w_obs']).sum()
    residual=full-direct
    alpha,_=f.select_alpha_one_se(c['source_obs'],residual,c['seed'])
    return Ridge(alpha=alpha,fit_intercept=True,solver='cholesky').fit(c['source_obs'],residual)

def _score(c,m,fit):
    direct=c['target']@c['w_obs']/np.abs(c['w_obs']).sum()
    return direct if m=='direct' else direct+fit.predict(c['target'])

def _block(c,m,fit,calls):
    start=time.perf_counter_ns()
    for _ in range(calls):last=_score(c,m,fit)
    return (time.perf_counter_ns()-start)/1e9,last

def _python_measure(c,m,phase,calls=None):
    gc.collect();start=time.perf_counter_ns();fit=None if m=='direct' else _fit(c)
    fit_s=0. if m=='direct' else (time.perf_counter_ns()-start)/1e9
    difference=float(np.max(np.abs(_score(c,m,fit)-c['expected'][m])))
    if difference>1e-7:raise RuntimeError('Repeated Python fitting changed scores')
    if phase=='warmup':
        calls=16
        while True:
            elapsed,_=_block(c,m,fit,calls)
            if elapsed>=.025 or calls>=131072:break
            calls*=2
        calls=int(min(131072,max(32,2**math.ceil(math.log2(.1/(elapsed/calls))))))
        return dict(prediction_calls=calls,max_abs_score_difference=difference)
    gc.collect();elapsed,last=_block(c,m,fit,calls)
    difference=float(np.max(np.abs(last-c['expected'][m])))
    return dict(fit_s=fit_s,prediction_ms=1000*elapsed/calls,prediction_calls=calls,max_abs_score_difference=difference,selected_alpha=None if fit is None else float(fit.alpha))

def summarize(records):
    case=records.groupby(['case_index','direction_unit_id','module','direction','method'],as_index=False)[['fit_s','prediction_ms']].median()
    rows=[]
    for method,g in case.groupby('method',sort=False):
        d=dict(method=method,case_count=len(g),formal_repeat_records=int((records.method==method).sum()))
        for metric in ['fit_s','prediction_ms']:
            d[metric+'_median']=float(g[metric].median())
            d[metric+'_q25']=float(g[metric].quantile(.25));d[metric+'_q75']=float(g[metric].quantile(.75))
        d['fit_status']='STRUCTURAL_ZERO_NO_FIT' if method=='direct' else 'MEASURED_INCLUDING_TUNING'
        d['summary_definition']='median across workload-specific repeat medians; IQR across workloads, not a confidence interval'
        rows.append(d)
    return case,pd.DataFrame(rows)

def run(config_path,workdir,repetitions=5,allow_version_drift=False):
    config,base=load_config(config_path);synthetic=config.get('data_kind')=='synthetic'
    if config.get('data_kind') not in {'synthetic','governed'}:raise ValueError('data_kind must be synthetic or governed')
    out=work_directory(workdir,synthetic)
    if (out/'RUN_STATUS.json').exists():raise FileExistsError('Choose a new timing work directory')
    repetitions=int(repetitions)
    if repetitions<1:raise ValueError('Positive repetition count required')
    environment=enforce_versions(allow_version_drift or synthetic)
    paths={k:resolve_path(v,base) for k,v in config['batches'].items()}
    if not synthetic:
        for p in paths.values():outside_repository(p)
        if config.get('enforce_paper_fingerprints',True):
            fingerprints=json.loads((ROOT/'metadata/expected_governed_inputs.json').read_text('utf-8'))
            for k,p in paths.items():
                if k not in fingerprints or sha256(p)!=fingerprints[k]['SHA256']:raise ValueError('Frozen manuscript matrix fingerprint differs')
    batches={k:read_matrix(v) for k,v in paths.items()}
    wp=resolve_path(config.get('weights',str(ROOT/'metadata/hmacs_module_weights.csv')),base)
    mp=resolve_path(config.get('masks',str(ROOT/'metadata/paper_coverage_masks.csv')),base)
    sp=resolve_path(config.get('seeds',str(ROOT/'metadata/paper_random_seeds.csv')),base)
    up=resolve_path(config.get('units',str(ROOT/'metadata/paper_direction_units.csv')),base)
    weights,masks,seeds,all_units=[pd.read_csv(p) for p in [wp,mp,sp,up]]
    if not synthetic and config.get('timing_use_paper_order',True):
        units=pd.read_csv(ROOT/'metadata/paper_timing_units.csv')
        if not set(units.direction_unit_id)<=set(all_units.direction_unit_id):raise ValueError('Timing units absent from the recovery design')
    else:
        selected=config.get('timing_modules',list(weights.module.unique()[:3]))
        units=all_units[all_units.module.isin(selected)].reset_index(drop=True).copy()
        units['case_index']=np.arange(1,len(units)+1)
    check_design(weights,masks,units,batches)
    f.atomic_csv(units,out/'timing_units.csv')
    seed_map=seeds[seeds.replicate==0].set_index('stream').seed.to_dict()
    cases=[]
    for u in units.itertuples():
        wd=weights[weights.module==u.module];w=wd.set_index('protein').weight
        obs=[p for p in masks.loc[masks.mask_id==u.mask_unit_id,'protein'] if p in w.index]
        source,target=batches[u.train_batch],batches[u.test_batch]
        c=dict(source_full=np.ascontiguousarray(f.feature_matrix(source,wd.protein.tolist())),source_obs=np.ascontiguousarray(f.feature_matrix(source,obs)),target=np.ascontiguousarray(f.feature_matrix(target,obs)),w_full=wd.weight.to_numpy(float),w_obs=w.loc[obs].to_numpy(float),seed=int(seed_map['fold::'+u.direction_unit_id]))
        cases.append(c)
    if not synthetic and config.get('timing_use_paper_order',True) and repetitions==5:
        order=pd.read_csv(ROOT/'metadata/paper_timing_order.csv')
    else:
        rng=np.random.default_rng(20261115);rows=[]
        for phase,rep in [('warmup',0)]+[('formal',r) for r in range(1,repetitions+1)]:
            pairs=[(i,m) for i in range(1,len(cases)+1) for m in METHODS]
            for j in rng.permutation(len(pairs)):
                i,m=pairs[j];rows.append(dict(sequence=len(rows)+1,phase=phase,repeat=rep,case_index=i,method=m))
        order=pd.DataFrame(rows)
    f.atomic_csv(order,out/'timing_order.csv')
    payload=dict(weights=str(wp),masks=str(mp),units=str(out/'timing_units.csv'),seeds=str(sp),batches={k:str(v) for k,v in paths.items()},output=str(out),strict_environment=not(synthetic or allow_version_drift))
    r=None;records=[];warmups=[];calls={}
    try:
        with threadpool_limits(limits=1):
            pools=threadpool_info()
            if any(p['num_threads']!=1 for p in pools):raise RuntimeError('Python numeric thread pool is not single-threaded')
            for c in cases:c['expected']={m:_score(c,m,None if m=='direct' else _fit(c)) for m in METHODS[:2]}
            r=RSession(payload,out,config.get('rscript'))
            if r.ready['cases']!=len(cases):raise RuntimeError('R/Python workload count mismatch')
            for command in order.itertuples():
                i=int(command.case_index);m=command.method;warm=command.phase=='warmup'
                phase='warmup' if warm else 'measure'
                if m=='CRNP_v0.2.0':
                    d=r.request(action=phase,case_index=i,**({} if warm else {'prediction_calls':calls[(i,m)]}))
                    d.pop('action',None)
                else:d=_python_measure(cases[i-1],m,phase,None if warm else calls[(i,m)])
                if d['max_abs_score_difference']>1e-7:raise RuntimeError('Timing score validation failed')
                u=units.iloc[i-1]
                d.update(sequence=int(command.sequence),repeat=int(order.loc[command.Index,'repeat']),case_index=i,method=m,direction_unit_id=u.direction_unit_id,module=u.module,direction=u.direction,fit_required=m!='direct')
                if warm:calls[(i,m)]=int(d['prediction_calls']);warmups.append(d)
                else:records.append(d)
                # Keep each completed observation if a later command fails.
                with (out/'timing_journal.jsonl').open('a',encoding='utf-8') as stream:stream.write(json.dumps(d,allow_nan=False)+'\n')
            r.request(action='stop')
    finally:
        if r is not None:r.close()
    measured=pd.DataFrame(records);case,overall=summarize(measured)
    for name,frame in [('formal_records.csv',measured),('warmup_records.csv',pd.DataFrame(warmups)),('case_medians.csv',case),('timing_overall.csv',overall)]:f.atomic_csv(frame,out/name)
    status=dict(status='PASS',data_kind='synthetic' if synthetic else 'governed',cases=len(cases),warmups=len(warmups),formal_records=len(measured),repetitions=repetitions,numerical_threads=1,max_abs_score_difference=float(measured.max_abs_score_difference.max()),randomized_across_languages=True,fit_includes_tuning=True,reference='independently fitted untimed predictions on identical inputs',runtime_scope='resident workers; preloaded matrices; fit includes score/residual construction; whole target batch scoring; excludes startup/IO/common preprocessing',environment=environment,r_environment=r.ready['environment'],manuscript_runtime_reproduced=False)
    save_json(status,out/'RUN_STATUS.json');return status
