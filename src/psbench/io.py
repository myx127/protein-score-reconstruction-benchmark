"""Input contracts and output boundaries; governed data stay outside the repo."""
from __future__ import annotations
import csv,gzip,hashlib,json,os,shutil,subprocess,tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from .frozen import Bundle

ROOT=Path(__file__).resolve().parents[2]

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def save_json(value,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=True,allow_nan=False)+'\n',encoding='utf-8')

def load_config(path):
    path=Path(path).expanduser().resolve()
    return json.loads(path.read_text('utf-8-sig')),path.parent

def resolve_path(value,base):
    p=Path(value).expanduser()
    return (p if p.is_absolute() else Path(base)/p).resolve()

def outside_repository(path):
    p=Path(path).expanduser().resolve()
    if p==ROOT or ROOT in p.parents:
        raise ValueError('Governed inputs and work directories must be outside the public repository')
    return p

def work_directory(path,synthetic=False):
    p=Path(path).expanduser().resolve()
    if not synthetic:outside_repository(p)
    p.mkdir(parents=True,exist_ok=True)
    return p

def read_matrix(path):
    path=Path(path)
    opener=gzip.open if path.suffix.lower()=='.gz' else open
    with opener(path,'rt',encoding='utf-8-sig',newline='') as f:header=next(csv.reader(f))
    if header[:2]!=['participant_id','record_id'] or len(header)!=len(set(header)):
        raise ValueError('Matrix needs unique columns: participant_id,record_id,then protein columns')
    d=pd.read_csv(path,dtype={'participant_id':str,'record_id':str})
    if d[['participant_id','record_id']].isna().any().any() or d.duplicated(['participant_id','record_id']).any():
        raise ValueError('Missing or duplicate participant/record keys')
    x=d.iloc[:,2:].apply(pd.to_numeric,errors='raise').to_numpy(float)
    if len(d)<63 or not np.isfinite(x).all():raise ValueError('At least 63 rows (outer-fold training >=50) and finite, preprocessed protein values are required')
    proteins=header[2:]
    return Bundle(d.participant_id.to_numpy(),d.record_id.to_numpy(),proteins,x,{p:i for i,p in enumerate(proteins)})

def check_design(weights,masks,units,batches):
    if not {'module','protein','weight'}<=set(weights):raise ValueError('Invalid weight schema')
    if weights.duplicated(['module','protein']).any() or weights.protein.isna().any():raise ValueError('Non-unique module members')
    if not np.isfinite(weights.weight.to_numpy(float)).all():raise ValueError('Nonfinite weights')
    if units.direction_unit_id.duplicated().any():raise ValueError('Duplicate direction units')
    for u in units.itertuples():
        w=weights[weights.module==u.module];obs=masks[masks.mask_id==u.mask_unit_id].protein.tolist()
        if len(obs)!=len(set(obs)) or not set(obs)<=set(w.protein) or not 2<=len(obs)<len(w):raise ValueError('Masks require 2 <= observed members < module size')
        if abs(w.weight).sum()<=0:raise ValueError('Zero module weight norm')
        source,target=batches[u.train_batch],batches[u.test_batch]
        if set(source.participant)&set(target.participant):raise ValueError('Source/test participants overlap')
        if not set(w.protein)<=set(source.proteins) or not set(w.protein)<=set(target.proteins):raise ValueError('Recovery evaluation needs complete reference protein columns')

def run_r_job(task,payload,workdir,rscript=None):
    executable=rscript or os.environ.get('RSCRIPT') or shutil.which('Rscript')
    if not executable:raise RuntimeError('Rscript is required for CRNP. Install R and see docs/INSTALL.md; no Python substitute is used.')
    # Only code and a job description enter this temporary directory. Data paths
    # refer to the caller's synthetic or governed work area, never the release.
    parent=Path(tempfile.gettempdir()).resolve()
    temp=Path(tempfile.mkdtemp(prefix='psbench_r_',dir=parent)).resolve()
    if temp.parent!=parent or not temp.name.startswith('psbench_r_'):raise RuntimeError('Invalid temporary code directory')
    try:
        shutil.copyfile(ROOT/'R/worker.R',temp/'worker.R')
        shutil.copytree(ROOT/'vendor/crnp/R',temp/'crnpR')
        payload=dict(payload,task=task,vendor_r='crnpR')
        save_json(payload,temp/'job.json')
        env=os.environ.copy()
        for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:env[key]='1'
        env['LC_ALL']='C';env['LANG']='C'
        p=subprocess.run([str(executable),'--vanilla','worker.R','job.json'],cwd=temp,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=86400)
        log=p.stdout.decode('utf-8','replace')
        (Path(workdir)/('r_'+task+'.log')).write_text(log,encoding='utf-8')
        if p.returncode:raise RuntimeError('R '+task+' failed; inspect the log in the selected work directory')
    finally:
        if temp.parent!=parent or not temp.name.startswith('psbench_r_'):raise RuntimeError('Invalid cleanup boundary')
        shutil.rmtree(temp)

def versions():
    import platform,scipy,sklearn
    return {'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'scikit-learn':sklearn.__version__}

def enforce_versions(allow_drift=False):
    expected=json.loads((ROOT/'environment/manuscript_versions.json').read_text('utf-8'))
    observed=versions();mismatch={k:{'expected':expected[k],'observed':v} for k,v in observed.items() if expected[k]!=v}
    if mismatch and not allow_drift:raise ValueError('Frozen analysis environment differs: '+json.dumps(mismatch)+'; use the lock file or explicitly allow version drift')
    return {'observed':observed,'mismatch':mismatch}
