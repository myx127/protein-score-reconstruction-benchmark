"""Paired-assay evaluation using governed source data and original assay inputs."""
from __future__ import annotations
import json
import pandas as pd
from .io import ROOT,enforce_versions,load_config,outside_repository,resolve_path,run_r_job,save_json,work_directory

def run(config_path,workdir,n_boot=None,allow_version_drift=False):
    config,base=load_config(config_path);synthetic=config.get('data_kind')=='synthetic'
    if config.get('data_kind') not in {'synthetic','governed'}:raise ValueError('data_kind must be synthetic or governed')
    out=work_directory(workdir,synthetic)
    if (out/'RUN_STATUS.json').exists():raise FileExistsError('Choose a new Emory work directory')
    source=resolve_path(config['source_matrix'],base)
    branches=[]
    for item in config['branches']:
        b=dict(item,path=str(resolve_path(item['path'],base)));branches.append(b)
        if not synthetic:outside_repository(b['path'])
    if not synthetic:outside_repository(source)
    if len(branches)!=6 or len({b['input_id'] for b in branches})!=6:raise ValueError('Exactly six unique platform-fluid branches are required')
    if {(b['platform'],b['biofluid']) for b in branches}!={(p,f) for p in ['SomaScan','Olink','TMT-MS'] for f in ['CSF','plasma']}:raise ValueError('Invalid six-branch design')
    n=int(n_boot if n_boot is not None else config.get('n_bootstrap',10000))
    if n<1:raise ValueError('Positive bootstrap count required')
    env=enforce_versions(allow_version_drift or synthetic)
    payload=dict(source_matrix=str(source),weights=str(resolve_path(config.get('weights',str(ROOT/'metadata/emory_module_weights.csv')),base)),mapping=str(resolve_path(config.get('mapping',str(ROOT/'metadata/emory_feature_mapping.csv')),base)),branches=branches,output=str(out),n_bootstrap=n,strict_paper_design=config.get('strict_paper_design',not synthetic),strict_environment=not(synthetic or allow_version_drift))
    run_r_job('emory',payload,out,config.get('rscript'))
    results=pd.read_csv(out/'emory_concordance.csv')
    status=dict(status='PASS',data_kind='synthetic' if synthetic else 'governed',participants=int(results.participants.iloc[0]),evaluable_cells=int(results[results.level=='all'].cells.iloc[0]),bootstrap_draws=n,endpoint='mean cellwise CRNP-minus-direct Fisher-z difference in Spearman correlations',environment=env,manuscript_estimates_reproduced=False if synthetic else None)
    save_json(status,out/'RUN_STATUS.json');return status
