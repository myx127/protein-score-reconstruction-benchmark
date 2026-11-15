"""Independent toy data; no HMACS values, fitted covariances or IDs are used."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from .frozen import fold_seed
from .io import save_json

def generate(directory,n_bootstrap=50,module_sizes=(8,12,16),n_rows=(160,180),seed=20261115):
    if n_bootstrap<1 or not 1<=len(module_sizes)<=3 or min(module_sizes)<4 or len(n_rows)!=2 or min(n_rows)<125:raise ValueError('Use positive bootstrap count, 1-3 modules of >=4 proteins, and two batches of >=125 rows')
    out=Path(directory).resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'SYNTHETIC_DATA_MANIFEST.json').exists():raise FileExistsError('Choose a new synthetic input directory')
    rng=np.random.default_rng(seed);modules=['toy_'+name for name in ['black','green','turquoise'][:len(module_sizes)]]
    weights=[];members={};position=0
    for m,p in zip(modules,module_sizes):
        proteins=[f'SYN_P{position+j:04d}' for j in range(p)];position+=p;members[m]=proteins
        w=rng.uniform(.3,1,p)*rng.choice([-1,1],p,p=[.25,.75])
        weights.extend(dict(source_system='synthetic',module=m,protein=protein,weight=float(value)) for protein,value in zip(proteins,w))
    weights=pd.DataFrame(weights);weights.to_csv(out/'weights.csv',index=False,float_format='%.17g')
    matrices={}
    for batch,n in zip(['SYN_Batch1','SYN_Batch2'],n_rows):
        columns=[]
        for m,p in zip(modules,module_sizes):
            latent=rng.normal(size=(n,2));load=rng.uniform(.45,.9,p)
            x=latent[:,[0]]*load+latent[:,[1]]*rng.normal(0,.15,p)+rng.normal(size=(n,p))*.55
            x=(x-x.mean(axis=0))/x.std(axis=0,ddof=1);columns.append(x)
        x=np.column_stack(columns);d=pd.DataFrame(x,columns=weights.protein)
        d.insert(0,'record_id',[f'{batch}_RECORD_{i:04d}' for i in range(n)])
        d.insert(0,'participant_id',[f'{batch}_PERSON_{i:04d}' for i in range(n)])
        d.to_csv(out/(batch+'.csv'),index=False,float_format='%.17g');matrices[batch]=d
    masks=[];units=[]
    for m,proteins in members.items():
        for label,fraction in [('toy25',.25),('toy50',.5),('toy70',.7),('toy90',.9)]:
            mask_id='SYN|'+m+'|'+label
            k=max(2,min(len(proteins)-1,round(len(proteins)*fraction)))
            observed=rng.choice(proteins,k,replace=False).tolist()
            masks.extend(dict(mask_id=mask_id,source_system='synthetic',module=m,target_dataset='synthetic',platform_branch=label,protein=p) for p in observed)
            for direction,source,target in [('B1_to_B2','SYN_Batch1','SYN_Batch2'),('B2_to_B1','SYN_Batch2','SYN_Batch1')]:
                units.append(dict(direction_unit_id=mask_id+'|'+direction,mask_unit_id=mask_id,module=m,platform_label=label,direction=direction,train_batch=source,test_batch=target,full_module_proteins=len(proteins),observed_proteins=k))
    pd.DataFrame(masks).to_csv(out/'masks.csv',index=False);units=pd.DataFrame(units);units.to_csv(out/'units.csv',index=False)
    seeds=[]
    for batch in matrices:
        for rep in range(1,n_bootstrap+1):seeds.append(dict(stream=batch+'_target',replicate=rep,seed=fold_seed(seed,batch+'|'+str(rep)),role='primary_target_bootstrap'))
    for u in units.itertuples():seeds.append(dict(stream='fold::'+u.direction_unit_id,replicate=0,seed=fold_seed(seed,'fold::'+u.direction_unit_id),role='source_and_target_fold_assignment'))
    pd.DataFrame(seeds).to_csv(out/'seeds.csv',index=False)
    config=dict(data_kind='synthetic',batches={k:k+'.csv' for k in matrices},weights='weights.csv',masks='masks.csv',units='units.csv',seeds='seeds.csv',n_bootstrap=n_bootstrap,enforce_paper_fingerprints=False)
    save_json(config,out/'recovery.json')
    # Independently simulated, same-person measurements for six assay-fluid
    # branches. Their toy features and weights are not real assay annotations.
    ew=weights.rename(columns={'module':'ModuleColor','protein':'ProteinID','weight':'kME_ownModule'})
    ew[['ModuleColor','ProteinID','kME_ownModule']].to_csv(out/'emory_weights.csv',index=False,float_format='%.17g')
    source=matrices['SYN_Batch1'].drop(columns='record_id').rename(columns={'participant_id':'sample_id'})
    source.to_csv(out/'definition_matrix.csv',index=False,float_format='%.17g')
    n=35;ids=[f'SYN_EMORY_{i:03d}' for i in range(n)]
    true=rng.normal(size=(n,len(weights)));common=rng.normal(size=(n,len(modules)))
    for j,m in enumerate(modules):
        cols=[weights.index[weights.protein==p][0] for p in members[m]];true[:,cols]=.7*common[:,[j]]+.4*true[:,cols]
    mapping=[];branches=[]
    for platform,short,fraction in [('SomaScan','soma',.75),('Olink','olink',.65),('TMT-MS','tmt',.55)]:
        for fluid in ['CSF','plasma']:
            input_id=short+'_'+fluid.lower();observed=[]
            for m,proteins in members.items():observed.extend(rng.choice(proteins,max(2,min(len(proteins)-1,round(len(proteins)*fraction))),replace=False).tolist())
            data={}
            for j,p in enumerate(observed):
                f=input_id+'_F'+str(j);k=int(weights.index[weights.protein==p][0])
                value=true[:,k]+rng.normal(0,.3,n)+(0 if fluid=='CSF' else rng.normal(0,.2,n))
                value=(2**(10+.3*value) if platform=='SomaScan' else 2**(.3*value) if platform=='TMT-MS' else 10+.3*value)
                if j%5==0:value[rng.choice(n,2,replace=False)]=np.nan
                data[f]=value;mapping.append(dict(input_id=input_id,feature_id=f,protein_id=p,conservative_eligible=True,mapping_route='independent_synthetic_one_to_one'))
            d=pd.DataFrame(data);d.insert(0,'sample_id',ids);d.to_csv(out/(input_id+'.csv'),index=False,float_format='%.17g')
            branches.append(dict(input_id=input_id,platform=platform,biofluid=fluid,path=input_id+'.csv',input_scale='raw'))
    pd.DataFrame(mapping).to_csv(out/'emory_mapping.csv',index=False)
    save_json(dict(data_kind='synthetic',source_matrix='definition_matrix.csv',weights='emory_weights.csv',mapping='emory_mapping.csv',branches=branches,n_bootstrap=n_bootstrap,strict_paper_design=False),out/'emory.json')
    save_json(dict(data_kind='synthetic',seed=seed,modules=len(modules),proteins=len(weights),settings=len(units),source_rows=n_rows[0],target_rows=n_rows[1],paired_assay_participants=n,HMACS_inputs_used=False),out/'SYNTHETIC_DATA_MANIFEST.json')
    return out/'recovery.json',out/'emory.json'
