"""Check public file types, text, CSV schemas and optional SHA256 manifest."""
from __future__ import annotations
import argparse,csv,hashlib,json,re,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
EXCLUDED={'.git','.venv','venv','__pycache__','outputs','build','dist','.pytest_cache'}
SERIALIZED={'.rds','.rda','.rdata','.pkl','.pickle','.joblib','.h5','.hdf5','.parquet','.gz'}
TEXT={'.py','.r','.rd','.md','.json','.csv','.toml','.txt','.cff','.yaml','.yml'}
MANIFEST='SHA256SUMS.json'
REPORT='validation/public_release_audit.json'

def files():
    result=[]
    for p in ROOT.rglob('*'):
        rel=p.relative_to(ROOT)
        if any(part in EXCLUDED or part.endswith('.egg-info') for part in rel.parts):continue
        if p.is_symlink():raise ValueError('Symlinks are not permitted in the public release: '+rel.as_posix())
        if p.is_file() and p.suffix.lower() not in {'.zip','.pyc','.pyo'}:result.append(p)
    return sorted(result,key=lambda p:p.relative_to(ROOT).as_posix())

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def audit(paths):
    errors=[];csv_schemas=[];scanned=0
    for p in paths:
        rel=p.relative_to(ROOT).as_posix()
        if rel in {MANIFEST,REPORT}:continue
        scanned+=1
        if p.suffix.lower() in SERIALIZED:errors.append(rel+': quantitative matrix/archive/serialized-object extension is not allowed')
        if p.suffix.lower() in TEXT or p.name in {'LICENSE','NAMESPACE','DESCRIPTION','.gitignore'}:
            try:text=p.read_text('utf-8-sig')
            except UnicodeDecodeError:errors.append(rel+': expected UTF-8 text');continue
            if re.search(r'(?<![A-Za-z0-9_])[A-Za-z]:[\\/]',text):errors.append(rel+': absolute drive path')
            if re.search(r'gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}',text):errors.append(rel+': possible access token')
            if re.search(r'-{5}BEGIN [A-Z ]*PRIVATE KEY-{5}',text):errors.append(rel+': possible private key')
        if p.suffix.lower()=='.csv':
            with p.open('r',encoding='utf-8-sig',newline='') as stream:
                reader=csv.reader(stream);header=next(reader);rows=sum(1 for _ in reader)
            sensitive={'participant_id','record_id','sample_id','subject_id','RID','PTID'}&set(header)
            if sensitive:errors.append(rel+': participant-key columns in public CSV')
            if not (rel.startswith('metadata/') or rel.startswith('results/published/')):errors.append(rel+': CSV is outside approved metadata/aggregate folders')
            csv_schemas.append(dict(file=rel,rows=rows,columns=header))
    return dict(status='PASS' if not errors else 'FAIL',files_scanned=scanned,public_csv_tables=len(csv_schemas),csv_inventory=csv_schemas,errors=errors,scope='mechanical file/schema/text audit; runtime outputs/environments/.git excluded; not a semantic data-access authorization review')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--write-manifest',action='store_true');args=parser.parse_args()
    paths=files();report=audit(paths)
    if report['errors']:
        print(json.dumps(report,indent=2));return 1
    if args.write_report:(ROOT/REPORT).write_text(json.dumps(report,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
    if args.write_manifest:
        entries={p.relative_to(ROOT).as_posix():digest(p) for p in files() if p.relative_to(ROOT).as_posix()!=MANIFEST}
        (ROOT/MANIFEST).write_text(json.dumps({'algorithm':'SHA256','files':entries},indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
    if (ROOT/MANIFEST).exists():
        expected=json.loads((ROOT/MANIFEST).read_text('utf-8'))['files']
        observed={p.relative_to(ROOT).as_posix():digest(p) for p in files() if p.relative_to(ROOT).as_posix()!=MANIFEST}
        changed=sorted(k for k in expected.keys()|observed.keys() if expected.get(k)!=observed.get(k))
        if changed:
            print(json.dumps({'status':'FAIL','manifest_mismatches':changed},indent=2));return 1
    print(json.dumps({'status':'PASS','files_scanned':report['files_scanned'],'public_csv_tables':report['public_csv_tables'],'manifest_checked':(ROOT/MANIFEST).exists()}));return 0

if __name__=='__main__':sys.exit(main())
