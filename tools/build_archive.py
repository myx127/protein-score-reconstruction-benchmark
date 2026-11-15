"""Build a verified public ZIP outside the repository, using its manifest."""
import argparse,hashlib,json,sys,zipfile
from pathlib import Path
from verify_release import ROOT,MANIFEST,audit,files,digest

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);args=p.parse_args()
    target=Path(args.output).expanduser().resolve()
    if target==ROOT or ROOT in target.parents:raise ValueError('Place the archive outside the repository')
    if target.suffix.lower()!='.zip':raise ValueError('Output must end in .zip')
    if target.exists():raise FileExistsError('Archive already exists; choose a new destination')
    if audit(files())['errors']:raise ValueError('Public release audit failed')
    expected=json.loads((ROOT/MANIFEST).read_text('utf-8'))['files']
    actual={q.relative_to(ROOT).as_posix():digest(q) for q in files() if q.name!=MANIFEST}
    if actual!=expected:raise ValueError('Manifest does not match; run verify_release.py first')
    target.parent.mkdir(parents=True,exist_ok=True)
    prefix='protein-score-reconstruction-benchmark-v1.0.0'
    with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name in sorted([*expected,MANIFEST]):archive.write(ROOT/name,prefix+'/'+name)
    checksum=hashlib.sha256(target.read_bytes()).hexdigest()
    target.with_suffix('.zip.sha256').write_text(checksum+'  '+target.name+'\n',encoding='ascii')
    print(json.dumps({'status':'PASS','archive':str(target),'files':len(expected)+1,'bytes':target.stat().st_size,'SHA256':checksum},ensure_ascii=True))

if __name__=='__main__':main()
