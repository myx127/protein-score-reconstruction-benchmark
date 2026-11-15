"""Run all scoring methods on independent toy data without HMACS access."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from psbench.synthetic import generate
from psbench.recovery import run as recovery
from psbench.emory import run as emory

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default='outputs/synthetic_demo');p.add_argument('--bootstrap',type=int,default=50);p.add_argument('--skip-emory',action='store_true');p.add_argument('--with-timing',action='store_true');args=p.parse_args()
    out=Path(args.output).resolve();r,e=generate(out/'inputs',n_bootstrap=args.bootstrap)
    status={'data_kind':'synthetic','recovery':recovery(r,out/'recovery',n_boot=args.bootstrap)}
    if not args.skip_emory:status['emory']=emory(e,out/'emory',n_boot=args.bootstrap)
    if args.with_timing:
        from psbench.timing import run
        status['timing']=run(r,out/'timing',repetitions=5)
    (out/'DEMO_STATUS.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(status,indent=2))

if __name__=='__main__':main()
