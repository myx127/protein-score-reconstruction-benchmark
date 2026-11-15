"""Check CRNP projection invariants with synthetic data; requires R."""
import argparse,sys
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from psbench.io import run_r_job

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default='outputs/r_contracts');args=p.parse_args()
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
    run_r_job('contracts',{'output':str(out),'strict_environment':False},out)
    print((out/'R_CONTRACTS.json').read_text('utf-8'))
