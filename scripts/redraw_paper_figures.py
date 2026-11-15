"""Recreate numerical evidence figures using public aggregate tables only."""
import argparse,sys
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from psbench.plots import redraw

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default='outputs/paper_figures');args=p.parse_args()
    for path in redraw(args.output):print(path)
