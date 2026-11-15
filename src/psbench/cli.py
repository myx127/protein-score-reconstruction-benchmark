"""Command-line interfaces; real-data outputs require an external work folder."""
from __future__ import annotations
import argparse,json

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    for command in ['recovery','emory','timing']:
        q=sub.add_parser(command);q.add_argument('--config',required=True);q.add_argument('--private-workdir',required=True)
        q.add_argument('--allow-version-drift',action='store_true')
        if command=='timing':q.add_argument('--repetitions',type=int,default=5)
        else:q.add_argument('--bootstrap',type=int)
    args=p.parse_args()
    if args.command=='recovery':
        from .recovery import run
        result=run(args.config,args.private_workdir,args.bootstrap,args.allow_version_drift)
    elif args.command=='emory':
        from .emory import run
        result=run(args.config,args.private_workdir,args.bootstrap,args.allow_version_drift)
    else:
        from .timing import run
        result=run(args.config,args.private_workdir,args.repetitions,args.allow_version_drift)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
