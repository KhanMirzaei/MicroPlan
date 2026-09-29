import argparse
import json
from . import __version__
from .planner import run


def main(argv=None):
    parser=argparse.ArgumentParser(prog='microplan',description='Microbiome processing allocation and design checks; no power estimates or treatment assignment.')
    parser.add_argument('--version',action='version',version=__version__)
    sub=parser.add_subparsers(dest='command',required=True)
    command=sub.add_parser('plan',help='Generate sample sheets, design checks, cost and attrition reports')
    command.add_argument('--config',required=True,help='JSON study configuration')
    command.add_argument('--out',required=True,help='New or empty output directory')
    args=parser.parse_args(argv)
    try:print(json.dumps(run(args.config,args.out),indent=2))
    except (ValueError,KeyError,OSError,TypeError) as e:parser.exit(2,f'microplan: {e}\n')
