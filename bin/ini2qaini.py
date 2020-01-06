#!/usr/bin/env python
"""
Transform an input file to a QA job
===================================

A Python script to take a DMFT MatDeLab input file and modify it to
one that is suitable for regression testing. That means the resulting 
input file will run quickly and allows one to test the mathematical 
"correctness" and the numerical stability of the code. However, physically
the results must be expected to be INACCURATE to the point of being useless.

The approach in this script is trivially simple. We read the whole input 
file. Then we process it line by line in that we either echo the input line
unmodified to the output, or for selected input lines we emit lines with
alternative parameter settings to the output.
"""
import sys
import re

versionstr = '%(prog)s version 0.0'

# Every line that contains the key from the dictionary is replaced with the
# line associated with this key. 
#
replace_with = { "CONTROL"  : "CONTROL   iter_dft=%3d  iter_hf=%3d  iter_gw=%3d iter_qp=%3d",
                 "nproc_tau": "          nproc_tau=  12 nproc_k=   2" }

def parse_arguments():
    """
    Parse the command line arguments for the ini2qaini.py script when run
    from the command line.
    """
    from argparse import ArgumentParser
    prs = ArgumentParser()
    prs.add_argument("--version",action='version',
                     version=versionstr)
    prs.add_argument("--outputfile","-o",
                     help="A MatDeLab QA quality input file")
    prs.add_argument("--dft-it",dest="dftit",default=80,type=int,
                     help="The number of GW iterations")
    prs.add_argument("--hf-it",dest="hfit",default=0,type=int,
                     help="The number of GW iterations")
    prs.add_argument("--gw-it",dest="gwit",default=4,type=int,
                     help="The number of GW iterations")
    prs.add_argument("--qp-it",dest="qpit",default=4,type=int,
                     help="The number of GW iterations")
    prs.add_argument("inputfile",help="A regular MatDeLab input file")
    args = prs.parse_args()
    if not args.outputfile:
        tmp = args.inputfile.rsplit(".",1)
        tmp[0] = tmp[0]+"-qa"
        if len(tmp) == 1:
            args.outputfile = tmp[0]
        elif len(tmp) == 2:
            args.outputfile = tmp[0]+"."+tmp[1]
        else:
            print("ERROR: something bad happened")
            print("inputfile=",args.inputfile)
            print("tmp      =",tmp)
    return args

def execute_with_arguments(args):
    f = open(args.inputfile,'r')
    lines = f.readlines()
    f.close()
    f = open(args.outputfile,'w')
    keys = replace_with.keys()
    for line in lines:
        for key in keys:
            if re.search(key,line):
                if key == "CONTROL":
                    line = (replace_with[key] % (args.dftit,args.hfit,args.gwit,args.qpit))+"\n"
                else:
                    line = replace_with[key]+"\n"
        f.write(line)
    f.close()

def main():
    execute_with_arguments(parse_arguments())

if __name__ == "__main__":
    main()
