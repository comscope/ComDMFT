#!/usr/bin/env python
"""
A module to read the joblist file. The joblist file contains one test job
per line. Additionally the file may contain comments that start with a '#' 
character.
"""
import re

def read_joblist(filelist):
    """
    Read all the files in the filelist, strip all the comments, and return
    the list of test jobs. Each test job is identified by a string.
    """
    joblist = []
    for filename in filelist:
        fp = open(filename,"r")
        lines = fp.readlines()
        fp.close()
        for line in lines:
            mobj = re.search("#",line)
            if mobj:
                # Remove everything after a "#" character
                line = line[:mobj.start()]
            # Delete all leading and trailing whitespace
            line = line.strip()
            # If there is anything left add it to the joblist
            if len(line) > 0:
                joblist.append(line)
    return joblist

def parse_arguments():
    """
    Parse command line arguments for the test_joblist.py script when run
    from the command line. In this case each argument is interpreted as 
    filename of a file containing job names. The list of names is passed to
    read_joblist and the resulting joblist is printed.
    """
    from argparse import ArgumentParser
    prs = ArgumentParser(description="Generate a joblist from the jobs listed in the specified files.")
    prs.add_argument("files",help="joblist files",type=str,nargs='+')
    args = prs.parse_args()
    return args

if __name__ == "__main__":
    filelist = parse_arguments()
    joblist = read_joblist(filelist.files)
    for job in joblist:
        print(job)
