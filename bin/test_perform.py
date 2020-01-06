#!/usr/bin/env python
"""
A script to perform a set of tests on the MatDeLab code.

The script takes a set of filenames where each of the files lists test case
directories. In each test case directory there is a script 'sub.sh' that runs 
the test. This script has job submission parameters as well as the command
to run. Hence it can be submitted to a job queue, if there is one, or run
interactively, if there is no job queue.

The script submits (or runs) all jobs and waits for them to finish. After
that the results are checked and a report is generated.
"""

from argparse import ArgumentParser
import test_joblist
import test_batch
import os
import sys
import time
import subprocess

def parse_arguments():
    """
    Parse command line arguments for the test_perform.py script. In this case
    each argument is interpreted as filename of a file containing job names.
    """
    prs = ArgumentParser(description="Generate a joblist from the jobs listed in the specified files.")
    prs.add_argument("-v","--validate",dest="validate",
        help="Validate the current set of outputs only, do not run the calculations.",
        const=True,default=False,action='store_const')
    prs.add_argument("files",help="joblist files",type=str,nargs='+')
    args = prs.parse_args()
    return args

filelist = parse_arguments()
testlist = test_joblist.read_joblist(filelist.files)
joblist  = []
if not filelist.validate:
    queue    = test_batch.batch()
    queue.guess()
    for test in testlist:
        os.chdir(test)
        joblist.append(queue.submit("./sub.sh"))
        os.chdir("..")

    status = queue.batch_running
    while status == queue.batch_running:
        time.sleep(300)
        status = queue.status(joblist)

exit_code = 0
toldiff = []
toldiff.append("toldiff.py")
toldiff.append("--summary")
toldiff.append("OK:OK:Fail")
toldiff.append("--output")
toldiff.append("summary")
toldiff.append("--exit")
toldiff.append("0:0:0")
toldiff.append("ref")
toldiff.append("log")
for test in testlist:
    os.chdir(test)
    sys.stdout.write(test+" ... ")
    result = subprocess.check_output(toldiff)
    if sys.version_info.major==3:
        result = result.decode()
    if result != "OK":
        exit_code = 1
    sys.stdout.write(result+"\n")
    os.chdir("..")
sys.exit(exit_code)
