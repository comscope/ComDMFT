#!/usr/bin/env python
"""
A module to deal with batch job systems.

Important functions include:
- detect what, if any, batch system is installed
- submit jobs
- check job status
If no batch system is installed job submission will just run the calculations
interactively. 
"""
# We will use a number of system calls to drive this thing
import subprocess
import re
import sys
from argparse import ArgumentParser

class batch:
    """
    The batch class provides methods to interact with the batch system.
    In addition some information on the type of batch system is kept
    so that methods can make decisions efficiently.
    """
    # Constant identifiers for the various batch systems
    batch_none  = 0
    batch_pbs   = 1
    batch_sge   = 2
    batch_slurm = 3
    # Constant identifiers for the various job states
    batch_running    = 100
    batch_terminated = 101

    def __init__(self):
        """
        Initialize an instance of the batch class which currently simply
        means to state that no batch system has been detected yet.
        """
        self.batch_actual = None

    def guess(self):
        """
        Guess which batch system is available and set the corresponding
        variable to the associated identifier.
        """
        try:
            subprocess.check_call(["which","sbatch"])
            self.batch_actual = self.batch_slurm
        except subprocess.CalledProcessError:
            # the sbatch command was not found
            pass

        try:
            subprocess.check_call(["which","qsub"])
            self.batch_actual = self.batch_pbs
        except subprocess.CalledProcessError:
            # the sbatch command was not found
            pass

        if not self.batch_actual:
            self.batch_actual = self.batch_none


    def submit(self,jobScript):
        """
        Submit the named job script to the job queue and return the job
        identifier. If no batch system is available run the job 
        interactively and return None.
        """
        cmd = ""
        if self.batch_actual == self.batch_slurm:
            cmd = "sbatch"
        elif self.batch_actual == self.batch_pbs:
            cmd = "qsub"
        if self.batch_actual == self.batch_none:
            try:
                subprocess.check_call([jobScript])
            except subprocess.CalledProcessError:
                pass
            rtn = None
        else:
            rtn = subprocess.check_output([cmd,jobScript])
            if sys.version_info.major==3:
                rtn = rtn.decode()
            rtn = rtn.replace("\n","")
        return rtn

    def status(self,jobList):
        """
        Given a list of job identifiers check whether any of them are still
        in the queue. A list of jobs is used to avoid hammering the batch
        system with job queries.

        Return batch_running if some jobs are still in the queue, and
        return batch_terminated otherwise.

        If there is no batch system this function will always return
        batch_terminated.
        """
        cmd = ""
        args = []
        if self.batch_actual == self.batch_slurm:
            cmd = "squeue"
        elif self.batch_actual == self.batch_pbs:
            cmd = "qstat"
        if self.batch_actual == self.batch_none:
            rtn = self.batch_terminated
        else:
            rtn = self.batch_terminated
            for job in jobList:
                try:
                    subprocess.check_call([cmd]+[job])
                    rtn = self.batch_running
                except subprocess.CalledProcessError:
                    pass
        return rtn

def parse_arguments():
    """
    Parse command line arguments for the test_batch.py script when run from
    the command line. In this case each argument is interpreted as a command
    and hence the same argument may be specified multiple times and it is
    important that the list of arguments returned preserves the order of 
    the arguments specified. This requires using special features of the
    argparse module.
    """
    prs = ArgumentParser(description="Interact with the batch job system by issueing a sequence of commands.")
    prs.add_argument("-q","--queue",help="queue a job",type=str,nargs='*',action="append")
    prs.add_argument("-s","--status",help="check whether any of the queued jobs are still running STATUS times",nargs=1,type=int,default=[5])
    prs.add_argument("-w","--wait",help="wait WAIT seconds",type=int,nargs=1,action="store",default=[10])
    args = prs.parse_args()
    return args

if __name__ == "__main__":
    import time
    jobs = batch()
    jobs.guess()
    if jobs.batch_actual == jobs.batch_slurm:
        print("Batch system looks like SLURM")
    elif jobs.batch_actual == jobs.batch_pbs:
        print("Batch system looks like PBS")
    elif jobs.batch_actual == jobs.batch_sge:
        print("Batch system looks like SGE")
    elif jobs.batch_actual == jobs.batch_none:
        print("Batch system looks like interactive")
    else:
        print("Batch system looks like we have a problem")
    list = []
    args = parse_arguments()
    print(args)
    print(vars(args))
    print(type(args.queue[0]))
    if type(args.queue[0]) is type([0,1]):
        queue = args.queue[0]
    else:
        queue = args.queue
    for job in queue:
        list.append(jobs.submit(job))
    print(list)
    for ii in range(0,args.status[0]):
        time.sleep(args.wait[0])
        status = jobs.status(list)
        if status == jobs.batch_running:
            print("Batch jobs still running")
        elif status == jobs.batch_terminated:
            print("Batch jobs still terminated")
        else:
            print("Batch jobs are messed up")
