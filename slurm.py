#!/usr/bin/env python
from optparse import OptionParser
import os, sys, subprocess, tempfile, time

################################################################################
# slurm.py
#
# Methods to run jobs on SLURM.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='go', default=False, action='store_true', help='Don\'t wait for the job to finish [Default: %default]')

    parser.add_option('-o', dest='out_file')
    parser.add_option('-e', dest='err_file')

    parser.add_option('-J', dest='job_name')

    parser.add_option('-q', dest='queue', default='general')
    parser.add_option('-n', dest='cpu', default=1, type='int')
    parser.add_option('-m', dest='mem', default=None, type='int')
    parser.add_option('-t', dest='time', default=None)

    (options,args) = parser.parse_args()

    cmd = args[0]

    main_job = Job(cmd, name=options.job_name, out_file=options.out_file, err_file=options.err_file, queue=options.queue, cpu=options.cpu, mem=options.mem, time=options.time)
    main_job.launch()

    if options.go:
        time.sleep(1)

        # find the job
        if not main_job.update_status:
            time.sleep(1)

        # delete sbatch
        main_job.clean()

    else:
        time.sleep(10)

        # find the job
        if not main_job.update_status():
            time.sleep(10)

        # wait for it to complete
        while main_job.update_status() and main_job.status in ['PENDING','RUNNING']:
            time.sleep(30)

        print('%s %s' % (main_job.name, main_job.status), file=sys.stderr)

        # delete sbatch
        main_job.clean()


################################################################################
# multi_run
#
#
# Launch and manage multiple SLURM jobs in parallel, using only one 'sacct'
# call per
################################################################################
def multi_run(jobs, max_proc=None, verbose=False, sleep_time=30):
    total = len(jobs)
    finished = 0
    running = 0
    active_jobs = []

    if max_proc is None:
        max_proc = len(jobs)

    while finished + running < total:
        # launch jobs up to the max
        while running < max_proc and finished+running < total:
            # launch
            jobs[finished+running].launch()
            if verbose:
                print(jobs[finished+running].name, jobs[finished+running].cmd, file=sys.stderr)

            # find it
            time.sleep(5)
            if not jobs[finished+running].update_status():
                time.sleep(10)

            # save it
            active_jobs.append(jobs[finished+running])
            running += 1

        # sleep
        time.sleep(sleep_time)

        # update all statuses
        multi_update_status(active_jobs)

        # update active jobs
        active_jobs_new = []
        for i in range(len(active_jobs)):
            if active_jobs[i].status in ['PENDING', 'RUNNING']:
                active_jobs_new.append(active_jobs[i])
            else:
                if verbose:
                    print('%s %s' % (active_jobs[i].name, active_jobs[i].status), file=sys.stderr)

                running -= 1
                finished += 1

        active_jobs = active_jobs_new


    # wait for all to finish
    while active_jobs:
        # sleep
        time.sleep(sleep_time)

        # update all statuses
        multi_update_status(active_jobs)

        # update active jobs
        active_jobs_new = []
        for i in range(len(active_jobs)):
            if active_jobs[i].status in ['PENDING', 'RUNNING']:
                active_jobs_new.append(active_jobs[i])
            else:
                if verbose:
                    print('%s %s' % (active_jobs[i].name, active_jobs[i].status), file=sys.stderr)

                running -= 1
                finished += 1

        active_jobs = active_jobs_new


################################################################################
# multi_update_status
#
# Update the status for multiple jobs at once.
################################################################################
def multi_update_status(jobs):
    # reset all
    for j in jobs:
        j.status = None

    # try multiple times because sometimes it fails
    attempt = 0
    while attempt < 3 and [j for j in jobs if j.status == None]:
        if attempt > 0:
            time.sleep(10)

        sacct_str = subprocess.check_output('sacct', shell=True)
        sacct_str = sacct_str.decode('UTF-8')

        # split into job lines
        sacct_lines = sacct_str.split('\n')
        for line in sacct_lines[2:]:
            a = line.split()

            try:
                line_id = int(a[0])
            except:
                line_id = None

            # check call jobs for a match
            for j in jobs:
                if line_id == j.id:
                    j.status = a[5]

        attempt += 1


class Job:
    ''' class to manage SLURM jobs.

    Notes:
     -Since we have two types of machines in the GPU queue, I'm asking
      for the machine type as "queue", and the "launch" method will handle it.
    '''
    ############################################################
    # __init__
    ############################################################
    def __init__(self, cmd, name, out_file=None, err_file=None, queue='general', cpu=1, mem=None, time=None, gpu=0):
        self.cmd = cmd
        self.name = name
        self.out_file = out_file
        self.err_file = err_file
        self.queue = queue
        self.cpu = cpu
        self.mem = mem
        self.time = time
        self.gpu = gpu

        self.id = None
        self.status = None


    ############################################################
    # launch
    #
    # Make an sbatch file, launch it, and save the job id.
    ############################################################
    def launch(self):
        # make sbatch script
        sbatch_tempf = tempfile.NamedTemporaryFile()
        sbatch_out = open(sbatch_tempf.name, 'w')

        print('#!/bin/sh\n', file=sbatch_out)
        if self.gpu > 0:
            print('#SBATCH -p gpu', file=sbatch_out)
            print('#SBATCH --gres=gpu:%s:%d\n' % (self.queue, self.gpu), file=sbatch_out)
        else:
            print('#SBATCH -p %s' % self.queue, file=sbatch_out)
        print('#SBATCH -n 1', file=sbatch_out)
        print('#SBATCH -c %d' % self.cpu, file=sbatch_out)
        if self.name:
            print('#SBATCH -J %s' % self.name, file=sbatch_out)
        if self.out_file:
            print('#SBATCH -o %s' % self.out_file, file=sbatch_out)
        if self.err_file:
            print('#SBATCH -e %s' % self.err_file, file=sbatch_out)
        if self.mem:
            print('#SBATCH --mem %d' % self.mem, file=sbatch_out)
        if self.time:
            print('#SBATCH --time %s' % self.time, file=sbatch_out)
        print(self.cmd, file=sbatch_out)

        sbatch_out.close()

        # launch it; check_output to get the id
        launch_str = subprocess.check_output('sbatch %s' % sbatch_tempf.name, shell=True)

        # e.g. "Submitted batch job 13861989"
        self.id = int(launch_str.split()[3])


    ############################################################
    # update_status
    #
    # Use 'sacct' to update the job's status. Return True if
    # found and False if not.
    ############################################################
    def update_status(self):
        status = None

        attempt = 0
        while attempt < 3 and status == None:
            if attempt > 0:
                time.sleep(10)

            sacct_str = subprocess.check_output('sacct', shell=True)
            sacct_str = sacct_str.decode('UTF-8')

            sacct_lines = sacct_str.split('\n')
            for line in sacct_lines[2:]:
                a = line.split()

                try:
                    line_id = int(a[0])
                except:
                    line_id = None

                if line_id == self.id:
                    status = a[5]

            attempt += 1

        if status == None:
            return False
        else:
            self.status = status
            return True


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
