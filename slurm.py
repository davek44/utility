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

    main_job = Job(cmd, job_name=options.job_name, out_file=options.out_file, err_file=options.err_file, queue=options.queue, cpu=options.cpu, mem=options.mem, time=options.time)
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

        print >> sys.stderr, '%s %s' % (main_job.job_name, main_job.status)

        # delete sbatch
        main_job.clean()


################################################################################
# multi_run
#
#
# Launch and manage multiple SLURM jobs in parallel, using only one 'sacct'
# call per 
################################################################################
def multi_run(jobs, max_proc=None, verbose=False):
    total = len(jobs)
    finished = 0
    running = 0
    active_jobs = []


    while finished + running < total:
        # launch jobs up to the max
        while running < max_proc and finished+running < total:            
            # launch
            jobs[finished+running].launch()
            if verbose:
                print >> sys.stderr, jobs[finished+running].job_name, jobs[finished+running].cmd

            # find it
            time.sleep(5)
            if not jobs[finished+running].update_status():
                time.sleep(10)

            # save it
            active_jobs.append(jobs[finished+running])
            running += 1

        # sleep
        time.sleep(30)

        # update all statuses
        multi_update_status(active_jobs)

        # update active jobs
        active_jobs_new = []
        for i in range(len(active_jobs)):
            if active_jobs[i].status in ['PENDING', 'RUNNING']:
                active_jobs_new.append(active_jobs[i])
            else:
                if verbose:
                    print >> sys.stderr, '%s %s' % (active_jobs[i].job_name, active_jobs[i].status)

                running -= 1
                finished += 1

        active_jobs = active_jobs_new


    # wait for all to finish
    while active_jobs:
        # sleep
        time.sleep(30)

        # update all statuses
        multi_update_status(active_jobs)

        # update active jobs
        active_jobs_new = []
        for i in range(len(active_jobs)):
            if active_jobs[i].status in ['PENDING', 'RUNNING']:
                active_jobs_new.append(active_jobs[i])
            else:
                if verbose:
                    print >> sys.stderr, '%s %s' % (active_jobs[i].job_name, active_jobs[i].status)

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
    ############################################################
    # __init__
    ############################################################
    def __init__(self, cmd, job_name, out_file=None, err_file=None, queue='general', cpu=1, mem=None, time=None):
        self.cmd = cmd
        self.job_name = job_name
        self.out_file = out_file
        self.err_file = err_file
        self.queue = queue
        self.cpu = cpu
        self.mem = mem
        self.time = time

        self.id = None
        self.status = None
        self.sbatch_file = None
        self.sbatch_fd = None


    ############################################################
    # clean
    ############################################################
    def clean(self):
        os.close(self.sbatch_fd)
        os.remove(self.sbatch_file)


    ############################################################
    # launch
    #
    # Make an sbatch file, launch it, and save the job id.
    ############################################################
    def launch(self):
        # make sbatch script
        self.sbatch_fd, self.sbatch_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        sbatch_out = open(self.sbatch_file, 'w')

        print >> sbatch_out, '#!/bin/sh'
        print >> sbatch_out, ''
        print >> sbatch_out, '#SBATCH -p %s' % self.queue
        print >> sbatch_out, '#SBATCH -n %d' % self.cpu
        if self.job_name:
            print >> sbatch_out, '#SBATCH -J %s' % self.job_name
        if self.out_file:
            print >> sbatch_out, '#SBATCH -o %s' % self.out_file
        if self.err_file:
            print >> sbatch_out, '#SBATCH -e %s' % self.err_file
        if self.mem:
            print >> sbatch_out, '#SBATCH --mem %d' % self.mem
        if self.time:
            print >> sbatch_out, '#SBATCH --time %s' % self.time
        print >> sbatch_out, ''
        print >> sbatch_out, self.cmd

        sbatch_out.close()

        # launch it; check_output to get the id
        launch_str = subprocess.check_output('sbatch %s' % self.sbatch_file, shell=True)

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
