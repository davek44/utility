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
    parser.add_option('-o', dest='out_file')
    parser.add_option('-e', dest='err_file')
    parser.add_option('-q', dest='queue', default='general')
    parser.add_option('-c', dest='cpu', default=1, type='int')
    parser.add_option('-m', dest='mem', default=None, type='int')    
    (options,args) = parser.parse_args()

    cmd = args[0]

    main_job = Job(cmd, out_file=options.out_file, err_file=options.err_file, queue=options.queue, cpu=options.cpu, mem=options.mem)
    main_job.launch()

    time.sleep(10)

    while main_job.update_status() and main_job.status in ['PENDING','RUNNING']:
        time.sleep(10)

    main_job.clean()
    

class Job:
    ############################################################
    # __init__
    ############################################################
    def __init__(self, cmd, out_file=None, err_file=None, queue='general', cpu=1, mem=None):
        self.cmd = cmd
        self.out_file = out_file
        self.err_file = err_file
        self.queue = queue
        self.cpu = cpu
        self.mem = mem
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
        if self.out_file:
            print >> sbatch_out, '#SBATCH -o %s' % self.out_file
        if self.err_file:
            print >> sbatch_out, '#SBATCH -e %s' % self.err_file
        if self.mem:
            print >> sbatch_out, '#SBATCH -mem %d' % self.mem
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
    # Use 'sacct' tp update the job's status. Return True if
    # found and False if not.
    ############################################################
    def update_status(self):
        status = None

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
