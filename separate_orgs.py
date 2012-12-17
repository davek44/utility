#!/usr/bin/env python
# encoding: utf-8
"""
separate_orgs.py

Created by Cole Trapnell on 2012-05-06.
Copyright (c) 2012 Cole Trapnell. All rights reserved.
"""

import sys
import getopt
import subprocess
import errno
import os
import warnings
import re
import signal
from datetime import datetime, date, time
from shutil import copy
import logging

help_message = '''
separate_orgs.py <input_file1.bam> <input_file2.bam>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

output_dir = "./separate_orgs_out/"
logging_dir = output_dir + "logs/"
run_log = None
tophat_log = None  #main log file handle
tophat_logger = None # main logging object
run_cmd = None
tmp_dir = output_dir + "tmp/"
bin_dir = sys.path[0] + "/"
use_zpacker = False # this is set by -z/--zpacker option (-z0 leaves it False)

use_BAM_Unmapped = False # automatically set to True for non-Solid reads, handles unmapped reads in BAM format

use_BWT_FIFO = False # can only be set to True if use_zpacker is True and only with -C/--color
# enabled by -X/-unmapped-fifo option (unless -z0)
unmapped_reads_fifo = None # if use_BWT_FIFO is True, this tricks bowtie into writing the
                           # unmapped reads into a compressed file

samtools_path = None
bowtie_path = None
fail_str = "\t[FAILED]\n"
gtf_juncs = None #file name with junctions extracted from given GFF file

def init_logger(log_fname):
    global tophat_logger
    tophat_logger = logging.getLogger('project')
    formatter = logging.Formatter('%(asctime)s %(message)s', '[%Y-%m-%d %H:%M:%S]')
    # level = logging.__dict__.get(options.loglevel.upper(),logging.DEBUG)
    tophat_logger.setLevel(logging.DEBUG)

    hstream = logging.StreamHandler(sys.stderr)
    hstream.setFormatter(formatter)
    tophat_logger.addHandler(hstream)
    #
    # Output logging information to file
    if os.path.isfile(log_fname):
        os.remove(log_fname)
    global tophat_log
    logfh = logging.FileHandler(log_fname)
    logfh.setFormatter(formatter)
    tophat_logger.addHandler(logfh)
    tophat_log=logfh.stream

def nonzeroFile(filepath):
  if os.path.exists(filepath):
     fpath, fname=os.path.split(filepath)
     fbase, fext =os.path.splitext(fname)
     if fext.lower() == ".bam":
        output = os.popen("samtools view %s | head | wc -l" % filepath).read()[:-1]
        if int(output) > 0:
          return True
     else:
        if os.path.getsize(filepath)>25:
          return True
  return False

def bamExists_and_NonEmpty(filepath):
  if os.path.exists(filepath):
      output = os.popen("samtools view %s | head | wc -l" % filepath).read()[:-1]
      if int(output) > 0:
          return True

  return False

# check if a file exists and has non-zero (or minimum) size
def fileExists(filepath, minfsize=1):
  if os.path.exists(filepath) and os.path.getsize(filepath)>=minfsize:
     return True
  else:
     return False


def getFileDir(filepath):
   #if fullpath given, returns path including the ending /
   fpath, fname=os.path.split(filepath)
   if fpath: fpath+='/'
   return fpath

def getFileBaseName(filepath):
   fpath, fname=os.path.split(filepath)
   fbase, fext =os.path.splitext(fname)
   fx=fext.lower()
   if (fx in ['.fq','.txt','.seq','.bwtout'] or fx.find('.fa')==0) and len(fbase)>0:
      return fbase
   elif fx == '.z' or fx.find('.gz')==0 or fx.find('.bz')==0:
      fb, fext = os.path.splitext(fbase)
      fx=fext.lower()
      if (fx in ['.fq','.txt','.seq','.bwtout'] or fx.find('.fa')==0) and len(fb)>0:
         return fb
      else:
         return fbase
   else:
     if len(fbase)>0:
        return fbase
     else:
        return fname

# Returns the current time in a nice format
def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

# The TopHat logging formatter
def th_log(out_str):
  #print >> sys.stderr, "[%s] %s" % (right_now(), out_str)
  if tophat_logger:
       tophat_logger.info(out_str)

def th_logp(out_str=""):
  print >> sys.stderr, out_str
  if tophat_log:
        print >> tophat_log, out_str

def die(msg=None):
  if msg is not None:
    th_logp(msg)
  sys.exit(1)

# Ensures that the output, logging, and temp directories are present. If not,
# they are created
def prepare_output_dir():
    global output_dir
    global logging_dir
    global tmp_dir
    #th_log("Preparing output location "+output_dir)
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)

    if os.path.exists(logging_dir):
        pass
    else:
        os.mkdir(logging_dir)

    if os.path.exists(tmp_dir):
        pass
    else:
        try:
          os.makedirs(tmp_dir)
        except OSError, o:
          die("\nError creating directory %s (%s)" % (tmp_dir, o))


# to be added as preexec_fn for every subprocess.Popen() call:
# see http://bugs.python.org/issue1652
def subprocess_setup():
 # Python installs a SIGPIPE handler by default, which causes
 # gzip or other de/compression pipes to complain about "stdout: Broken pipe"
   signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Retrive a tuple containing the system's version of samtools.  Parsed from
# `samtools`
def get_samtools_version():
   try:
       # Launch Bowtie to capture its version info
       proc = subprocess.Popen(samtools_path, stderr=subprocess.PIPE)
       samtools_out = proc.communicate()[1]

       # Find the version identifier
       version_match = re.search(r'Version:\s+(\d+)\.(\d+).(\d+)([a-zA-Z]?)', samtools_out)
       samtools_version_arr = [int(version_match.group(x)) for x in [1,2,3]]
       if version_match.group(4):
           samtools_version_arr.append(version_match.group(4))
       else:
           samtools_version_arr.append(0)

       return version_match.group(), samtools_version_arr
   except OSError, o:
      errmsg=fail_str+str(o)+"\n"
      if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
          errmsg+="Error: samtools not found on this system"
      die(errmsg)

# Make sure the SAM tools are installed and are recent enough to be useful
def check_samtools():
   th_log("Checking for Samtools")
   global samtools_path
   samtools_path=prog_path("samtools")
   samtools_version_str, samtools_version_arr = get_samtools_version()
   if samtools_version_str == None:
       die("Error: Samtools not found on this system")
   elif  samtools_version_arr[1] < 1 or samtools_version_arr[2] < 7:
       die("Error: TopHat requires Samtools 0.1.7 or later")
   th_logp("\t\tSamtools version:\t %s" % ".".join([str(x) for x in samtools_version_arr]))

# Format a DateTime as a pretty string.
# FIXME: Currently doesn't support days!
def formatTD(td):
   days = td.days
   hours = td.seconds // 3600
   minutes = (td.seconds % 3600) // 60
   seconds = td.seconds % 60

   if days > 0:
       return '%d days %02d:%02d:%02d' % (days, hours, minutes, seconds)
   else:
       return '%02d:%02d:%02d' % (hours, minutes, seconds)


# rough equivalent of the 'which' command to find external programs
# (current script path is tested first, then PATH envvar)
def which(program):
   def is_executable(fpath):
       return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
   fpath, fname = os.path.split(program)
   if fpath:
       if is_executable(program):
           return program
   else:
       progpath = os.path.join(bin_dir, program)
       if is_executable(progpath):
          return progpath
       for path in os.environ["PATH"].split(os.pathsep):
          progpath = os.path.join(path, program)
          if is_executable(progpath):
             return progpath
   return None

def prog_path(program):
   progpath=which(program)
   if progpath == None:
       die("Error locating program: "+program)
   return progpath

# FIXME: this should get set during the make dist autotools phase of the build
def get_version():
  return "2.0.0"

def mlog(msg):
 print >> sys.stderr, "[DBGLOG]:"+msg

def test_input_file(filename):
   try:
       test_file = open(filename, "r")
   except IOError, o:
       die("Error: Opening file %s" % filename)
   return

def sort_by_read(input_name, sorted_out_name=None):
    if sorted_out_name == None:
        sorted_out_name = tmp_dir + "%s_sortedbyread" % input_name  
        
    bamsort_cmd = ["sort",
                   "-k",
                   "1,1",
                   input_name]
                          
    th_logp("  Executing: " + " ".join(bamsort_cmd))
    ret = 0
    
    sam_file = open(sorted_out_name, "w")
    ret = subprocess.call(bamsort_cmd, stdout=sam_file)
    if ret != 0:
        die("Could not sort file by read")
    return sorted_out_name
    
def sort_by_coord(input_name, sorted_out_name=None):
    if sorted_out_name == None:
        sorted_out_name = tmp_dir + "%s_sorted" % input_name  

    bamsort_cmd = [samtools_path,
                   "sort",
                   input_name,
                   sorted_out_name]

    th_logp("  Executing: " + " ".join(bamsort_cmd))
    ret = 0
    ret = subprocess.call(bamsort_cmd)
    if ret != 0:
        die("Could not sort file by coord")
    sorted_out_name += ".bam"
    return sorted_out_name

def convert_bam_to_sam(bam_filename, sam_filename=None):
    if sam_filename == None:
        sam_filename = "%s.sam" % bam_filename  
      
    header = []
    
    bamtosam_cmd = [samtools_path,
                   "view",
                   "-H",
                   bam_filename]
                          
    th_logp("  Executing: " + " ".join(bamtosam_cmd))
    ret = 0
    
    sam_file = open(tmp_dir + "/header.txt", "w")
    
    ret = subprocess.call(bamtosam_cmd, stdout=sam_file)
    if ret != 0:
        die("Could not convert file to sam")
    
    sam_file = open(tmp_dir + "/header.txt", "r")
    header = get_header(sam_file)
      
    bamtosam_cmd = [samtools_path,
                   "view",
                   bam_filename]
                          
    th_logp("  Executing: " + " ".join(bamtosam_cmd))
    ret = 0
    
    sam_file = open(sam_filename, "w")
    
    ret = subprocess.call(bamtosam_cmd, stdout=sam_file)
    if ret != 0:
        die("Could not convert file to sam")

    return (sam_filename, header)
    
def convert_sam_to_bam(sam_filename, bam_filename=None):
    if bam_filename == None:
        bam_filename = "%s.sam" % sam_filename  

    samtobam_cmd = [samtools_path,
                   "view",
                   "-b",
                   "-h",
                   "-S",
                   sam_filename]

    th_logp("  Executing: " + " ".join(samtobam_cmd))
    ret = 0

    bam_file = open(bam_filename, "wb")

    ret = subprocess.call(samtobam_cmd, stdout=bam_file)
    if ret != 0:
        die("Could not convert file to bam")

    return bam_filename

def get_header(fp):
    header_lines = []
    last_pos = fp.tell()
    line = fp.readline()
    while line != '':
        if line[0] != '@':
            fp.seek(last_pos)
            break
        line = line.strip()
        header_lines.append(line)
        last_pos = fp.tell()
        line = fp.readline()
    return header_lines
    
def get_next_align_group(sam_file):
    alignments = []
    last_pos = sam_file.tell()
    line = sam_file.readline()
    if line == '':
        return alignments
    cols = line.split()
    read_name = cols[0]
    while line != '':
        cols = line.split()
        if cols[0] != read_name:
            sam_file.seek(last_pos)
            break
        alignments.append(cols)
        last_pos = sam_file.tell()
        line = sam_file.readline()
    #print alignments
    return alignments

def print_alignments(alignments, fp_out):
    for al_cols in alignments:
        print >> fp_out, "\t".join(al_cols)

def append_header(header_lines, fp_out):
    for line in header_lines:
        print >> fp_out, line

def rank_alignments(alignments):
    best_ranking = [99,99, False]
    for al_cols in alignments:
        if al_cols[6] == "=":
            best_ranking[2] = True # alignment is paired
        if len(al_cols) >= 12:
            for i in range(11, len(al_cols)):
                attr_cols = al_cols[i].split(":")
                flag = int(al_cols[1])
                if flag & 0x40:
                    mate_idx = 0
                elif flag & 0x80:
                    mate_idx = 1
                else:
                    #print "warning: bad mate flag"
                    # ... or unmated reads
                    mate_idx = 0
                    
                if len(attr_cols) == 3 and attr_cols[0] == "NM" and int(attr_cols[2]) < best_ranking[mate_idx]:
                    best_ranking[mate_idx] = int(attr_cols[2])
    return best_ranking

def separate_streams(left_sam_sortedbyread_filename, 
                     right_sam_sortedbyread_filename,
                     left_header,
                     right_header,
                     left_org_sam_sortedbyread_filename,
                     right_org_sam_sortedbyread_filename):
    left_file = open(left_sam_sortedbyread_filename, "r")
    right_file = open(right_sam_sortedbyread_filename, "r")

    left_alignments = get_next_align_group(left_file)
    right_alignments = get_next_align_group(right_file)

    #print left_alignments
    #print right_alignments
    
    left_org_sam_sortedbyread_file = open(left_org_sam_sortedbyread_filename,"w")
    right_org_sam_sortedbyread_file = open(right_org_sam_sortedbyread_filename,"w")

    append_header(left_header, left_org_sam_sortedbyread_file)
    append_header(right_header, right_org_sam_sortedbyread_file)
    
    suppressed_left = 0
    suppressed_right = 0
    
    while True:
        if left_alignments == []:
            if right_alignments == []:
                break # We're done with both streams 
            else:
                print_alignments(right_alignments, right_org_sam_sortedbyread_file)
                right_alignments = get_next_align_group(right_file)
        elif right_alignments == []:
            if left_alignments == []:
                break # We're done with both streams
            else:
                print_alignments(left_alignments, left_org_sam_sortedbyread_file)
                left_alignments = get_next_align_group(left_file)
        else:
            left_read_id = left_alignments[0][0]
            right_read_id = right_alignments[0][0]
            if left_read_id < right_read_id:
                # Just advance the left stream
                print_alignments(left_alignments, left_org_sam_sortedbyread_file)
                left_alignments = get_next_align_group(left_file)
            elif right_read_id < left_read_id:
                # Just advance the right stream
                print_alignments(right_alignments, right_org_sam_sortedbyread_file)
                right_alignments = get_next_align_group(right_file)
            else:
                #print_alignments(left_alignments, left_org_sam_sortedbyread_file)
                #print_alignments(right_alignments, right_org_sam_sortedbyread_file)
                left_rank = rank_alignments(left_alignments)
                right_rank = rank_alignments(right_alignments)
                
                if left_rank[2] == False and right_rank[2] == True:
                    suppressed_left += len(left_alignments)
                    print_alignments(right_alignments, right_org_sam_sortedbyread_file)
                elif left_rank[2] == True and right_rank[2] == False:
                    suppressed_right += len(right_alignments)
                    print_alignments(left_alignments, left_org_sam_sortedbyread_file)
                else: 
                    # Both are paired (or unpaired), pick the one with fewer mismatches
                    if left_rank[0] + left_rank[1] > right_rank[0] + right_rank[1]:
                        suppressed_left += len(left_alignments)
                        print_alignments(right_alignments, right_org_sam_sortedbyread_file)
                    elif left_rank[0] + left_rank[1] < right_rank[0] + right_rank[1]:
                        suppressed_right += len(right_alignments)
                        print_alignments(left_alignments, left_org_sam_sortedbyread_file)
                    else:
                        #print "left is better"
                        print_alignments(left_alignments, left_org_sam_sortedbyread_file)
                        print_alignments(right_alignments, right_org_sam_sortedbyread_file)  
                            
                # Advance both streams                
                left_alignments = get_next_align_group(left_file)
                right_alignments = get_next_align_group(right_file)
    print "Suppressed %d alignments from left stream" % suppressed_left
    print "Suppressed %d alignments from right stream" % suppressed_right
                
def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
    
        global output_dir
        global logging_dir
        global tmp_dir
        
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                global output_dir
                output_dir = value + "/"
                logging_dir = output_dir + "/logs/"
                tmp_dir = output_dir + "/tmp/"
        
        left_bam_name = args[0]
        right_bam_name = args[1]
        
        start_time = datetime.now()
        prepare_output_dir()
        init_logger(logging_dir + "tophat.log")

        th_logp()
        th_log("Beginning TopHat run (v"+get_version()+")")
        th_logp("-----------------------------------------------")

        global run_log
        run_log = open(logging_dir + "run.log", "w", 0)
        global run_cmd
        run_cmd = " ".join(argv)
        print >> run_log, run_cmd

        check_samtools()
        
        (left_sam_filename, left_header) = convert_bam_to_sam(left_bam_name, tmp_dir + "left.sam")
        left_sam_sortedbyread_filename = sort_by_read(left_sam_filename, tmp_dir + "left.byread.sam")
        
        (right_sam_filename, right_header) = convert_bam_to_sam(right_bam_name, tmp_dir + "right.sam")
        right_sam_sortedbyread_filename = sort_by_read(right_sam_filename, tmp_dir + "right.byread.sam")
        
        left_org_sam_sortedbyread_filename = tmp_dir + "left_org.sam"
        right_org_sam_sortedbyread_filename = tmp_dir + "right_org.sam"
        
        separate_streams(left_sam_sortedbyread_filename, 
                         right_sam_sortedbyread_filename,
                         left_header,
                         right_header,
                         left_org_sam_sortedbyread_filename,
                         right_org_sam_sortedbyread_filename)
                         
        left_sep_bam_sortedby_read = convert_sam_to_bam(left_org_sam_sortedbyread_filename, tmp_dir + "right.sep.byread.bam")
        sort_by_coord(left_sep_bam_sortedby_read, output_dir + "/left")
        
        right_sep_bam_sortedby_read = convert_sam_to_bam(right_org_sam_sortedbyread_filename, tmp_dir + "right.sep.byread.bam")
        sort_by_coord(right_sep_bam_sortedby_read, output_dir + "/right")
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
