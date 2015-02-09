#!/usr/bin/python

"""
Perform file maintenance operations on a directory
Dependencies:
  - python 2.7+ (tested with 2.7.3)
"""

__version__ = "1.0"
__author__ = "Albert Vilella (avilella@aptamex.com)"
__license__ = "GPL"

import argparse
import sys
import os
import re

samtools_exe         = "/opt/ontoitsoftware/aptamex/thirdparty/samtools-1.1/samtools"
sge_bashrc           = "/opt/ontoitsoftware/aptamex/thirdparty/sge/bash/dot_bashrc"
extract_variants_exe = "/opt/ontoitsoftware/aptamex/thirdparty/htsmanager/htsmanager_deps/extract_variants"
bgzip_exe            = "/opt/ontoitsoftware/aptamex/thirdparty/htsmanager/htsmanager_deps/bgzip"
bedtools_exe         = "/opt/ontoitsoftware/aptamex/thirdparty/htsmanager/htsmanager_deps/bedtools/bin/bedtools"
hg19_genome          = "/opt/ontoitsoftware/aptamex/thirdparty/htsmanager/htsmanager_deps/hg19.genome"

data_scratch         = "/data/aptamex/scratch/tmp/users"

args = None
list_operations = ['remove', 'move', 'rehead-move', 'rehead-delete', 'bam2cram', 'cram2bam', 'bamlet']
slop = 500

def main():
    """parse command line input and call appropriate functions"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, \
             description=("Perform file maintenance on a directory:\n"
            "This tool intends to complement directory-based file management capabilities\n"
                         "with finding files inside a directory and performing operations on them.\n"), epilog="""
 Examples:
     htsmanager.py -d /data/aptamex/project_id/analysis_id --name "*.bam" --minsize 10G --operation rehead-delete

         Will look for files with "*.bam" pattern in directory, and for those bigger than the minimum size, it will apply the 'rehead-delete' operation

     htsmanager.py -d /data/aptamex/project_id/analysis_id --maxdepth 1 --name "*.bam" --minsize 1M --operation bam2cram --dryrun

         Will look for files with "*.bam" pattern in depth 1 of directory, and for those bigger than 1M, in dryrun mode will print out the operations that would apply operation

    htsmanager.py -d /data/aptamex/project_id/analysis_id --operation move --qsubit --name "*.bam" 2>/dev/null 1>move.sh

         Will look for file with "*.bam" pattern, and print out qsubit command lines for the move operation. You can then do "bash move.sh" to submit each one to the grid.

    htsmanager.py -d /data/aptamex/project_id/analysis_id --maxdepth 5 --operation remove --not_name "*.bam" --not_name "*.vcf.gz" --name "*" --minsize 10M

         Will aggressively delete everything other than bam and vcf.gz files larger than 10M

 Operations:
     remove:         Delete selected files

     rehead-delete:  Create a "*.bam.header.bam" with only the bam header information and no reads, delete original bam and symlink.

                     Using 'rehead', the meta information of the bam file is kept but the reads are gone.
                     You can still run metrics with the empty bam and associated vcf, but the bam metrics will be void.

     rehead-move:    Create a "*.bam.header.bam" with only the bam header information and no reads, then move the selected bam files
                     to scratch tmp and symlinks them to the original place.

     move:           Moves the selected files to scratch tmp and symlinks them to the original place.

     bam2cram:       Convert from bam format to cram format. Average file reduction is around 50%.

     cram2bam:       Restore "*.bam.cram" files to bam format.

     bamlet:         Find "*.bam" files and respective "*.genome.vcf.gz" files and create a bamlet with only the reads
                     contained inside windows around variants, delete original bam and symlink.

 """)

    parser.add_argument('--version', '-V', action='version', version="%(prog)s " + __version__ + "\n")
    parser.add_argument("--debug",     help="debug", action="store_true")
    parser.add_argument("-v", "--verbose", help="verbose", action="store_true")
    parser.add_argument("--dryrun", "--simulate", help="simulate operations without applying them", action="store_true")
    parser.add_argument("--qsubit", "--qsubit", help="print out qsubit commands", action="store_true")
    parser.add_argument("-d", "--dir", "--directory", help="input directory where to search for files", required=True)
    parser.add_argument("--movedir", "--move_directory", help="directory where to move the files if operation move is used", required=False)
    parser.add_argument("--maxdepth",  help="find maxdepth (default: %(default)s)", default=3)
    parser.add_argument("--name", "--file_name",  help="UNIX find -name option", action='append', required=True)
    parser.add_argument("--type", "--find_type",  help="UNIX find -type option (default: %(default)s)", default="f", required=False)
    parser.add_argument("--notname", "--not_name",  help="UNIX find -not -name option", action='append')
    parser.add_argument("--samtools", "--samtools_exe",  help="samtools binary (default: %(default)s)", default=samtools_exe)
    parser.add_argument("--minsize", "--minimum_file_size", help="minimum file size to perform operation (default: %(default)s)", default="1G")
    parser.add_argument("--operation",  help="operation to perform on files (default: %(default)s)", default="remove", choices=list_operations)
    parser.add_argument("--slop", "--bamlet_extend",  help="extend upstream and downstream window around coordinate (default: %(default)s)", default=slop)
    global args
    args = parser.parse_args()
    global qsub_cmd
    global qsub_mem
    qsub_cmd = list()
    qsub_mem = '2G'
    if args.qsubit:
        args.dryrun = True

    if args.debug:
        print "# DEBUG=1"

    # Start program
    if args.verbose:
        print "dir = %s" % args.dir

    cmd = ''

    # Defining names
    list_name = []
    this_name = ''
    for elem in args.name:
        this_str = " -name \"%s\" -type %s " % (elem, args.type)
        list_name.append(this_str)
    this_name = " -or ".join(list_name)

    # Defining not names
    list_notname = []
    this_notname = ''
    if args.notname:
        for elem in args.notname:
            this_str = " -not -name \"%s\" " % elem
            list_notname.append(this_str)
        this_notname = " ".join(list_notname)

    cmd = "find %s -maxdepth %s %s %s | sort" % (args.dir, args.maxdepth, this_name, this_notname)
    if args.verbose: print >> sys.stderr, "# %s" % cmd
    files = os.popen(cmd).read()

    # Loop over each file and perform operation
    for this_file in files.split("\n"):
        if (0 == len(this_file)): continue
        if not (is_min_size(this_file)):
            print >> sys.stderr, "# file = %s" % this_file
            print >> sys.stderr, "File too small to perform operation"
            continue
        print >> sys.stderr, "# file = %s" % this_file
        print >> sys.stderr, "# operation = %s" % args.operation
        try_operation(this_file, args.operation)
    return

########################################


def try_operation(this_file, this_operation):
    global qsub_mem
    global qsub_cmd
    if args.dryrun:
        perform_operation(this_file, args.operation, True)
    else:
        if query_yes_no("Perform operation to file?", default="yes"):
            print "Performing operation %s" % args.operation
            perform_operation(this_file, args.operation, False)
        else:
            print "Skipping file %s" % this_file

    if args.qsubit:
        cmd = " && ".join(qsub_cmd)
        print "source %s" % (sge_bashrc)
        print "qsubit %s \"%s\"" % (qsub_mem, cmd)
        qsub_cmd = list()
    return


def perform_operation(this_file, this_operation, simulate):
    if 'remove' == this_operation:
        remove(this_file, simulate)
    elif 'move' == this_operation:
        move(this_file, simulate)
    elif 'rehead-move' == this_operation:
        rehead_move(this_file, simulate)
    elif 'rehead-delete' == this_operation:
        rehead_delete(this_file, simulate)
    elif 'bam2cram' == this_operation:
        bam2cram(this_file, simulate)
    elif 'cram2bam' == this_operation:
        cram2bam(this_file, simulate)
    elif 'bamlet' == this_operation:
        bamlet(this_file, simulate)
    else:
        print "Operation not defined in list:"

    return


def run_cmd(cmd, simulate):
    global qsub_cmd
    if (args.verbose and not simulate):
        print >> sys.stderr, "# %s" % cmd
    if args.qsubit:
        qsub_cmd.append(cmd)
        return
    if simulate:
        print "%s" % cmd
    else:
        ret = os.popen(cmd).read()
        if args.verbose:
            print >> sys.stderr, "# %s" % ret
    return


def remove(this_file, simulate):
    cmd = ''
    recursive = ''
    if "d" == args.type:
        recursive = "-r"
    cmd = "rm -f %s %s" % (recursive, this_file)
    run_cmd(cmd, simulate)
    return


def move(this_file, simulate):
    cmd = ''
    global qsub_mem
    qsub_mem = '1G'
    destroot = data_scratch
    if args.movedir:
        destroot = args.movedir
    userhome = os.path.expanduser('~')
    username = os.path.split(userhome)[-1]
    this_dir, this_filename = os.path.split(this_file)
    destdir = os.path.join(destroot, username) + this_dir
    if args.verbose:
        print "# Destdir " + destdir
    if not os.path.exists(destdir):
        os.makedirs(destdir)
    cmd = "mv %s %s" % (this_file, destdir)
    run_cmd(cmd, simulate)
    new_file = os.path.join(destdir, this_filename)
    cmd = "ln -s %s %s" % (new_file, this_file)
    run_cmd(cmd, simulate)
    return


def rehead_move(this_file, simulate):
    cmd = ''
    global qsub_mem
    qsub_mem = '1G'
    destroot = data_scratch
    if args.movedir:
        destroot = args.movedir
    userhome = os.path.expanduser('~')
    username = os.path.split(userhome)[-1]
    this_dir, this_filename = os.path.split(this_file)
    destdir = os.path.join(destroot, username) + this_dir
    if args.verbose:
        print "# Destdir " + destdir
    if not os.path.exists(destdir):
        os.makedirs(destdir)
    cmd = "%s view -Hb %s > %s.header.bam" % (args.samtools_exe, this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "mv %s %s" % (this_file, destdir)
    run_cmd(cmd, simulate)
    new_file = os.path.join(destdir, this_filename)
    cmd = "ln -s %s %s" % (new_file, this_file)
    run_cmd(cmd, simulate)
    return


def rehead_delete(this_file, simulate):
    cmd = ''
    cmd = "%s view -Hb %s > %s.header.bam" % (args.samtools_exe, this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "mv %s %s.full.bam" % (this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "ln -s %s.header.bam %s" % (this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "rm %s.full.bam" % (this_file)
    run_cmd(cmd, simulate)
    return


def bam2cram(this_file, simulate):
    cmd = ''
    cmd = "%s view -Hb %s > %s.header.bam" % (args.samtools_exe, this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "%s view -h -C %s > %s.cram" % (args.samtools_exe, this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "rm %s" % (this_file)
    run_cmd(cmd, simulate)
    cmd = "ln -s %s.header.bam %s" % (this_file, this_file)
    run_cmd(cmd, simulate)
    return


def cram2bam(this_file, simulate):
    cmd = ''
    this_file_bam = re.sub('\.cram$', '', this_file)
    cmd = "%s view -h -b %s > %s" % (args.samtools_exe, this_file, this_file_bam)
    run_cmd(cmd, simulate)
    cmd = "rm %s" % (this_file)
    run_cmd(cmd, simulate)
    return

def bamlet(this_file, simulate):
    cmd = ''
    this_file_gvcf = re.sub('\.bam$', '.genome.vcf.gz', this_file)
    this_file_bed_gz  = re.sub('\.bam$', '.slop.bed.gz', this_file)
    cmd = "gunzip -c %s | %s | %s merge -i stdin | %s slop -b %s -g %s -i stdin | %s merge | gzip -c > %s" % (this_file_gvcf, extract_variants_exe, bedtools_exe, bedtools_exe, args.slop, hg19_genome, bedtools_exe, this_file_bed_gz)
    run_cmd(cmd, simulate)
    cmd = "%s view -h -b -L %s %s > %s.bamlet.bam" % (args.samtools_exe, this_file_bed_gz, this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "%s index %s.bamlet.bam" % (args.samtools_exe, this_file)
    run_cmd(cmd, simulate)
    cmd = "mv %s %s.full.bam" % (this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "ln -s %s.bamlet.bam %s" % (this_file, this_file)
    run_cmd(cmd, simulate)
    cmd = "rm %s.full.bam" % (this_file)
    run_cmd(cmd, simulate)
    return

def is_min_size(this_file):
    this_file_size = os.path.getsize(this_file)
    if args.verbose:
        print("[%s] size = %s (minsize %s)" % (this_file, bytes2human(this_file_size), args.minsize))
    if (this_file_size < human2bytes(args.minsize)):
        ret = False
    else:
        ret = True
    return ret


def bytes2human(n, format='%(value).1f %(symbol)s', symbols='customary'):
    """
    Convert n bytes into a human readable string based on format.
    symbols can be either "customary", "customary_ext", "iec" or "iec_ext",
    see: http://goo.gl/kTQMs

      >>> bytes2human(0)
      '0.0 B'
      >>> bytes2human(0.9)
      '0.0 B'
      >>> bytes2human(1)
      '1.0 B'
      >>> bytes2human(1.9)
      '1.0 B'
      >>> bytes2human(1024)
      '1.0 K'
      >>> bytes2human(1048576)
      '1.0 M'
      >>> bytes2human(1099511627776127398123789121)
      '909.5 Y'

      >>> bytes2human(9856, symbols="customary")
      '9.6 K'
      >>> bytes2human(9856, symbols="customary_ext")
      '9.6 kilo'
      >>> bytes2human(9856, symbols="iec")
      '9.6 Ki'
      >>> bytes2human(9856, symbols="iec_ext")
      '9.6 kibi'

      >>> bytes2human(10000, "%(value).1f %(symbol)s/sec")
      '9.8 K/sec'

      >>> # precision can be adjusted by playing with %f operator
      >>> bytes2human(10000, format="%(value).5f %(symbol)s")
      '9.76562 K'
    """

    SYMBOLS = {
        'customary'     : ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'),
        'customary_ext' : ('byte', 'kilo', 'mega', 'giga', 'tera', 'peta', 'exa',
                           'zetta', 'iotta'),
        'iec'           : ('Bi', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi'),
        'iec_ext'       : ('byte', 'kibi', 'mebi', 'gibi', 'tebi', 'pebi', 'exbi',
                           'zebi', 'yobi'),
    }

    n = int(n)
    if n < 0:
        raise ValueError("n < 0")
    symbols = SYMBOLS[symbols]
    prefix = {}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i+1)*10
    for symbol in reversed(symbols[1:]):
        if n >= prefix[symbol]:
            value = float(n) / prefix[symbol]
            return format % locals()
    return format % dict(symbol=symbols[0], value=n)

def human2bytes(s):
    """
    Attempts to guess the string format based on default symbols
    set and return the corresponding bytes as an integer.
    When unable to recognize the format ValueError is raised.

      >>> human2bytes('0 B')
      0
      >>> human2bytes('1 K')
      1024
      >>> human2bytes('1 M')
      1048576
      >>> human2bytes('1 Gi')
      1073741824
      >>> human2bytes('1 tera')
      1099511627776

      >>> human2bytes('0.5kilo')
      512
      >>> human2bytes('0.1  byte')
      0
      >>> human2bytes('1 k')  # k is an alias for K
      1024
      >>> human2bytes('12 foo')
      Traceback (most recent call last):
          ...
      ValueError: can't interpret '12 foo'
    """


    SYMBOLS = {
        'customary'     : ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'),
        'customary_ext' : ('byte', 'kilo', 'mega', 'giga', 'tera', 'peta', 'exa',
                           'zetta', 'iotta'),
        'iec'           : ('Bi', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi'),
        'iec_ext'       : ('byte', 'kibi', 'mebi', 'gibi', 'tebi', 'pebi', 'exbi',
                           'zebi', 'yobi'),
    }

    init = s
    num = ""
    while s and s[0:1].isdigit() or s[0:1] == '.':
        num += s[0]
        s = s[1:]
    num = float(num)
    letter = s.strip()
    for name, sset in SYMBOLS.items():
        if letter in sset:
            break
    else:
        if letter == 'k':
            # treat 'k' as an alias for 'K' as per: http://goo.gl/kTQMs
            sset = SYMBOLS['customary']
            letter = letter.upper()
        else:
            raise ValueError("can't interpret %r" % init)
    prefix = {sset[0]:1}
    for i, s in enumerate(sset[1:]):
        prefix[s] = 1 << (i+1)*10
    return int(num * prefix[letter])


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

########################################
if __name__ == "__main__":
    main()

1
