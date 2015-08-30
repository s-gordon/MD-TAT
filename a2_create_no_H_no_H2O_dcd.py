#!/usr/bin/env python
# AUTHOR:   Shane Gordon
# FILE:     a1_other_analyses.py
# ROLE:     TODO (some explanation)
# CREATED:  2015-06-16 21:46:32
# MODIFIED: 2015-06-22 08:29:07

import os
import sys
import logging
import argparse
import subprocess
import shutil
from glob import glob

# Argparse


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser(
            description="""
                        Analysis script part 2. Concatenates DCD trajectory
                        files under ./MainJob_dir into a single, reduced
                        trajectory file. Files are named sequentially from
                        no_water_<n>.dcd where <n> is the replicate number.
                        The bulk of the work is delegated to CATDCD, found
                        under ../Scripts/Tools/.""",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-v', '--verbose', action="store_true",
                    help="Increase verbosity")
parser.add_argument('-s', '--stride', type=int, default=1,
                    help="""
                         Factor to sub-sample the trajectories by. e.g. For -s
                         10, we would only save one of every 10 frames.
                         """)

result = parser.parse_args()

"""
If verbosity set, change logging to debug.
Else leave at info
"""
if result.verbose:
    logging.basicConfig(format='%(levelname)s:\t%(message)s',
                        level=logging.DEBUG)
else:
    logging.basicConfig(format='%(levelname)s:\t%(message)s',
                        level=logging.INFO)


def check_dir(dir):
    """
    Check whether directory dir exists.
    If true continue. Else exit.
    """
    if not os.path.isdir(dir):
        logging.error('Path %s not found', dir)
        logging.error('Aborting')
        sys.exit()


def check_file(file):
    """
    Check whether directory dir exists.
    If true continue. Else exit.
    """
    if not os.path.isfile(file):
        logging.error('Path %s not found', file)
        logging.error('Aborting')
        sys.exit()


def delete_dir(dir):
    """
    Check whether directory dir exists.
    If true delete and remake.
    """
    if os.path.exists(dir):
        shutil.rmtree(dir)
    os.makedirs(dir)


def delete_file(file):
    """
    Check whether file exists.
    If true delete.
    """
    if os.path.isfile(file):
        os.remove(file)


def parse_num_frames(prefix, i, catdcd):
        nf = subprocess.check_output('{catdcd} -num {p}_{iter}.dcd |\
                grep "Total frames:"|\
                awk \'{{print $3}}\''.format(
                    p=prefix, catdcd=catdcd, iter=i), shell=True)
        nf = nf.replace('\n', '')
        f = open('number_frames_{iter}.txt'.format(iter=i), 'w')
        f.write(nf)
        f.close()


def run_command(command):
        p = subprocess.Popen(command, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        return p.communicate()


def check_cmd(cmd):
    try:
        subprocess.check_call(['%s' % cmd], shell=True)
    except subprocess.CalledProcessError:
        pass  # handle errors in the called executable
    except OSError:
        logging.error('Command %s not found' % cmd)
        sys.exit()


def make_dir(dir):
    """
    Check whether dir exists.
    If true make it.
    """
    if not os.path.exists(dir):
        os.makedirs(dir)

dcdfile_list = glob.glob('dcdfile_list_*.txt')
catdcd = '../Scripts/Tools/catdcd'

dir_list = []
with open('../.dir_list.txt') as dir:
    for line in dir:
        dir_list.append(line)

for l in sorted(dcdfile_list):
        dir = l.rstrip('\n')
        i = subprocess.check_output(
                        "echo {0} | sed 's/.*_//' | sed 's/\.*//' | \
                                        sed 's/\.[^.]*$//'".format(dir),
                        shell=True)
        i = i.replace('\n', '')
        iter = 0
        with open(l) as f:
                for dcd in f:
                        dcd = dcd.replace('\n', '')
                        logging.debug('\tProcessing dcd: ' + dcd)
                        try:
                                run_command('{catdcd} -otype dcd -i no_water.text \
                                        -o {0}_temp_{1:04d}.dcd {dcd}'.format(
                                            i, iter, dcd=dcd, catdcd=catdcd))
                        except OSError as e:
                                logging.error(e)
                                logging.error("failed")
                                sys.exit()
                        iter = iter+1
        try:
                run_command('{0} -otype dcd -stride {2} -o no_water_{1}.dcd \
                        -dcd {1}_temp_????.dcd'.format(
                            catdcd, i, result.stride))
                if glob.glob('{0}_temp_*.dcd'.format(i)):
                        for idcd in glob.glob('{0}_temp_*.dcd'.format(i)):
                                delete_file(idcd)
                try:
                        parse_num_frames('no_water', i, catdcd)
                except OSError as e:
                        logging.error(e)
                        logging.error("Failed to work out the number of \
                                frames in the resulting trajectory")
                        sys.exit()
        except OSError as e:
                logging.error(e)
                logging.error("failed")
                sys.exit()
