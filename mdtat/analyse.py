#!/usr/bin/env python
# AUTHOR:       Shane Gordon
# CREATED:      2015-06-16 21:46:32

# Note: This package is fairly poorly organised. I should do something about
# it, but I'm putting it in the TODO list for now

# Set matplotlib to have non-interactive backend
try:
    import matplotlib as mpl
except ImportError:
    pass
else:
    mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import argparse
import subprocess
import shutil
from glob import glob
from mdtat.analysis.rmsd import compute_rmsd
from mdtat.analysis.rg import compute_rg
from mdtat.analysis.sasa import compute_sasa
import time
plt.style.use('ggplot')


class DataDirs:
    """
    Class wrapper for assigning input/output directory structure.
    """
    def __init__(self, tcl, data, plots):
        self.tcl = tcl
        self.data = data
        self.plots = plots


def benchmark(function, *args):
    start = time.time()
    function(*args)
    end = time.time()
    t = end - start
    return t


def init_plotting():

    plt.rcParams['figure.figsize'] = (8, 3)
    plt.rcParams['figure.autolayout'] = True


def check_dict(dict):
    """
    If any dictionary items have been set, return True. Else False.
    """
    for i, j in dict.iteritems():
        if i is True:
            return True


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
        logging.warn("Directory {0} found.\nRemoving {0}".format(dir))
        os.makedirs(dir)


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
    If directory does not exist, make it
    """
    if not os.path.exists(dir):
        logging.debug('Directory %s not found. Creating.' % dir)
        os.makedirs(dir)
    else:
        logging.warn('Directory %s already exists!' % dir)


def command_catch_error(command):
    """
    Wrapper for shell commands. Stdin/stdout returned.
    """
    try:
        p = subprocess.Popen(command, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        p.wait()
        return p.communicate()
    except OSError as e:
        """
        If this fails, report the error.
        """
        logging.error(e)
        logging.error("Command {com} failed. Please troubleshoot this and try \
                      again".format(com=command))
        sys.exit()


def autocorr(x):
    "Compute an autocorrelation with numpy"
    x = x - np.mean(x)
    result = np.correlate(x, x, mode='full')
    result = result[result.size//2:]
    return result / result[0]


def basic_plot(data, xlabel, ylabel, output):
    logging.debug('Attempting to plot to %s' % output)
    f, ax = plt.subplots(2)
    count = np.arange(0, len(data))
    for row, i in zip(data, count):
        ax[0].set_xlabel(xlabel)
        ax[0].set_ylabel(ylabel)
        ax[0].plot(row, label=i)
        # ax[0].set_ylim(0,)
        ax[0].legend()
        ax[1].set_xlabel(xlabel)
        ax[1].set_ylabel('ACF')
        ax[1].semilogx(autocorr(row), label=i)
        ax[1].legend()
    plt.savefig(output, format='pdf')
    plt.close()


def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(
        description="""
Batch analysis script. Calculates several observables including RMSD, and
radius of gyration.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v', '--verbose',  action="store_true", help="""
Increase verbosity.
                        """)

    parser.add_argument('-d', '--raw-data',  action="store_true",
                        default='data', help="""
Directory to output raw data, such as RMSF, RMSD, etc.
                        """)

    parser.add_argument('-p', '--plot-dir',  action="store_true",
                        default='plots', help="""
Directory to output plots (if generated), such as RMSF, RMSD, etc.
                        """)

    parser.add_argument('--plot',  action="store_true", help="""
Optionally, plot output data using matplotlib into output folder 'plot'.
                        """)

    parser.add_argument('--sel', type=str, nargs='+', help="""
Selection text.
                        """)

    parser.add_argument('--rmsd',  action="store_true", default=False, help="""
Calculate and plot root-mean squared deviation of selection backbone over the
time course of the simulation. Useful for measuring structure equilibration.
                        """)

    parser.add_argument('--rg',  action="store_true", default=False,
                        help="""
Calculate and plot radius of gyration of selection over the
time course of the trajectory. Useful for measuring
structural changes
                        """)

    parser.add_argument('--trajfiles',  nargs='+', metavar="FILE",
                        help="""
Trajfiles.
                        """)

    parser.add_argument('--step', type=int, default=1,
                        help="""
Optional step factor for analysis functions.
                        """)

    # Argument strings are accessed through the argparser 'args'
    args = vars(parser.parse_args())

    # Return argparse help if no arguments are specified. Catches situations
    # where only verbosity is set and does the same thing.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2:
        if args['verbose'] is True:
            parser.print_help()
            sys.exit(1)

    # Logger-control of feed-out
    # Normal output is prescribed to info
    # Debug output (accessed through -v/--verbose) is prescribed to debug
    """
    If verbosity set, change logging to debug.
    Else leave at info
    """
    if args['verbose'] is True:
        logging.basicConfig(format='%(levelname)s:\t%(message)s',
                            level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s:\t%(message)s',
                            level=logging.INFO)

    # Define data output directories
    # Access using out_dir.data, out_dir.plots
    out_dir = DataDirs('tcl_scipts', args['raw_data'], args['plot_dir'])

    # File containing a list of simulation replicates. We use this to work out
    # how many times other subroutines need to run. This could be better.
    # dir_list = "../.dir_list.txt"

    if args['trajfiles'] is not None:
        tfile = args['trajfiles']
        # Make sure that the filenames match existing files
        for f in tfile:
            check_file(f)
    else:
        tfile = sorted(glob('traj*.dcd'))

    rmsd, rg, sasa = ([] for i in range(3))

    top = 'reduced.pdb'

    if any([
            args['rmsd'],
            args['rg'],
            args['sasa']
    ]) is True:

        sel = " ".join(args['sel'])
        step = args['step']
        make_dir(out_dir.data)
        o = out_dir.data

        if args['rmsd'] is True:
            rmsd = []
            outfile = '{dir}/rmsd.txt'.format(dir=o)
            for traj in tfile:
                start_time = time.time()
                rmsd.append(compute_rmsd(traj, top, sel, step))
                benchmark_time = time.time() - start_time
                logging.debug('Benchmark time for {} was {:.2f} seconds'
                              .format(traj, benchmark_time))
            with open(outfile, 'w') as f:
                # This seems to break down when the array has only a single
                # dimension
                f.write("\n".join(" ".join(map(str, x)) for x in (rmsd)))

        if args['rg'] is True:
            rg = []
            outfile = '{dir}/rg.txt'.format(dir=o)
            for traj in tfile:
                data = compute_rg(traj, top, step=step)
                rg.append(data)
            with open(outfile, 'w') as f:
                # This seems to break down when the array has only a single
                # dimension
                f.write("\n".join(" ".join(map(str, x)) for x in (rg)))

        if args['sasa'] is True:
            sasa = []
            outfile = '{dir}/sasa.txt'.format(dir=o)
            for traj in tfile:
                data = compute_sasa(traj, top)
                sasa.append(data)
            with open(outfile, 'w') as f:
                # This seems to break down when the array has only a single
                # dimension
                f.write("\n".join(" ".join(map(str, x)) for x in (sasa)))

    plot = {
        'rg': {
            'arg': 'rg',
            'array': rg,
            'xlabel': 'Frame No.',
            'ylabel': 'R$_g$ (nm)',
            'ymin': 0,
            'ofile': 'rg.pdf'
        },
        'rmsd': {
            'arg': 'rmsd',
            'array': rmsd,
            'xlabel': 'Frame No.',
            'ylabel': 'RMSD (nm)',
            'ymin': 0,
            'ofile': 'rmsd.pdf'
        },
        'sasa': {
            'arg': 'sasa',
            'array': sasa,
            'xlabel': 'Frame No.',
            'ylabel': 'SASA (nm$^3$)',
            'ymin': 0,
            'ofile': 'sasa.pdf'
        }

    }

    if args['plot'] is False:
        logging.debug("Optional plotting defaulting to False." +
                      "Plot output using the '--plot' flag.")
    elif args['plot'] is True:
        logging.debug("Optional plotting set to True.")
        make_dir(out_dir.plots)
        d = out_dir.plots
        for key, value in plot.iteritems():
            if args[value['arg']] is True:
                data = value['array']
                xlabel = value['xlabel']
                ylabel = value['ylabel']
                f = '{root}/{f}'.format(root=d, f=value['ofile'])
                basic_plot(data, xlabel, ylabel, f)


if __name__ == '__main__':
    main()
