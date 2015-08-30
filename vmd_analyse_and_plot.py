#!/usr/bin/env python
# AUTHOR:       Shane Gordon
# FILE:         a5_other_analyses.py
# ROLE:         TODO (some explanation)
# CREATED:      2015-06-16 21:46:32

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import argparse
import subprocess
import shutil
from glob import glob

# Pyplot formatting bits


def init_plotting():

    plt.rcParams['figure.figsize'] = (8, 3)
    plt.rcParams['figure.autolayout'] = True
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['lines.linewidth'] = 0.5
    plt.rcParams['legend.loc'] = 'center right'
    plt.rcParams['legend.fontsize'] = 6
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')


def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(
        description="""Batch analysis script. Takes output from a2 and allows
for interpretation of RMSD, RMSF, SASA, secondary structure and more.
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

    parser.add_argument('--rmsd',  action="store_true", default=False, help="""
Calculate and plot root-mean squared deviation of selection backbone over the
time course of the simulation. Useful for measuring structure equilibration.
                        """)

    parser.add_argument('-s', '--selection', default='protein', help="""
Protein selection to use in analyses. Must be a valid selection. At present,
there are no explicit checks for making sure what you pass is valid. Very
fragile!
                        """)

    parser.add_argument('--align_sel', default='protein and name CA',
                        help="""
Optional selection to use for protein alignment, if
applicable. Currently there are no sanity checks for this,
so assume that it is very fragile.Syntax documentation
for VMD selections can be found here:
ks.uiuc.edu/Research/vmd/vmd-1.2/ug/vmdug_node137.html
                        """)

    parser.add_argument('--rmsf', action="store_true", default=False, help="""
RMSF
                        """)

    parser.add_argument('--sasa',  action="store_true", default=False, help="""
Calculate and plot solvent-accessible surface area (SASA)
of selection over the time course of the trajectory. Useful
for measuring changes in protein interfaces.
Computationally intensive.
                        """)

    parser.add_argument('--rsasa',  action="store_true", default=False,
                        help="""
Residue-level SASA.
                        """)

    parser.add_argument('--rg',  action="store_true", default=False,
                        help="""
Calculate and plot radius of gyration of selection over the
time course of the trajectory. Useful for measuring
structural changes
                        """)

    parser.add_argument('--dccm',  action="store_true", default=False,
                        help="""
Calculate and plot dynamic cross-correlation matrices
(DCCMs) for each trajectory. Useful for measuring
interaction networks within a structure (both positive and
negative).
                        """)

    parser.add_argument('--phi',  action="store_true", default=False,
                        help="""
Calculate and plot residue phi angles for each trajectory.
useful for measuring structural perturbations at the
residue level.
                        """)

    parser.add_argument('--psi',  action="store_true", default=False,
                        help="""
Calculate and plot residue psi angles for each trajectory.
useful for measuring structural perturbations at the
residue level.
                        """)

    # Argument strings are accessed through the argparser 'args'
    # result = parser.parse_args()
    args = vars(parser.parse_args())

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

    # Directory where raw data is placed, if generated
    data_dir = args['raw_data']

    # File containing a list of simulation replicates. We use this to work out
    # how many times other subroutines need to run. This could be better.
    dir_list = "../.dir_list.txt"

    # Path to directory containing VMD tcl scripts we need to use later
    script_dir = "../Scripts/Analysis_Scripts"

    # shell command to call vmd in CLI mode
    vmd_cmd = 'vmd -dispdev text'

    # }}}

    # Run VMD analyses
    res = {
        args['rmsd']: "{0}/analysis_rmsd.tcl".format(script_dir),
        args['rmsf']: "{0}/analysis_rmsf.tcl".format(script_dir),
        args['sasa']: "{0}/analysis_sasa.tcl".format(script_dir),
        args['rsasa']: "{0}/analysis_rsasa.tcl".format(script_dir),
        args['rg']: "{0}/analysis_rg.tcl".format(script_dir),
        args['dccm']: "{0}/analysis_dccm.tcl".format(script_dir),
        args['phi']: "{0}/analysis_phi.tcl".format(script_dir),
        args['psi']: "{0}/analysis_psi.tcl".format(script_dir)
    }

    sub_dirs = get_dirs(dir_list)

    # Checks
    for f in [dir_list]:
        check_file(f)
    make_dir(data_dir)
    for r in sub_dirs:
        dname = '{parent}/sim_{rep}'.format(parent=data_dir, rep=r)
        make_dir(dname)

    for r, a in res.iteritems():
        if r is True:
            for replicate in sub_dirs:
                command_catch_error(
                    '{vmd} -e {script} -args no_water_{rep}.dcd {data_dir} \
                    {align_sel}'.format(vmd=vmd_cmd, script=a,
                                        rep=replicate, data_dir=data_dir,
                                        align_sel=args['align_sel']))

    # If plot is set to true, attempt to plot output data using matplotlib.
    if args['plot'] is False:
        logging.debug("""
Optional plotting defaulting to False. Plot output
using the '--plot' flag.
                      """)
    elif args['plot'] is True:
        logging.debug("""
Optional plotting set to True.
                      """)

        # Directory where processed (plotted) data is placed, if generated by
        # matplotlib
        plot_dir = args['plot_dir']

        make_dir(plot_dir)
        for r in sub_dirs:
            dname = '{parent}/sim_{rep}'.format(parent=plot_dir, rep=r)
            make_dir(dname)

        # Load matplotlib params
        init_plotting()

        # Iterate over each replicate simulation and plot data is found
        for replicate in sub_dirs:
            out_d = "{0}/sim_{1}".format(plot_dir, replicate)

            try:
                rg_dict = {
                    'fname': '\
                    {0}/sim_{1}/protein_radius_gyration_{1}.txt'.format(
                        data_dir, replicate),
                    'xlabel': 'Simulation time (ns)',
                    'ylabel': 'R$_g$ ($\AA$)',
                    'ymin': 0,
                    'ofile': '{0}/rg_plot_{1}'.format(out_d, replicate)
                }
                rmsd_dict = {
                    'fname':
                    '{0}/sim_{1}/rmsd_protein_{1}.txt'.format(
                        data_dir, replicate),
                    'xlabel': 'Simulation time (ns)',
                    'ylabel': 'RMSD ($\AA$)',
                    'ymin': 0,
                    'ofile': '{0}/rmsd_plot_{1}'.format(out_d, replicate)
                }
                rmsf_dict = {
                    'fname': '{0}/sim_{1}/rmsf_protein_backbone'.format(
                        data_dir, replicate),
                    'all_avg':
                    '{0}/sim_{1}/rmsf_all_protein_backbone'.format(
                        data_dir, replicate),
                    'xlabel': 'Residue No.',
                    'ylabel': 'RMSF ($\AA$)',
                    'ymin': 0,
                    'ymax': 10,
                    'ofile': '{0}/rmsf_plot_{1}'.format(out_d, replicate)
                }
                sasa_dict = {
                    'fname': '{0}/sim_{1}/protein_sasa_{1}.txt'.format(
                        data_dir, replicate),
                    'xlabel': 'Simulation time (ns)',
                    'ylabel': 'Solvent-accessible surface area ($\AA^2$)',
                    'ymin': 0,
                    'ofile': '{0}/sasa_plot_{1}'.format(out_d, replicate)
                }
                rsasa_dict = {
                    'fname': '{0}/sim_{1}/protein_sasa_{1}.txt'.format(
                        data_dir, replicate),
                    'xlabel': 'Simulation time (ns)',
                    'ylabel': 'Solvent-accessible surface area ($\AA^2$)',
                    'ymin': 0,
                    'ofile': '{0}/sasa_plot_{1}'.format(out_d, replicate)
                }

                # If true, plot rmsd, sasa, and radius of gyrations data
                for p in [rmsd_dict, sasa_dict, rg_dict, rsasa_dict]:
                    basic_plot(p)

                # For each replicate, plot fractional rmsf values and total
                # rmsf values
                rmsf_shared_axis(rmsf_dict, out_d)

            except OSError as e:
                logging.error(e)


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
                logging.debug("Directory {0} found.\nRemoving {0}".format(dir))
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
        logging.debug('Directory %s already exists!' % dir)


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


def get_dirs(dirlist):
    sub_dirs = []
    with open(dirlist) as f:
        for line in f:
            line = line.rstrip('\n')
            i = subprocess.check_output(
                "echo {0} | sed 's/.*_//' | sed 's/\.*//'".format(line),
                shell=True)
            i = i.replace('\n', '')
            sub_dirs.append(i)
            logging.debug('Subdirectory indices are %s' % i)
    return sub_dirs


def basic_plot(proc):
    if os.path.isfile(proc['fname']):
        logging.debug("Plotfile %s found. Attempting to plot using\
                      matplotlib" % proc['fname'])
        data = np.loadtxt(proc['fname'])
        plt.subplot(111)
        plt.xlabel(proc['xlabel'])
        plt.ylabel(proc['ylabel'])
        plt.plot(data[:, 0]/10, data[:, 1, ], lw=2)
        plt.ylim(proc['ymin'])
        plt.savefig('{ofile}.pdf'.format(ofile=proc['ofile']))
        plt.close()


def rmsf_shared_axis(d, out):
    if glob('{fn}*'.format(fn=d['fname'])):
        a = [i for i in sorted(glob('{fn}*'.format(fn=d['fname'])))]
        n = len(a)
        count = 1
        color = iter(plt.cm.Blues(np.linspace(0, 1, n)))
        for f in a:
            if os.path.isfile(f):
                c = next(color)
                o = os.path.splitext(f)[0]
                data = np.loadtxt('{o}.txt'.format(o=o))
                path, prefix = os.path.split('{o}.txt'.format(o=o))
                ax = plt.subplot(111)
                plt.plot(data[:, 0], data[:, 1, ], c=c,
                         label='Fraction {count} of 5'.format(count=count))
                count = count+1
    if glob('{fn}*'.format(fn=d['all_avg'])):
        if os.path.isfile(f):
            o = os.path.splitext(f)[0]
            data = np.loadtxt('{o}.txt'.format(o=o))
            path, prefix = os.path.split('{oprefix}.txt'.format(oprefix=o))
            ax = plt.subplot(111)
            ax.legend(loc='upper left')
    plt.xlabel('{label}'.format(label=d['xlabel']))
    plt.ylabel('{label}'.format(label=d['ylabel']))
    plt.ylim(d['ymin'], d['ymax'])
    plt.savefig('{oprefix}/rmsf.pdf'.format(oprefix=out))
    plt.close()


if __name__ == '__main__':
    main()
