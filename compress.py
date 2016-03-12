#!/usr/bin/env python
# AUTHOR:   Shane Gordon
# CREATED:  Sat 31 Oct 2015 15:58:17 AEDT

import mdtraj as md
import logging
import argparse
import sys
import numpy as np
import os
import time
from glob import glob


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser(description="""
TODO
                  """,
                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--verbose', action="store_true",
                    help="Increase verbosity")
parser.add_argument('-i', '--input', nargs='+', required=True,
                    help="Files. TODO")
parser.add_argument('--stride', type=int, default=1,
                    help="Optional striding of trajectories.")
parser.add_argument('--topology', type=str, required=True,
                    help="Topology file.")
parser.add_argument('-sel', '--selection', type=str, required=True,
                    help="""
                    Reduced selection text. Must conform to MDTraj rules.
                    """)

# Argument strings are accessed through the argparser 'args'
args = vars(parser.parse_args())

"""
If verbosity set, change logger to debug.
Else leave at info
"""
if args['verbose']:
    logging.basicConfig(format='%(levelname)s:\t%(message)s',
                        level=logging.DEBUG)
else:
    logging.basicConfig(format='%(levelname)s:\t%(message)s',
                        level=logging.INFO)

logger = logging.getLogger(__name__)


def check_dir(dir):
    """
    Check whether directory dir exists.
    If true continue. Else exit.
    """
    if not os.path.isdir(dir):
        logger.warn('Path %s not found', dir)
        logger.error('Aborting')
        sys.exit(2)
    elif os.path.isdir(dir):
        logger.debug('Path %s exists and is a directory', dir)


def delete_file(file):
    """
    Check whether file exists.
    If true delete.
    """
    if os.path.isfile(file):
        logging.warning('File {} found.'.format(file))
        logging.warning('Removing {}.'.format(file))
        os.remove(file)


def get_indices(topology, sel):
    t = md.load(topology)
    try:
        indices = t.topology.select(sel)
        return indices
    except ValueError:
        logging.error('Atom selection invalid.\nTry something more sensible.')
        sys.exit(2)


def reduced_topology(topology, indices, outfile=None):
    t = md.load(topology, atom_indices=indices)
    if outfile is not None:
        t.save_pdb(outfile)
    return t


def main():

    dirs = args['input']
    topology = args['topology']
    sel = args['selection']
    stride = args['stride']
    indices = get_indices(topology, sel)
    top = reduced_topology(topology, indices, 'reduced.pdb')

    for dir, traj in zip(dirs, np.arange(0, len(dirs))):
        combined = 'traj{:0>4d}.dcd'.format(traj)
        delete_file(combined)
        logging.info('Concatenating trajectory files under {} into {}.'
                     .format(dir, combined))

        check_dir(dir)
        traj_list = sorted(glob('{}/OutputFiles/[0-9]*.dcd'.format(dir)))

        # if len(traj_list) is not 0:
        #     logging.info('Concatenating trajectory files under {} into {}.'
        #                  .format(dir, combined))
        #     for traj in traj_list:
        #         start_time = time.time()
        #         try:
        #             t = md.load(traj, top=top, stride=stride,
        #                         atom_indices=indices)
        #             if not os.path.isfile(combined):
        #                 try:
        #                     t.save_dcd(combined)
        #                     logger.debug('Combined output file {} not found.'
        #                                  .format(combined))
        #                     logger.debug('Creating {}.'.format(combined))
        #                 except Exception:
        #                     sys.exit(2)
        #             elif os.path.isfile(combined):
        #                 # This would be a lot faster if we didn't have to keep
        #                 # loading the combined trajectory for each traj
        #                 c = md.load(combined, top=top)
        #                 c = c.join(t)
        #                 c.save_dcd(combined)
        #         except:
        #             pass
        #         elapsed_time = time.time() - start_time
        #         logger.debug(('Elapsed time for {} was {:.2f} seconds').
        #                      format(traj.split('/')[-1], elapsed_time))
        # elif len(traj_list) is 0:
        #     logger.warn('Directory {} does not contain any trajectory files.'
        #                 .format(dir))
        #     logger.warn('Skipping {}'.format(dir))

        if len(traj_list) is not 0:
            logging.info('Concatenating trajectory files under {} into {}.'
                         .format(dir, combined))
            start_time = time.time()
            try:
                t = md.load(traj_list, top=top, stride=stride,
                            atom_indices=indices)
                if not os.path.isfile(combined):
                    try:
                        t.save_dcd(combined)
                        logger.debug('Combined output file {} not found.'
                                        .format(combined))
                        logger.debug('Creating {}.'.format(combined))
                    except Exception:
                        sys.exit(2)
                elif os.path.isfile(combined):
                    # This would be a lot faster if we didn't have to keep
                    # loading the combined trajectory for each traj
                    c = md.load(combined, top=top)
                    c = c.join(t)
                    c.save_dcd(combined)
            except:
                pass
            elapsed_time = time.time() - start_time
            logger.debug(('Elapsed time for {} was {:.2f} seconds').
                            format(traj_list, elapsed_time))
                            # format(traj.split('/')[-1], elapsed_time))

        elif len(traj_list) is 0:
            logger.warn('Directory {} does not contain any trajectory files.'
                        .format(dir))
            logger.warn('Skipping {}'.format(dir))

if __name__ == '__main__':
    main()
