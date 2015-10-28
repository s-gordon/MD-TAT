#!/usr/bin/env python
# AUTHOR:       Shane Gordon
# CREATED:      2015-06-16 21:46:32

import mdtraj as md
import numpy as np


def compute_rmsd(fname, topname):
    rmsd = []
    sel = 'name CA'
    atom_indices = md.load(topname).topology.select(sel)
    for chunk in md.iterload(fname, top=topname, chunk=2):
        rmsd.append(md.rmsd(chunk, md.load(topname), 0,
                            atom_indices=atom_indices))
    rmsd = np.concatenate(rmsd)
    return rmsd
