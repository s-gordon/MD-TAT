#!/usr/bin/env python
# AUTHOR:       Shane Gordon
# CREATED:      2015-06-16 21:46:32

import mdtraj as md
import numpy as np


def compute_rmsd(fname, topname, sel="name CA", step=1):
    rmsd = []
    atom_indices = md.load(topname).topology.select(sel)
    top = md.load(topname)
    for chunk in md.iterload(fname, top=top, stride=step):
        rmsd.append(md.rmsd(chunk, top, 0,
                            atom_indices=atom_indices))
    rmsd = np.concatenate(rmsd)
    return rmsd
