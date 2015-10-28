#!/usr/bin/env python
# AUTHOR:       Shane Gordon
# CREATED:      2015-06-16 21:46:32

import mdtraj as md
import numpy as np


def compute_rg(fname, topname):
    rg = []
    for chunk in md.iterload(fname, top=topname, chunk=2):
        rg.append(md.compute_rg(chunk))
    rg = np.concatenate(rg)
    return rg
