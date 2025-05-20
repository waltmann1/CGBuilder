from __future__ import division
from CGBuilder.CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class GrimeLipid(CGMolyAbs):
    
    def __init__(self, lower=True):
        super(GrimeLipid, self).__init__()
        mass = [300 for _ in range(4)]
        types = [3,2,2,1]
        self.site_indexes = [[] for _ in range(4)]
        self.f_weights = [1 for _ in range(4)]
        self.x_weights = [10 for _ in range(4)]
        self.name = "GrimeLipid"
        if lower:
            self.positions = [[0.0, 0.0, -1 * 7 *i] for i in range(4)]
        else:
            self.positions= [[0.0, 0.0, 7 * (i+1)] for i in range(4)]
        self.atom_types = [[types[i], mass[i]] for i in range(4)]
        ave_pos = np.mean(self.positions, axis=0)
        self.positions = np.subtract(self.positions, ave_pos)
        self.get_bonds(self.abs_path("../../data/lipid_bonds.txt"))
        self.angles = [[i, i+1, i+2] for i in range(1,3)]
        theta = 180
        k = 0.5961
        self.angle_types = [[i+1, theta, k] for i in range(len(self.angles))]