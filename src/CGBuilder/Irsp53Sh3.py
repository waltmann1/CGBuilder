from __future__ import division
from CGBuilder.CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class Irsp53Sh3(CGMolyAbs):

    def __init__(self):
        super(Irsp53Sh3, self).__init__()
        domains = [[1, 2], [45, 46, 47, 60, 61, 62], [10, 11, 12, 13, 34, 35, 36, 67, 68], [71], [14, 15, 16, 17, 30, 31, 32, 33], [63, 64, 65, 66], [25, 26, 27, 57, 58, 59], [3, 4], [21, 22, 23, 24], [50, 51, 52, 53, 54, 55, 56], [7, 8, 9], [40, 41, 42], [69, 70], [37, 38, 39, 48, 49], [18, 19, 20, 28, 29], [5, 6], [72], [43, 44]]
        mass = [110 * len(lizt) for lizt in domains]
        self.site_indexes = [[] for _ in range(18)]
        self.f_weights = [[1] for _ in range(18)]
        self.x_weights = [[10] for _ in range(18)]
        self.name = "Irsp53Sh3"
        self.positions = [[0.0, 0.0, 0.0] for _ in range(len(domains))]
        self.atom_types = [[i+1, mass[i]] for i in range(18)]
        ave_pos = np.mean(self.positions, axis=0)
        self.positions = np.subtract(self.positions, ave_pos)
        self.get_bonds(self.abs_path("../../data/sh3_cghenm.txt"))