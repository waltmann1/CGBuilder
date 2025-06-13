from __future__ import division
from CGBuilder.CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class CADimer(CGMolyAbs):

    def __init__(self):
        super(CADimer, self).__init__()
        mass = self.read_mass(self.abs_path("../../data/BottomUpCA_mass.txt"))
        self.site_indexes = [[] for _ in range(462)]
        self.f_weights = [[1] for _ in range(462)]
        self.x_weights = [[10] for _ in range(462)]
        self.name = "CADimer"
        self.positions = self.read_positions(self.abs_path("../../data/BottomUpCA_pos.txt"))
        self.atom_types = [[i+1, mass[i]] for i in range(231)]
        self.atom_types.extend(self.atom_types)
        ave_pos = np.mean(self.positions, axis=0)
        self.positions = np.subtract(self.positions, ave_pos)
        #self.positions = np.add(self.positions, [50, 50, 50])
        self.get_bonds(self.abs_path("../../data/BottomUpCA_henm.txt"))
        more_bonds = []
        for bond in self.bonds:
            new = cp.deepcopy(bond)
            new[0] += 231
            new[1] += 231
            more_bonds.append(new)
        np.array(list(self.bonds).extend(more_bonds))
        self.bonds = np.concatenate((self.bonds, more_bonds), axis=0)
        self.bond_types = np.concatenate((self.bond_types, self.bond_types), axis=0)
        temp = cp.deepcopy(self.bonds)
        temp_types = cp.deepcopy(self.bond_types)

        self.get_bonds(self.abs_path("../../data/BottomUpCA_inter_henm.txt"))
        #print(temp)
        for t in self.bond_types:
            t[0] += len(np.unique(temp_types, axis=0))
        self.bonds = np.concatenate((temp, self.bonds), axis=0)
        self.bond_types = np.concatenate((temp_types, self.bond_types), axis=0)


