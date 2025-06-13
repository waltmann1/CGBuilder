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
        #self.positions = [[0.0, 0.0, 0.0] for _ in range(len(domains))]
        self.positions = self.read_positions(self.abs_path("../../data/sh3_positions.txt"))
        self.atom_types = [[i+1, mass[i]] for i in range(18)]
        ave_pos = np.mean(self.positions, axis=0)
        self.positions = np.subtract(self.positions, ave_pos)
        self.get_bonds(self.abs_path("../../data/sh3_cghenm.txt"))

    def add_linker(self, linker_length=16, uvec=[1,0,0]):
        linker_mass = [550 for _ in range(linker_length)]
        linker_site_indexes = [[] for _ in range(linker_length)]
        linker_f_weights = [[1] for _ in range(linker_length)]
        self.name += "WithLinker"
        linker_positions = [np.multiply(uvec, 5.0 * i) for i in range(linker_length)]
        linker_end = np.multiply(uvec, 5.0 * linker_length)
        linker_positions = np.subtract(linker_positions, linker_end)
        linker_positions = list(np.add(linker_positions, self.positions[0]))
        linker_positions.extend(self.positions)
        self.positions = linker_positions
        linker_atom_types = [[1, linker_mass[i]] for i in range(linker_length)]
        linker_bond_matrix = [[i + 1, i + 2] for i in range(linker_length)]
        r0 = 5
        k = 5
        linker_bonds = [[1, r0, k] for _ in range(len(linker_bond_matrix))]
        self.site_indexes.extend(linker_site_indexes)
        self.f_weights.extend(linker_f_weights)
        self.bonds = np.add(self.bonds, linker_length)
        linker_bond_matrix.extend(self.bonds)
        self.bonds = linker_bond_matrix
        for i in range(len(self.bond_types)):
            self.bond_types[i][0] = int(self.bond_types[i][0] + 1)
        for i in range(len(self.atom_types)):
            self.atom_types[i][0] = int(self.atom_types[i][0] + 1)
        linker_bonds.extend(self.bond_types)
        linker_atom_types.extend(self.atom_types)
        self.bond_types = linker_bonds
        self.atom_types = linker_atom_types


class Irsp53Sh3WithLinker(Irsp53Sh3):

    def __init__(self, linker_length=16, uvec= [1,0,0]):
        super(Irsp53Sh3WithLinker, self).__init__()
        self.add_linker(linker_length=linker_length, uvec=uvec)


