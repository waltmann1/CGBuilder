from __future__ import division
from CGBuilder.CGMolyAbs import CGMolyAbs
from CGBuilder.Irsp53Sh3 import Irsp53Sh3WithLinker
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class IBarIrsp53(CGMolyAbs):

    def __init__(self):
        super(IBarIrsp53, self).__init__()
        n_beads = 120
        connection1 = 1
        connection2 = 120
        self.connection_indices = [connection1, connection2]

        mass = [550 for _ in range(n_beads)]
        self.site_indexes = [[] for _ in range(n_beads)]
        self.f_weights = [[1] for _ in range(n_beads)]
        self.x_weights = [[10] for _ in range(n_beads)]
        self.name = "IBarIrsp53"
        #self.positions = [[0.0, 0.0, 0.0] for _ in range(n_beads)]
        self.positions = self.read_positions(self.abs_path("../../data/IBar_positions.txt"))
        self.atom_types = [[i+1, mass[i]] for i in range(n_beads)]
        ave_pos = np.mean(self.positions, axis=0)
        self.positions = np.subtract(self.positions, ave_pos)
        self.get_bonds(self.abs_path("../../data/IBar_cghenm.txt"))
        self.ibar_atom_types = n_beads
        self.ibar_bond_types = len(self.get_unique_bonds())

    def add_Sh3s(self, connection_indices, linker_length=16, uvecs=None):

        if uvecs is None:
            uvecs = [[1,0,0] for _ in range(len(connection_indices))]
        else:
            uvecs = [np.divide(u, np.linalg.norm(u)) for u in uvecs]
        for ind, connect_index in enumerate(connection_indices):
            sh3_link = Irsp53Sh3WithLinker(linker_length=linker_length, uvec=uvecs[ind])
            self.connect_sh3(cp.deepcopy(sh3_link), connect_index, uvec=uvecs[ind])


    def connect_sh3(self, sh3, connect_index, uvec=[1,0,0]):

        connect_pos = np.add(self.positions[connect_index], np.multiply(uvec, 5.0))
        vec = np.subtract(sh3.positions[0], connect_pos)
        sh3.positions = list(np.add(sh3.positions, -vec))
        self.site_indexes.extend(sh3.site_indexes)
        self.f_weights.extend(sh3.f_weights)
        self.x_weights.extend(sh3.x_weights)
        self.bond_types = list(self.bond_types)
        self.bonds = list(self.bonds)
        self.positions = list(self.positions)
        self.bonds.extend([[connect_index, len(self.positions) + 1]])
        self.bond_types.extend([sh3.bond_types[0]])
        #print(sh3.bond_types[0])
        #quit()
        sh3.bonds = np.add(sh3.bonds, len(self.positions))
        self.positions.extend(sh3.positions)
        self.name += "LinkedToSh3"
        for i in range(len(sh3.atom_types)):
            sh3.atom_types[i][0] += self.ibar_atom_types
        for i in range(len(sh3.bond_types)):
            sh3.bond_types[i][0] += self.ibar_bond_types
        self.atom_types.extend(sh3.atom_types)
        self.bond_types.extend(sh3.bond_types)
        self.bonds.extend(sh3.bonds)

class FullIBarIrsp53(IBarIrsp53):

    def __init__(self, connection_indices=[28, 90], uvecs= [[1, 1, 2], [-1, -1, 2]]):
        super(FullIBarIrsp53, self).__init__()
        self.add_Sh3s(connection_indices, uvecs=uvecs)






