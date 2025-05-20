import numpy.linalg as la
import copy as cp
from CGBuilder.BasicPolymer import BasicPolymer
import math
class PolymerWithAngles(BasicPolymer):

    def __init__(self, nbeads, bond_length, mass=10, bond_k=10, separation=None, lattice_points=False, angle_k=None, angle_theta=None):
        super(BasicPolymer, self).__init__(nbeads, bond_length, mass=mass, k=bond_k, separation=separation, lattice_points=lattice_points)

        self.angles = [[i, i+1, i+2] for i in range(1, nbeads - 1)]
        theta = angle_theta
        k = angle_k
        self.angle_types = [[1, theta, k] for _ in range(len(self.angles))]


class RNAForGag(PolymerWithAngles):

    def __init__(self, nbeads=3000, lattice_points=False, mass=330):
        super(PolymerWithAngles, self).__init__(nbeads, 6.8, mass=mass, bond_k=2.5, separation=None, lattice_points=lattice_points,
                                                angle_theta=145, angle_k=1.0)

