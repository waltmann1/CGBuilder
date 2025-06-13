import yaml 
import numpy as np
import sys
from CGBuilder.Gag import Gag
from CGBuilder.CGMolyComplex import CGMolyComplex
import copy as cp
from CGBuilder.GrimeLipid import GrimeLipid
from CGBuilder.IBarIrsp53 import FullIBarIrsp53
from CGBuilder.IBarIrsp53 import IBarIrsp53
from CGBuilder.OneBeadMoly import OneBeadMoly
from CGBuilder.BasicPolymer import BasicPolymer


fibar = FullIBarIrsp53()
gag = Gag()
ip6 = OneBeadMoly()
rna = BasicPolymer(100, 6.8)
gl =GrimeLipid()

CA_one = CGMolyComplex([gag, ip6,gl,rna ],numbers=[162, 41,5408 ,1],  box = [1500, 1500, 1500])

CA_one.get_master_positions_from_lammpstrj("dump.lammpstrj")

gag_positions = [CA_one.get_moly_positions(0,i) for i in range(162)]
ip6_positions = [CA_one.get_moly_positions(1,i) for i in range(41)]
gl_positions = [CA_one.get_moly_positions(2,i) for i in range(5408)]
rna_positions = [CA_one.get_moly_positions(3,i) for i in range(1)]
elements = [gag_positions, ip6_positions, gl_positions, rna_positions]

for e in elements:
    for g_pos in e:
       for apos in g_pos:
           apos[2] = -apos[2]
           if apos[1] < 190:
               apos[1]+=CA_one.box[1]


gl = GrimeLipid()
spacing = 15
patch_size = spacing * 75
mid = patch_size/2
lattice = [[x,y] for x in range(0,patch_size, spacing) for y in range(0, patch_size, spacing)]


CA_one = CGMolyComplex([gag, ip6, gl,rna, fibar],numbers=[162,41,2*len(lattice), 1, 1],  box = [patch_size, patch_size, 400])

CA_one.reset_master_positions()

for i in range(162):
    CA_one.set_moly_positions(gag_positions[i], 0, i)
    CA_one.shift_moly(-gag_positions[0][0] , 0, i)
    CA_one.shift_moly([mid, mid, 95], 0, i)

for i in range(41):
    CA_one.set_moly_positions(ip6_positions[i], 1, i)
    CA_one.shift_moly(-gag_positions[0][0] , 1, i)
    CA_one.shift_moly([mid, mid, 95], 1, i)

for ind, point in enumerate(lattice):
        x,y = point
        CA_one.shift_moly([x,y,21], 2, ind)
        CA_one.shift_moly([x, y, 21], 2, ind + len(lattice))
        CA_one.apply_lipid_flip(2, ind + len(lattice))

for i in range(1):
    CA_one.set_moly_positions(rna_positions[i], 3, i)
    CA_one.shift_moly(-gag_positions[0][0] , 3, i)
    CA_one.shift_moly([mid, mid, 95], 3, i)

CA_one.shift_moly([mid * 0.5, mid*1.5, 76], 4, 0)


CA_one.write_data_file("test.data")
CA_one.write_lammpstrj("test.lammpstrj")
