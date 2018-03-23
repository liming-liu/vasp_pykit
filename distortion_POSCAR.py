#!/usr/bin/env python
# To draw atomic structure using python
# author: Liming Liu
# version: v_1
# date: 2018-1-26

import numpy as np
import math as m
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as axes3d
from ase.io import read, write

# define colors
ele_color = {'C':'gray', 'H':'yellow', 'O':'red', 'Ti':'cyan'} # elemental colors
arr_color = {'C':'gray', 'H':'balck', 'O':'red', 'Ti':'cyan'} # bond colors

# read in POSCARS
pos_0 = read('POSCAR_0', format='vasp')
pos = read('POSCAR_1', format='vasp')
ele = pos.get_chemical_symbols() # all element name in POSCAR_1
sys_color = [ele_color[xx] for xx in ele] # sort element-color pair
num = pos.get_number_of_atoms()
ij = [(i, j)  for i in range(num) for j in range(i)   if (i != j) and
        pos.get_distance(i, j) < 2.3] # finde the adjacent atom pairs

# figure instance
fig = plt.figure(figsize=(8.0, 6.0))
# 3d axis instance
ax = plt.subplot(111, projection='3d')

# plot bonds
for i, j in ij:
    pos_i = pos.positions[i]
    pos_j = pos.positions[j]
    mid = (pos_i + pos_j)/2

    ele_i = pos.get_chemical_symbols()[i]
    ele_j = pos.get_chemical_symbols()[j]

    ax.plot((pos_i[0], mid[0]), (pos_i[1], mid[1]), (pos_i[2], mid[2]),
            color=ele_color[ele_i], lw=0.5)
    ax.plot((pos_j[0], mid[0]), (pos_j[1], mid[1]), (pos_j[2], mid[2]),
            color=ele_color[ele_j], lw=0.5)

# atom
p = pos.positions
ax.scatter(p[:,0], p[:,1], p[:,2], color=sys_color)

# --- quiver ---
x, y, z = p[:,0], p[:,1], p[:,2]
d = pos.positions - pos_0.positions
u, v, w = d[:,0], d[:,1], d[:,2]
ax.quiver(x, y, z, u, v, w)

# set view angle
ax.view_init(90, 90)

# save figure
plt.savefig('fig.png',dpi=300)
#plt.show()
