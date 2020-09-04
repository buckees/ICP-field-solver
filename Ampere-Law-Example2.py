# -*- coding: utf-8 -*-
"""
ICP Field Solver
Ampere's Law
"""

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm

from Constants import MU0
from Mesh import MESHGRID

#set infinite-long wire to position (wx, wy) with infinitesimal radius
I = 1.0 # wire current in A
# According to Ampere's law in integral form
# B(r|r>r0) = mu0*I/(2*pi*r)
#The earth's magnetic field is about 0.5 gauss. 
width, height, nx, ny = 10.0, 10.0, 101, 101
mesh = MESHGRID(width, height, nx, ny)
mesh.init_mesh()

def calc_bf(position, I):
    dist, vecx, vecy = mesh.calc_dist(position, I)
    bf = MU0*I/(2.0*pi)
    dist_min = min(width/(nx-1), height/(ny-1))
    bf_max = np.ones_like(dist)*bf/dist_min
    
    bf = np.divide(bf, dist, where=dist>dist_min, out=bf_max)
    bf = abs(bf)
                 
    print('B field min = %.2e max = %.2e' % (bf.min(), bf.max()))
    
#    fig, ax = plt.subplots(figsize=(3,3))
#    ax.plot(mesh.posx, mesh.posy, '.k',
#            marker='.', markersize=3,
#            color='black', linestyle='None')
    #fmt = ticker.LogFormatterMathtext()
    #fmt.create_dummy_axis()
    #cs = ax.contour(mesh.posx, mesh.posy, bf,
    #                locator=ticker.LogLocator(subs=range(1,6)), 
    #                cmap=cm.plasma)
    # Alternatively, you can manually set the levels
    # and the norm:
#    lev_exp = np.arange(np.floor(np.log10(bf.min())),
#                        np.ceil(np.log10(bf.max())), 0.1)
#    levs = np.power(10, lev_exp)
#    cs = ax.contour(mesh.posx, mesh.posy, bf, levs, norm=colors.LogNorm())
    #ax.clabel(cs, cs.levels)
#    fig.colorbar(cs)
#    ax.quiver(mesh.posx, mesh.posy, vecx, vecy)
#    ax.plot(position[0], position[1],
#            color='red', marker='o', markersize=15)
    return bf, vecx, vecy

pos1, pos2 = (-1.5, 0.0), (1.5, 0.0)
bf1, vx1, vy1 = calc_bf(pos1, I)
bf2, vx2, vy2 = calc_bf(pos2, I)

vx = np.multiply(bf1, vx1) + np.multiply(bf2, vx2)

vy = np.multiply(bf1, vy1) + np.multiply(bf2, vy2)
bf = np.sqrt(np.power(vx, 2) + np.power(vy, 2))
print('B field min = %.2e max = %.2e' % (bf[np.nonzero(bf)].min(),
                                         bf.max()))
vx, vy = np.divide(vx, bf), np.divide(vy, bf)

fig, ax = plt.subplots(figsize=(3,3))
#ax.plot(mesh.posx, mesh.posy, '.k',
#        marker='.', markersize=3,
#        color='black', linestyle='None')
# Alternatively, you can manually set the levels
# and the norm:
lev_exp = np.arange(np.floor(np.log10(bf[np.nonzero(bf)].min())),
                    np.ceil(np.log10(bf.max())), 0.1)
levs = np.power(10, lev_exp)
#levs = np.linspace(bf.min(), bf.max(), 50)
cs = ax.contour(mesh.posx, mesh.posy, bf, levs, norm=colors.LogNorm())
#ax.clabel(cs, cs.levels)
#fig.colorbar(cs)
#ax.quiver(mesh.posx, mesh.posy, vx, vy)
#ax.plot(pos1[0], pos1[1],
#        color='red', marker='o', markersize=15)
#ax.plot(pos2[0], pos2[1],
#        color='red', marker='o', markersize=15)

#fig, ax = plt.subplots(2, 1, figsize=(6,6))
#ax[0].plot(mesh.posx[int(nx/2), :], bf[int(nx/2), :])
#ax[1].plot(mesh.posy[:, int(ny/2)], bf[:, int(ny/2)])


