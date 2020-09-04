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
from Analytic_Solution import wire_bf, vec_norm_mat
#
width, height, nx, ny = 10.0, 10.0, 21, 21
mesh = MESHGRID(width, height, nx, ny)
mesh.init_mesh()

#
wire1 = wire_bf(np.array([-3.0, 0.0]), 1.0)   
wire2 = wire_bf(np.array([+3.0, 0.0]), -1.0)   

wire1.bf, wire1.uvecx, wire1.uvecy = \
    wire1.calc_bf_mat(mesh.posx, mesh.posy)
#wire1.plot_bf(mesh.posx, mesh.posy)

wire2.bf, wire2.uvecx, wire2.uvecy = \
    wire2.calc_bf_mat(mesh.posx, mesh.posy)
#wire2.plot_bf(mesh.posx, mesh.posy)

wire3 = wire_bf()
wire3_x = np.multiply(wire1.bf, wire1.uvecx) + \
          np.multiply(wire2.bf, wire2.uvecx)
wire3_y = np.multiply(wire1.bf, wire1.uvecy) + \
          np.multiply(wire2.bf, wire2.uvecy)


wire3.bf, wire3.uvecx, wire3.uvecy = vec_norm_mat(wire3_x, 
                                                  wire3_y)

wire3.plot_bf_contour(mesh.posx, mesh.posy)
wire3.plot_bf_stream(mesh.posx, mesh.posy)


#fig, ax = plt.subplots(1, 1, figsize=(3.5,3))
#ax.quiver(mesh.posx, mesh.posy, wire3.uvecx, wire3.uvecy)
#
#ax.plot(wire1.loc[0], wire1.loc[1],
#        marker='o', markersize=5, color='black')
#ax.plot(wire2.loc[0], wire2.loc[1],
#        marker='o', markersize=5, color='black')

#posn = np.random.uniform(-5.0, 5.0, 2)
#bf1, uvec1 = wire1.calc_bf(posn)
#bf2, uvec2 = wire2.calc_bf(posn)
#
#bf3_vec = bf1*uvec1 + bf2*uvec2
#bf3 = np.sqrt(np.sum(np.power(bf3_vec, 2)))
#uvec3 = bf3_vec/bf3
#
#fig, ax = plt.subplots(1, 1, figsize=(3.5,3))
#ax.plot([wire1.loc[0], posn[0]], [wire1.loc[1], posn[1]], 
#        marker='o', markersize=5, color='red', linestyle='-')
#ax.plot([wire2.loc[0], posn[0]], [wire2.loc[1], posn[1]], 
#        marker='o', markersize=5, color='blue', linestyle='-')
#
#ax.plot(wire1.loc[0], wire1.loc[1],
#        marker='o', markersize=5, color='black')
#ax.plot(wire2.loc[0], wire2.loc[1],
#        marker='o', markersize=5, color='black')
#ax.quiver(posn[0], posn[1], uvec1[0], uvec1[1],
#          color='red')
#ax.quiver(posn[0], posn[1], uvec2[0], uvec2[1], 
#          color='blue')
#
#ax.quiver(posn[0], posn[1], uvec3[0], uvec3[1], 
#          color='green')
#
#ax.set_xlim(-5.0, 5.0)
#ax.set_ylim(-5.0, 5.0)
#
#
#fig, ax = plt.subplots(1, 1, figsize=(3.5,3))
#ax.plot(wire1.loc[0], wire1.loc[1],
#        marker='o', markersize=5, color='black')
#ax.plot(wire2.loc[0], wire2.loc[1],
#        marker='o', markersize=5, color='black')
#ax.streamplot(mesh.posx, mesh.posy, wire3.uvecx, wire3.uvecy,
#              color='k', linewidth=0.8, density=1.3,
#              minlength=0.9, arrowstyle='->')
