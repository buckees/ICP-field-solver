# -*- coding: utf-8 -*-
"""
ICP Field Solver
Ampere's Law
"""

import numpy as np
from math import pi, sin, cos
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm

from Constants import MU0
from Mesh import MESHGRID

class wire_bf(object):
    """B-field object for infinite wire in 2D"""
    def __init__(self, loc=np.array([0.0, 0.0]), I=1.0):
        self.loc = loc # wire location
        self.curr = abs(I) # wire current
        self.sign = np.sign(I)
    
    def calc_bf(self, posn):
        vec = posn - self.loc
        dist = np.sqrt(np.sum(np.power(vec, 2)))
        uvec = vec/dist
        uvec = vec_rotate(uvec, self.sign*pi/2.0)
        bf = MU0*self.curr/(2.0*pi*dist)
        return bf, uvec
    
    def calc_bf_mat(self, X, Y):
        distx = X - self.loc[0]
        disty = Y - self.loc[1]
        dist, uvecx, uvecy = vec_norm_mat(distx, disty)
        uvecx, uvecy = mat_rotate(uvecx, uvecy,
                                  self.sign*pi/2.0)
        bf = MU0*self.curr/(2.0*pi)    
        bf = np.divide(bf, dist, where=dist>0, 
                       out=np.zeros_like(dist))
        return bf, uvecx, uvecy

    def plot_bf_contour(self, X, Y):
        fig, ax = plt.subplots(1, 1, figsize=(3.5,3))
        ax.plot(X, Y, marker='.', markersize=1,
                color='grey', linestyle='None')
        # you can manually set the levels
        lev_exp = np.arange(np.floor(np.log10(self.bf[np.nonzero(self.bf)].min())),
                            np.ceil(np.log10(self.bf.max())), 0.2)
        levs = np.power(10, lev_exp)
        cs = ax.contour(X, Y, self.bf,
                           levs, norm=colors.LogNorm())
#        ax.clabel(cs, cs.levels)
        fig.colorbar(cs)
        ax.quiver(X, Y, self.uvecx, self.uvecy)
    
    def plot_bf_stream(self, X, Y):
        fig, ax = plt.subplots(1, 1, figsize=(3,3))
        ax.plot(X, Y, marker='.', markersize=1,
                color='grey', linestyle='None')
        ax.streamplot(X, Y, self.uvecx, self.uvecy)

def vec_norm_mat(x, y):        
    dist = np.sqrt(np.power(x, 2) + np.power(y, 2))
    vecx = np.divide(x, dist, where=dist > 0.0,
                      out=np.zeros_like(dist))
    vecy = np.divide(y, dist, where=dist > 0.0,
                      out=np.zeros_like(dist))
    return dist, vecx, vecy


def vec_rotate(vec, theta):
    """
    counter-clock wise rotating vect by angle
    """
    rot_M = np.array([[cos(theta), -sin(theta)],
                      [sin(theta),  cos(theta)]])
    return np.matmul(rot_M, vec)

def mat_rotate(vecx, vecy, angle):
    """
    Rotate a mat counterclockwise by a given angle
    The angle should be given in radians.
    """
    vecx_rot = cos(angle)*vecx - sin(angle)*vecy
    vecy_rot = sin(angle)*vecx + cos(angle)*vecy
    return vecx_rot, vecy_rot

if __name__ == '__main__':
    wire1 = wire_bf(np.array([-3.0, 0.0]), 1.0)   
    width, height, nx, ny = 10.0, 10.0, 51, 51
    mesh = MESHGRID(width, height, nx, ny)
    mesh.init_mesh()
    
    wire1.bf, wire1.uvecx, wire1.uvecy = \
        wire1.calc_bf_mat(mesh.posx, mesh.posy)
    wire1.plot_bf_contour(mesh.posx, mesh.posy)
    wire1.plot_bf_stream(mesh.posx, mesh.posy)
    
    
