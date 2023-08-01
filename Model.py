# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 11:00:49 2021

@author: Mark.Illingworth
"""

import Global as Glob
import numpy as np

# Use global defaults
kdiss = 0.002 # Glob.kdiss
MAlumina = Glob.MAlumina
CE = Glob.CE
Faraday = Glob.Faraday
z = Glob.z
m = Glob.m
r = 0.6 # Glob.r
dcon = Glob.dcon
daccum = Glob.daccum


def StateEquations(dt):

    F = np.array([[1.0, kdiss*dt, 0.0], [0.0, 1.0-kdiss*dt, 0.0], [0.0, 0.0, 1.0]])
    G = np.array([[dt, 0.0, 0.0], [0.0, dt, 0.0], [0.0, 0.0, dt]])
    H = np.array([1.0])
    B = np.array([[-MAlumina*CE*dt/(10*Faraday*z*m), 100*(1-r)/m, 0.0],
                 [0.0, 100*r/m, 0.0], [(dcon-daccum)*dt, 0.0, 1.0]])
    return (F, G, H, B)
