# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 08:33:10 2021

# File to explore the sensitivity y of VCell to state variable
# perturbation both ddirectly and via transformation matrix
@author: Mark.Illingworth
"""

# Do imports first
# import Global as Glob
import Model as Mod
import numpy as np
# import InputSimulation as InpSim
import ObservationSimulation as ObSim
# import StateSimulation as StateSim
# import UnscentedKalmanFilterSim as UKFSim
import matplotlib.pyplot as plt


'''
# Establish Simulation and trnasformation matrices
'''

# only going to do a single timestep, so make it decent one
dt = 10
# examine the increments on 60 increments
nt = 60
cdmax = 8    # maximum alumina content of 8%
cdmin = 2    # minimum alumina content of 2%
cunmax = 6   # maximum undisolved alumina content of 6%
cunmin = 0   # minimum undisolved alumina content of 0%


(F, G, H, B) = Mod.StateEquations(dt)

# initial vectors - note the intention is to study x so
# we will use a fixed/constant u
x0 = np.array([[3.5], [0.5], [2.9193]])
u0 = np.array([[126000], [0.0], [0.0]])

# initialise u array to be a copy of u0
dims = np.shape(u0)
udim = dims[0]
uconst = np.zeros((nt+1, udim, 1))
uconst[0] = u0
# Initialise x arrays to explore the x1 and x2 characteristics
dims = np.shape(x0)
xdim = dims[0]
# independent flex over x1 for default x2 and vice versa
x1low = np.copy(x0)
x1low[0] = cdmin
x2low = np.copy(x0)
x2low[1] = cunmin
x1 = np.zeros((nt+1, xdim, 1))
x1dt = np.zeros((nt+1, xdim, 1))
x1[0] = x1low
x1dt[0] = F @ x1[0] + B @ uconst[0]
x2 = np.zeros((nt+1, xdim, 1))
x2dt = np.zeros((nt+1, xdim, 1))
x2[0] = x2low
x2dt[0] = F @ x2[0] + B @ uconst[0]
# independent increments over x1 and x2 span
x1inc = np.array([[(cdmax-cdmin)/nt], [0.0], [0.0]])
x2inc = np.array([[0.0], [(cunmax-cunmin)/nt], [0.0]])
for i in range(nt):
    # for each increment generate span in x1 or x2
    x1[i+1] = x1[i] + x1inc
    x2[i+1] = x2[i] + x2inc
    # hold inputs constant across the spans
    uconst[i+1] = u0
    # apply transformation for one timestep across each span
    x1dt[i+1] = F @ x1[i+1] + B @ uconst[i+1]
    x2dt[i+1] = F @ x2[i+1] + B @ uconst[i+1]

# transform state-spans to z with no noise and no plot
z1 = ObSim.obs(x1, uconst, H, 0, dt, False)
z1dt = ObSim.obs(x1dt, uconst, H, 0, dt, False)
z2 = ObSim.obs(x2, uconst, H, 0, dt, False)
z2dt = ObSim.obs(x2dt, uconst, H, 0, dt, False)
# differential on one timestep
z1delta = z1dt-z1
z2delta = z2dt-z2

fig, dissolved = plt.subplots()
dissolved.plot(x1[:, 0, 0], z1[:, 0, 0], label='Instantaneous')
dissolved.plot(x1dt[:, 0, 0], z1dt[:, 0, 0], label='60 sec timestep')
dissolved.set_xlabel('Dissolved Alumina Content (wt%)')
dissolved.set_ylabel('Vcell (Volts)')
dissolved.set_title("VCell Dependency on Dissolved Alumina (Undissolved 0.5%)")
difcolour = 'tab:red'
disdelta = dissolved.twinx()
disdelta.plot(x1[:, 0, 0], z1delta[:, 0, 0], label='Differential', color = difcolour)
dissolved.legend(loc='center')
disdelta.legend(loc='center right')

fig, undissolved = plt.subplots()
undissolved.plot(x2[:, 1, 0], z2[:, 0, 0], label='Instantaneous')
undissolved.plot(x2dt[:, 1, 0], z2dt[:, 0, 0], label='60 sec timestep')
undissolved.set_xlabel('UnDissolved Alumina Content (wt%)')
undissolved.set_ylabel('Vcell (Volts)')
undissolved.set_title("VCell Dependency - UnDissolved Alumina (Dissolved 3.5%)")
difcolour = 'tab:red'
undisdelta = undissolved.twinx()
undisdelta.plot(x2[:, 1, 0], z2delta[:, 0, 0], label='Differential', color = difcolour)
undissolved.legend(loc='upper center')
undisdelta.legend(loc='upper right')
