# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 17:37:51 2021

@author: Mark.Illingworth
"""

# Do imports first
import Global as Glob
import Model as Mod
import numpy as np
import InputSimulation as InpSim
import ObservationSimulation as ObSim
import StateSimulation as StateSim
import UnscentedKalmanFilterSim as UKFSim


'''
# Establish Simulation and trnasformation matrices
'''

# 20 minutes simulation at 5 sec intervals
dt = 5
nt = 240
(F, G, H, B) = Mod.StateEquations(dt)

# Define the UKF noise standard deviations
q1dstd = 0.0002    # alumina concentration variation
q1ustd = 0.0002    # undisolved alumina concentration variation
q2std = 0.0000001  # acd varation needs to be much lower
# array of state variation terms
qstd = np.array([[q1dstd], [q1ustd], [q2std]])
rstd = 0.001      # measurement variation (volts)
Istd = 100         # Line Current input variation (amps)
Shotstd = 0.1      # Shot mass standard deviation (kg)
# print(qstd, rstd, Istd)

# start with initial concentrations of:
# 3.5% dissolved
# 0.5% undissolved
# starting ACD of 29.2 mm
x0 = np.array([[3.5], [0.5], [2.9193]])
u0 = np.array([[126000], [0.0], [0.0]])

# simulate inputs
ut = InpSim.uNoFeeding(u0, Glob.AlDump, Shotstd, Istd, dt, nt)
# simulate states, plot the results
xt = StateSim.xtrue(F, G, B, qstd, x0, ut, dt, True)

# Simulate the measurements and plot the results
zt = ObSim.obs(xt, ut, H, rstd, dt, True)

# Run the UKF SImulation
xest, Pest, xpred, Ppred, dz, S = UKFSim.ukf(F, G, H, B, qstd, rstd, zt, ut, x0, 2)

UKFSim.plot_ukf(xt, zt, dz, xest, xpred, dt, S)
