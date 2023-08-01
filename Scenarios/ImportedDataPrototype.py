# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 19:21:11 2021

@author: Mark.Illingworth
"""

# Do imports first
import Global as Glob
import Model as Mod
import numpy as np
import pandas as pd
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
rstd = 0.001       # measurement variation (volts)
Istd = 2.00        # Not used
Shotstd = 0.1      # Shot mass standard deviation (kg)
ShotOffset = -0.0  # Discrepancy between nominal and actual shot mass (kg)
OFmass = 24        # 25 kg overfeed mass in excess of consumption
OFrate = 1.5       # Overfeed delivered at 1.5 times nominal rate
NomCycles = 5      # The number of nominal feed cycles to simulate
AdaptiveCycles = 5  # The number of OF-Nom-UF cycles to simulate

# start with initial concentrations at end of previous underfeed of:
# 2.8% dissolved
# 0.05% undissolved
# starting ACD of 29.2 mm
x0 = np.array([[2.5], [0.1], [2.5193]])


def LoadData(filename: str):
    # This is a quick loader - not finalised
    celldata = pd.read_csv(filename, sep=',')
    npcelldata = celldata.to_numpy()
    dims = np.shape(npcelldata)
    tsteps = dims[0]
    udim = 3
    zdim = 1
    u = np.zeros((tsteps, udim, 1))
    z = np.zeros((tsteps, zdim, 1))
    for i in range(tsteps):
        VCell = npcelldata[i, 0]
        Iline = npcelldata[i, 1]
        feed = npcelldata[i, 3]
        beammove = npcelldata[i, 4]
        uobs = np.array([[Iline], [feed], [beammove]])
        u[i] = np.copy(uobs)
        zobs = np.array([[VCell]])
        z[i] = np.copy(zobs)
    return u, z

# load in actual data for input and observation
ut, zt = LoadData("C:/Python Files/AluminaConcEstimator/DataImport/excel_prototype.csv")
print(ut)
print(zt)

# simulate states with the true inputs, plot the results
xsim = StateSim.xtrue(F, G, B, qstd, x0, ut, dt, True)

# Run the UKF SImulation using the nominal inputs
xest, Pest, xpred, Ppred, dz, S = UKFSim.ukf(F, G, H, B, qstd, rstd, zt, ut, x0, 0)

UKFSim.plot_ukf(xsim, zt, dz, xest, xpred, dt, S)
