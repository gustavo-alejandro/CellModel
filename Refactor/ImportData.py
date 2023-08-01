# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 19:21:11 2021

@author: Mark.Illingworth
"""

# Do imports first
import numpy as np
import pandas as pd
from Refactor.Properties import ValidProperty, CellProperties
from Refactor.Model import VariableVector, ControlMatrix, clone_vector
from Refactor.CellVoltage import VIcell
from Refactor.UKFClass import DiscreteModel, UKF, ukf_unpack, ukf_zpack


'''
# Establish Simulation and trnasformation matrices
'''


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
