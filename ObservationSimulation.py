# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:07:53 2021

@author: Mark.Illingworth
"""

# Do imports first
import Global as Glob
import numpy as np
import random as rd
import CellVoltage as Meas
import matplotlib.pyplot as plt


# plot simulation of z true
def plot_zsim(z, dt, InvalidObs):

    dims = np.shape(z)
    tsteps = dims[0]
    t = np.linspace(0, dt*tsteps, tsteps)
    fig, axvolt = plt.subplots()
    axvolt.plot(t, z[:, 0, 0], label='Vcell')
    axvolt.set_xlabel('time')
    axvolt.set_ylabel('Voltage (V)')
    if InvalidObs > 0.0:
        axvolt.set_title("Simulated Measured Cell Voltage - Invalid Results")
        axvolt.axvline(InvalidObs)
    else:
        axvolt.set_title("Simulated Measured Cell Voltage")
    axvolt.legend()


# simulation of z
def obs(x, u, H, R, dt, plotresult):
    # initialise the observation array to the correct size
    dims = np.shape(x)
    tsteps = dims[0]
    dims = np.shape([H])
    zdim = dims[0]
    z = np.zeros((tsteps, zdim, 1))
    # Valid simulation - check there are no invalid flags returned
    ValidSim = True
    InvalidObs = 0
    # simulation
    for i in range(tsteps):
        # unpack each state observation/simulation
        obs = x[i]
        inp = u[i]
        # calculate measurement noise
        noise = rd.gauss(0, R)
        # non-linear function Vcell(x,u,...) multiplied by a nominal H array
        # H is a single element array [[1.0]]
        V_cell, valid = Meas.Vcell(obs, inp, Glob.BathConstVector)
        if (valid is False):
            ValidSim = True
            if (InvalidObs == 0):
                InvalidObs = i*dt
        z[i] = H @ np.array([[V_cell]]) + np.array([[noise]])
    if (plotresult is True) or (ValidSim is False):
        plot_zsim(z, dt, InvalidObs)
    if (ValidSim is False):
        print("One or more simulated observations was out of range")
    return(z)
