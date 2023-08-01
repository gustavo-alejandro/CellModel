# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 14:55:42 2021

@author: Mark.Illingworth
"""

import numpy as np
import random as rd
import matplotlib.pyplot as plt


# plot simulation of x true
def plot_xsim(x, dt):

    dims = np.shape(x)
    tsteps = dims[0]
    t = np.linspace(0, dt*tsteps, tsteps)
    fig, axconc = plt.subplots()
    axconc.plot(t, x[:, 0, 0], label='C_d')
    axconc.plot(t, x[:, 1, 0], label='C_un')
    axconc.set_xlabel('time')
    axconc.set_ylabel('wt% alumina')
    axconc.set_title("Simulated Dissolved and Undissolved Alumina")
    axconc.legend()
    fig, axacd = plt.subplots()
    axacd.plot(t, x[:, 2, 0], label='D')
    axacd.set_xlabel('time')
    axacd.set_ylabel('ACD (cm)')
    axacd.set_title("Simulated ACD")
    axacd.legend()


# Simulation of x
def xtrue(F, G, B, Q, x0, u, dt, plotresult):

    # set up the state vector as a 3d array
    dims = np.shape(x0)
    xdim = dims[0]
    # the number of timesteps will match the input array
    dims = np.shape(u)
    tsteps = dims[0]
    x = np.zeros((tsteps, xdim, 1))
    x[0] = x0
    # set up dimensions for noise vector
    dims = np.shape(Q)
    qdim = dims[0]
    # simulation
    for i in range(tsteps - 1):
        # for each step, initialise noise vector
        noise = np.zeros((qdim, 1))
        for n in range(qdim):
            # calculate guassian for each element of noise vector
            noise[n, 0] = rd.gauss(0, Q[n, 0])
        x[i+1] = F @ x[i] + B @ u[i] + G @ noise
        for n in range(xdim):
            # none of the state variables can actually be < 0.0
            if x[i+1, n] < 0.0:
                x[i+1, n] = 0.0
    if plotresult is True:
        plot_xsim(x, dt)
    return(x)
