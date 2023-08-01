# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 16:16:26 2021

@author: Mark.Illingworth
"""

# Do imports first
import math as ma
import Global as Glob
import numpy as np
import random as rd


# Simulation of in the inputs for zero feed underfeed
def uNoFeeding(u0, Shotmass, Shotstd, Istd, dt, tsteps):

    # set up the input vector as a 3d array
    dims = np.shape(u0)
    udim = dims[0]
    u = np.zeros((tsteps, udim, 1))
    u[0] = u0
    for i in range(tsteps - 1):
        nI = rd.gauss(0, Istd)
        noise = np.array([[nI], [0], [0]])
        u[i+1] = u0 + noise
    return(u)


# Simulation of in the inputs for nominal feeding
def uNominalFeeding(u0, Shotmass, Shotstd, Istd, dt, cycles):

    # establish nominal feed cycle time
    Inom = u0[0]
    NomCE = Glob.CE
    NomShot = Glob.AlDump
    nomfeedcycle = Glob.FeedCycleCalc(Inom, NomCE, NomShot)
    # calculate duration and number of tsteps required
    duration = cycles*nomfeedcycle
    rawsteps = duration/dt
    # need to add 1 step to account for step zero
    tsteps = int(ma.ceil(rawsteps))+1
    # set up the input vector as a 3d array
    dims = np.shape(u0)
    udim = dims[0]
    utrue = np.zeros((tsteps, udim, 1))
    unominal = np.zeros((tsteps, udim, 1))
    # do a small up movement at the start of each nominal feed
    beammove = np.array([[0], [0], [0.02]])
    utrue[0] = u0 + beammove
    unominal[0] = u0+ beammove
    feedcount = 0
    for i in range(tsteps - 1):
        # need to use i+2 as we assume that last feed happened at step 1 = -1
        remainder, rqst = ma.modf((i+2)*dt/nomfeedcycle)
        if rqst > feedcount:
            # if a feed is due then create a feed event
            feedcount = feedcount+1
            nShot = rd.gauss(0, Shotstd)
            # note the "true" feed uses the actual shot mass
            # but the filter will use the assumed/nominal shot ,mass
            feedtrue = np.array([[0], [Shotmass+nShot], [0]])
            feednominal = np.array([[0], [NomShot+nShot], [0]])
        else:
            # if a feed is not due then create a zero feed event
            feedtrue = np.array([[0], [0], [0]])
            feednominal = np.array([[0], [0], [0]])
        # generate the input noise for amperage
        nI = rd.gauss(0, Istd)
        noise = np.array([[nI], [0], [0]])
        # the next input vector for the "true" and "nominal" cases
        utrue[i+1] = u0 + noise + feedtrue
        unominal[i+1] = u0 + noise + feednominal
    return(utrue, unominal)

# Simulation of in the inputs for overfeeding
def uOverFeeding(u0, Shotmass, Shotstd, Istd, dt, OFMass, OFRate):

    # establish nominal feed cycle time
    Inom = u0[0]
    NomCE = Glob.CE
    NomShot = Glob.AlDump
    nomfeedcycle = Glob.FeedCycleCalc(Inom, NomCE, NomShot)
    # Calculate the overfeed cycle and number of cycles to deliver the mass
    ofcycle = nomfeedcycle/OFRate
    # theoretical duration if feeding is continuous
    contduration = (OFMass-NomShot)*nomfeedcycle/((OFRate-1)*NomShot)
    cycles = int(ma.floor(contduration/ofcycle))
    remaining = contduration - cycles*ofcycle
    massremaining = remaining*NomShot/nomfeedcycle*(OFRate-1)
    finalcycle = (NomShot-massremaining)/NomShot*nomfeedcycle
    # calculate duration and number of tsteps required
    duration = cycles*ofcycle + finalcycle
    rawsteps = duration/dt
    # need to add 1 step to account for step zero
    tsteps = int(ma.ceil(rawsteps))+1
    # print(contduration, cycles, remaining, massremaining, finalcycle, duration, tsteps)
    # set up the input vector as a 3d array
    # create two instances, one with the actual feed mass
    # that is unknown, versus one with the nominal or assumed
    # feed mass for the model
    dims = np.shape(u0)
    udim = dims[0]
    utrue = np.zeros((tsteps, udim, 1))
    unominal = np.zeros((tsteps, udim, 1))
    # overfeed starts with a feed
    nShot = rd.gauss(0, Shotstd)
    feedtrue = np.array([[0], [Shotmass+nShot], [0]])
    feednominal = np.array([[0], [NomShot+nShot], [0]])
    utrue[0] = u0 + feedtrue
    unominal[0] = u0 + feednominal
    feedcount = 0
    for i in range(tsteps - 1):
        # need to use i+2 as we assume that last feed happened at step 1 = -1
        remainder, rqst = ma.modf((i+2)*dt/ofcycle)
        if ((rqst > feedcount) and (feedcount < cycles)) or ((i+3) == tsteps):
            # if a feed is due then create a feed event
            # note that the last step is also a feed
            feedcount = feedcount + 1
            # print(i, feedcount)
            nShot = rd.gauss(0, Shotstd)
            feedtrue = np.array([[0], [Shotmass+nShot], [0]])
            feednominal = np.array([[0], [NomShot+nShot], [0]])
        else:
            feedtrue = np.array([[0], [0], [0]])
            feednominal = np.array([[0], [0], [0]])
        nI = rd.gauss(0, Istd)
        noise = np.array([[nI], [0], [0]])
        utrue[i+1] = u0 + noise + feedtrue
        unominal[i+1] = u0 + noise + feednominal
    return(utrue, unominal)