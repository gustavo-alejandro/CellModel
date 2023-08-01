# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 13:33:56 2021

@author: Mark.Illingworth
"""

import numpy as np
import random as rd
import math as ma
from Refactor.Properties import ValidProperty, CellProperties
from Refactor.Model import VariableVector, ControlMatrix, clone_vector
from Refactor.CellVoltage import VIcell
from Refactor.UKFClass import DiscreteModel, UKF, ukf_unpack, ukf_zpack
import matplotlib.pyplot as plt


def sim_transform(x: np.array, u: np.array, z: np.array,
                  Inom: float, ratefactor :float,
                  SimModel: DiscreteModel, UKFModel: DiscreteModel):

    # specific fucntion to transform a simulation by one model to
    # match the format of the UKF model
    # at this stage unique/specific to the 3 state model for simualtion
    # and 4 state model for UKF
    validtfr = True
    steps = np.shape(x)[0]
    xvar = np.size(UKFModel.x0.vector)
    uvar = np.size(UKFModel.u0.vector)
    zvar = np.size(UKFModel.z0.vector)
    xt = np.zeros((steps, xvar, 1), dtype=float)
    ut = np.zeros((steps, uvar, 1), dtype=float)
    zt = np.zeros((steps, zvar, 1), dtype=float)
    if ((np.shape(x)[0] == np.shape(u)[0] == np.shape(z)[0]) and
        (SimModel.order == 3) and (UKFModel.order == 4)):
        # criteria to do transform is met
        for key in UKFModel.x0.dictionary:
            keyfound = True
            ukf_idx = UKFModel.x0.dictionary[key]
            try:
                sim_idx = SimModel.x0.dictionary[key]
            except KeyError:
                keyfound = False
            if keyfound is True:
                # straight copy if variable matches
                xt[:, ukf_idx, :] = x[:, sim_idx, :]
            elif key == 'deltaI':
                # delta I is the variation from supplied nominal current
                sim_idx = SimModel.u0.dictionary['I']
                xt[:, ukf_idx, :] = (u[:, sim_idx, :] - Inom)
            else:
                # delta g rate term assumed to be left
                sim_idx = SimModel.x0.dictionary['uAl2O3']
                xt[:, ukf_idx, :] = (x[:, sim_idx, :] / ratefactor)
        for key in UKFModel.u0.dictionary:
            keyfound = True
            ukf_idx = UKFModel.u0.dictionary[key]
            try:
                sim_idx = SimModel.u0.dictionary[key]
            except KeyError:
                keyfound = False
            if keyfound is True:
                # straight copy if variable matches
                ut[:, ukf_idx, :] = u[:, sim_idx, :]
            elif key == 'I0':
                # insert nominal current for all time steps
                ut[:, ukf_idx, :] = Inom
            else:
                # map the feed mass as nominal feed mass g0
                sim_idx = SimModel.u0.dictionary['g']
                ut[:, ukf_idx, :] = u[:, sim_idx, :]
        for key in UKFModel.z0.dictionary:
            keyfound = True
            ukf_idx = UKFModel.z0.dictionary[key]
            try:
                sim_idx = SimModel.z0.dictionary[key]
            except KeyError:
                keyfound = False
            if keyfound is True:
                # straight copy if variable matches (Vcell)
                zt[:, ukf_idx, :] = z[:, sim_idx, :]
            else:
                # measured I is mapped from the input
                sim_idx = SimModel.u0.dictionary['I']
                zt[:, ukf_idx, :] = u[:, sim_idx, :]
    else:
        validtfr = False
    return validtfr, xt, ut, zt


# plot simulation of u true
def plot_usim(ulist, u, dt):

    dims = np.shape(u)
    tsteps = dims[0]
    order = dims[1]
    t = np.linspace(0, dt*tsteps, tsteps)
    for i in range(order):
        plt.figure(figsize=(8, 6), dpi=128)
        plt.plot(t, u[:, i, 0])
        plt.xlabel('time (sec)')
        plt.ylabel(ulist[i])
        plt.title(f'Simulated Input Variable: {ulist[i]}')


# plot simulation of x true
def plot_xsim(xlist, x, dt):

    dims = np.shape(x)
    tsteps = dims[0]
    order = dims[1]
    t = np.linspace(0, dt*tsteps, tsteps)
    for i in range(order):
        plt.figure(figsize=(8, 6), dpi=128)
        plt.plot(t, x[:, i, 0])
        plt.xlabel('time (sec)')
        plt.ylabel(xlist[i])
        plt.title(f'Simulated State Variable: {xlist[i]}')


# plot simulation of z true
def plot_zsim(zlist, z, dt, InvalidObs):

    dims = np.shape(z)
    tsteps = dims[0]
    order = dims[1]
    t = np.linspace(0, dt*tsteps, tsteps)
    for i in range(order):
        plt.figure(figsize=(8, 6), dpi=128)
        plt.plot(t, z[:, i, 0])
        plt.xlabel('time (sec)')
        plt.ylabel(zlist[i])
        title = f'Simulated Measurement Variable: {zlist[i]}'
        if InvalidObs > 0.0:
            title = title + f', Invalid at t = {InvalidObs} sec'
            plt.axvline(InvalidObs)
        plt.title(title)


# simulation of z
def sim_observation(dT, x, u, Cell: CellProperties, Model: DiscreteModel,
                    plotresult=True):
    sim_valid = Model.defined
    if sim_valid is True:
        # a zero or negative timestep is invalid
        dT_valid = (dT > 0.0)
        # Sim initialisation only valid if all conditions met
        sim_valid = (sim_valid and dT_valid)
    if sim_valid is False:
        if Model.defined is False:
            raise ValueError("Discrete Process Model not defined")
        else:
            raise ValueError("Time Step Not Valid")
    else:
        # initialise the observation array to the correct size
        dims = np.shape(x)
        tsteps = dims[0]
        dims = np.shape(Model.z0.vector)
        zdim = dims[0]
        z = np.zeros((tsteps, zdim, 1))
        # Valid simulation - check there are no invalid flags returned
        ValidSim = True
        InvalidObs = 0
        # set up dimensions for noise vector and extract standard
        # deviations from the noise variance matrix Q
        Rsd = np.sqrt(Model.R(dT))
        noisedim = np.shape(Rsd)[0]
        # simulation
        for i in range(tsteps):
            # unpack each state simulation and input into variable vector form
            xk = clone_vector(Model.x0, newvector=x[i])
            uk = clone_vector(Model.u0, newvector=u[i])
            cAl2O3, ACD, I0, deltaI = ukf_unpack(xk, uk)
            # calculate measurement noise
            znoise = np.zeros((noisedim, 1))
            for n in range(noisedim):
                # calculate guassian for each element of noise vector
                znoise[n, 0] = rd.gauss(0, Rsd[n, n])
            # non-linear function VIcell(x,u,...) as there is no H matrix
            Vcell, Icell, Vvalid = VIcell(cAl2O3, ACD, deltaI, I0, Cell)
            if (Vvalid is False):
                ValidSim = False
                if (InvalidObs == 0):
                    InvalidObs = i*dT
            # repack the simulated z in correct format
            zk = ukf_zpack(Vcell, Icell, Model.z0)
            # add the noise terms to reflect measurement error
            z[i] = zk.vector + znoise
        if (plotresult is True) or (ValidSim is False):
            zlist = list(Model.z0.dictionary.keys())
            plot_zsim(zlist, z, dT, InvalidObs)
    return(z)


def sim_validate_model(dT, Model: DiscreteModel,
                       u0: VariableVector, x0: VariableVector):

    sim_valid = Model.defined
    if sim_valid is True:
        # a zero or negative timestep is invalid
        dT_valid = (dT > 0.0)
        # The variable vector definition must match that used to set model
        state_dict_valid = (x0.dictionary == Model.x0.dictionary)
        # The variable vector definition must match that used to set model
        input_dict_valid = (u0.dictionary == Model.u0.dictionary)
        # Sim initialisation only valid if all conditions met
        sim_valid = (sim_valid and dT_valid and state_dict_valid and
                     input_dict_valid)
    if sim_valid is False:
        if Model.defined is False:
            raise ValueError("Discrete Process Model not defined")
        elif dT_valid is False:
            raise ValueError("Time Step Not Valid")
        else:
            raise ValueError("State or Input dictionary mismatch")
    return sim_valid


def sim_nofeed(dT, Tsteps, Cell: CellProperties, Model: DiscreteModel,
               u0: VariableVector, x0: VariableVector, Istd=0.0,
               plotresult=True, statenoise=True):

    sim_valid = sim_validate_model(dT, Model, u0, x0)
    if sim_valid is True:
        # set up the input vector as a 3d array
        ninput = np.shape(u0.vector)[0]
        order = np.shape(x0.vector)[0]
        u = np.zeros((Tsteps, ninput, 1))
        x = np.zeros((Tsteps, order, 1))
        u[0] = u0.vector
        x[0] = x0.vector
        # all models must have concentration and ACD
        IndexcAl2O3 = x0.dictionary['cAl2O3']
        IndexACD = x0.dictionary['ACD']
        try:
            IndexuAl2O3 = x0.dictionary['uAl2O3']
        except KeyError:
            # if the model doesnt have undissolved, ok to point it at dissolved
            IndexuAl2O3 = IndexcAl2O3
        try:
            indexI = u0.dictionary['I']
            InputI = True
        except KeyError:
            # the input vector doesnt contain total current so
            # assuming noise term is in a deltaI state
            InputI = False
            indexI = x0.dictionary['deltaI']
        # set up dimensions for noise vector and extract standard
        # deviations from the noise variance matrix Q
        Qsd = np.sqrt(Model.Q(dT))
        noisedim = np.shape(Qsd)[0]
        # simulation
        for i in range(Tsteps - 1):
            unoise = np.zeros((ninput, 1))
            if InputI is True:
                # Only impose noise on inputs if total current is input
                unoise[indexI, 0] = rd.gauss(0, Istd)
            u[i+1] = u[0] + unoise
            # for each step, initialise noise vector
            xnoise = np.zeros((noisedim, 1))
            for n in range(noisedim):
                # calculate guassian for each element of noise vector
                if statenoise is True:
                    xnoise[n, 0] = rd.gauss(0, Qsd[n, n])
                else:
                    # leave the noise as zero if statenoise is OFF
                    pass
            x[i+1] = (Model.F(dT) @ x[i] + Model.B(dT) @ u[i] +
                      Model.G(dT) @ xnoise)
            # none of the state variables can actually be < 0.0
            if x[i+1, IndexcAl2O3] < 0.0:
                x[i+1, IndexcAl2O3] = 0.0
            if x[i+1, IndexuAl2O3] < 0.0:
                x[i+1, IndexuAl2O3] = 0.0
            if x[i+1, IndexACD] < 0.0:
                x[i+1, IndexACD] = 0.0
            # for simulation purposes, state variable deltaI will be zeros
            # if we havent applied statenoise. If we want to impose gaussian
            # noise then we do it afterwards
            # if statenoise is False and InputI is False:
                # special case to impose variable current on
        if plotresult is True:
            xlist = list(x0.dictionary.keys())
            plot_xsim(xlist, x, dT)
            ulist = list(u0.dictionary.keys())
            plot_usim(ulist, u, dT)
    return (x, u)


# Simulation of in the inputs for nominal feeding
def sim_nominalfeed(dT, NomCycles, Cell: CellProperties, Model: DiscreteModel,
                    u0: VariableVector, x0: VariableVector, shotstd=0.0,
                    Istd=0.0, plotresult=True, statenoise=True):

    sim_valid = sim_validate_model(dT, Model, u0, x0)
    if sim_valid is True:
        # all models must have concentration and ACD and beam moves
        IndexcAl2O3 = x0.dictionary['cAl2O3']
        IndexACD = x0.dictionary['ACD']
        IndexMove = u0.dictionary['B']
        try:
            IndexuAl2O3 = x0.dictionary['uAl2O3']
        except KeyError:
            # if the model doesnt have undissolved, ok to point it at dissolved
            IndexuAl2O3 = IndexcAl2O3
        try:
            indexI = u0.dictionary['I']
            InputI = True
            # limitation of 3 state model - use initial I as nominal
            Inom = u0.variable('I')
        except KeyError:
            # the input vector doesnt contain total current so
            # assuming noise term is in a deltaI state
            InputI = False
            indexI = x0.dictionary['deltaI']
            Inom = u0.variable('I0')
        try:
            IndexFeed = u0.dictionary['g']
        except KeyError:
            IndexFeed = u0.dictionary['g0']
        # establish nominal feed cycle time
        NomShot = Cell.shotSetPoint.value
        TrueShot = Cell.shotMass.value
        NomFeedCycle, validcycle = Cell.calc_feed_cycle(Inom)
        if validcycle is False:
            # shouldnt get this far, will have error earlier
            raise ValueError("Value out of range")
        # calculate duration and number of tsteps required
        duration = NomCycles*NomFeedCycle
        rawsteps = duration/dT
        # need to add 1 step to account for step zero
        Tsteps = int(ma.ceil(rawsteps))+1
        # set up the input vector as a 3d array
        ninput = np.shape(u0.vector)[0]
        order = np.shape(x0.vector)[0]
        unom = np.zeros((Tsteps, ninput, 1))
        x = np.zeros((Tsteps, order, 1))
        u = np.zeros((Tsteps, ninput, 1))
        # do a small up movement at the start of each nominal feed period
        beammove = np.zeros((ninput, 1))
        beammove[IndexMove, 0] = 0.02
        uinit = u0.vector
        unom[0] = uinit + beammove
        u[0] = uinit + beammove
        x[0] = x0.vector
        # set up dimensions for noise vector and extract standard
        # deviations from the noise variance matrix Q
        Qsd = np.sqrt(Model.Q(dT))
        noisedim = np.shape(Qsd)[0]
        feedcount = 0
        for i in range(Tsteps - 1):
            # need to use i+2, assume last overfeed happened at Tstep ~ -1
            remainder, rqst = ma.modf((i+2)*dT/NomFeedCycle)
            feedtrue = np.zeros((ninput, 1))
            feednom = np.zeros((ninput, 1))
            if rqst > feedcount:
                # if a feed is due then create a feed event
                feedcount = feedcount+1
                nShot = rd.gauss(0, shotstd)
                # note the "true" feed uses the actual shot mass with noise
                # but the filter will use the assumed/nominal shot mass
                feedtrue[IndexFeed, 0] = TrueShot + nShot
                feednom[IndexFeed, 0] = NomShot
            # generate the input noise for amperage
            unoise = np.zeros((ninput, 1))
            if InputI is True:
                # Only impose noise on inputs if total current is input
                unoise[indexI, 0] = rd.gauss(0, Istd)
            # the next input vector for the "true" and "nominal" cases
            u[i+1] = uinit + unoise + feedtrue
            unom[i+1] = uinit + unoise + feednom
            # for each step, initialise state noise vector
            xnoise = np.zeros((noisedim, 1))
            for n in range(noisedim):
                # calculate guassian for each element of noise vector
                if statenoise is True:
                    xnoise[n, 0] = rd.gauss(0, Qsd[n, n])
                else:
                    # leave the noise as zero if statenoise is OFF
                    pass
            x[i+1] = (Model.F(dT) @ x[i] + Model.B(dT) @ u[i] +
                      Model.G(dT) @ xnoise)
            # none of these state variables can actually be < 0.0
            if x[i+1, IndexcAl2O3] < 0.0:
                x[i+1, IndexcAl2O3] = 0.0
            if x[i+1, IndexuAl2O3] < 0.0:
                x[i+1, IndexuAl2O3] = 0.0
            if x[i+1, IndexACD] < 0.0:
                x[i+1, IndexACD] = 0.0
            # for simulation purposes, state variable deltaI will be zeros
            # if we havent applied statenoise. If we want to impose gaussian
            # noise then we do it afterwards
            # if statenoise is False and InputI is False:
                # special case to impose variable current on
        if plotresult is True:
            xlist = list(x0.dictionary.keys())
            plot_xsim(xlist, x, dT)
            ulist = list(u0.dictionary.keys())
            plot_usim(ulist, u, dT)
    return (x, u, unom)


# Simulation of in the inputs and states for overfeeding feeding
def sim_overfeed(dT, OFMass, OFRate, Cell: CellProperties,
                 Model: DiscreteModel,
                 u0: VariableVector, x0: VariableVector, shotstd=0.0,
                 Istd=0.0, plotresult=True, statenoise=True):

    sim_valid = sim_validate_model(dT, Model, u0, x0)
    if sim_valid is True:
        # all models must have concentration and ACD and beam moves
        IndexcAl2O3 = x0.dictionary['cAl2O3']
        IndexACD = x0.dictionary['ACD']
        IndexMove = u0.dictionary['B']
        try:
            IndexuAl2O3 = x0.dictionary['uAl2O3']
        except KeyError:
            # if the model doesnt have undissolved, ok to point it at dissolved
            IndexuAl2O3 = IndexcAl2O3
        try:
            indexI = u0.dictionary['I']
            InputI = True
            # limitation of 3 state model - use initial I as nominal
            Inom = u0.variable('I')
        except KeyError:
            # the input vector doesnt contain total current so
            # assuming noise term is in a deltaI state
            InputI = False
            indexI = x0.dictionary['deltaI']
            Inom = u0.variable('I0')
        try:
            IndexFeed = u0.dictionary['g']
        except KeyError:
            IndexFeed = u0.dictionary['g0']
        # establish nominal feed cycle time
        NomShot = Cell.shotSetPoint.value
        TrueShot = Cell.shotMass.value
        NomFeedCycle, validcycle = Cell.calc_feed_cycle(Inom)
        if validcycle is False:
            # shouldnt get this far, will have error earlier
            raise ValueError("Value out of range")
        # calculate duration and number of tsteps required
        OFCycle = NomFeedCycle/OFRate
        # theoretical duration if feeding is continuous
        OFDrtnTheo = (OFMass-NomShot)*NomFeedCycle/((OFRate-1)*NomShot)
        cycles = int(ma.floor(OFDrtnTheo/OFCycle))
        remaining = OFDrtnTheo - cycles*OFCycle
        massremaining = remaining*NomShot/NomFeedCycle*(OFRate-1)
        finalcycle = (NomShot-massremaining)/NomShot*NomFeedCycle
        # calculate duration and number of tsteps required
        duration = cycles*OFCycle + finalcycle
        rawsteps = duration/dT
        # need to add 1 step to account for step zero
        Tsteps = int(ma.ceil(rawsteps))+1
        # set up the input vector as a 3d array
        ninput = np.shape(u0.vector)[0]
        order = np.shape(x0.vector)[0]
        unom = np.zeros((Tsteps, ninput, 1))
        x = np.zeros((Tsteps, order, 1))
        u = np.zeros((Tsteps, ninput, 1))
        unom[0] = u0.vector
        u[0] = u0.vector
        x[0] = x0.vector
        # set up dimensions for noise vector and extract standard
        # deviations from the noise variance matrix Q
        Qsd = np.sqrt(Model.Q(dT))
        noisedim = np.shape(Qsd)[0]
        feedcount = 0
        for i in range(Tsteps - 1):
            # need to use i+2 as we assume that first feed happens at Tstep 1
            remainder, rqst = ma.modf((i+2)*dT/OFCycle)
            feedtrue = np.zeros((ninput, 1))
            feednom = np.zeros((ninput, 1))
            if ((rqst > feedcount) and
                (feedcount < cycles)) or ((i+3) == Tsteps):
                # if a feed is due then create a feed event
                feedcount = feedcount+1
                nShot = rd.gauss(0, shotstd)
                # note the "true" feed uses the actual shot mass with noise
                # but the filter will use the assumed/nominal shot mass
                feedtrue[IndexFeed, 0] = TrueShot + nShot
                feednom[IndexFeed, 0] = NomShot
            # generate the input noise for amperage
            unoise = np.zeros((ninput, 1))
            if InputI is True:
                # Only impose noise on inputs if total current is input
                unoise[indexI, 0] = rd.gauss(0, Istd)
            # the next input vector for the "true" and "nominal" cases
            u[i+1] = u[0] + unoise + feedtrue
            unom[i+1] = u[0] + unoise + feednom
            # for each step, initialise state noise vector
            xnoise = np.zeros((noisedim, 1))
            for n in range(noisedim):
                # calculate guassian for each element of noise vector
                if statenoise is True:
                    xnoise[n, 0] = rd.gauss(0, Qsd[n, n])
                else:
                    # leave the noise as zero if statenoise is OFF
                    pass
            x[i+1] = (Model.F(dT) @ x[i] + Model.B(dT) @ u[i] +
                      Model.G(dT) @ xnoise)
            # none of these state variables can actually be < 0.0
            if x[i+1, IndexcAl2O3] < 0.0:
                x[i+1, IndexcAl2O3] = 0.0
            if x[i+1, IndexuAl2O3] < 0.0:
                x[i+1, IndexuAl2O3] = 0.0
            if x[i+1, IndexACD] < 0.0:
                x[i+1, IndexACD] = 0.0
            # for simulation purposes, state variable deltaI will be zeros
            # if we havent applied statenoise. If we want to impose gaussian
            # noise then we do it afterwards
            # if statenoise is False and InputI is False:
                # special case to impose variable current on
        if plotresult is True:
            xlist = list(x0.dictionary.keys())
            plot_xsim(xlist, x, dT)
            ulist = list(u0.dictionary.keys())
            plot_usim(ulist, u, dT)
    return (x, u, unom)
