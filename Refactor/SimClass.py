# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 13:14:00 2021

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


# plot simulation of u true
def plot_usim(ulist, u, dt, maxplot=100):

    dims = np.shape(u)
    tsteps = dims[0]
    order = dims[1]
    order = min(order, maxplot)
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


class SimCell:

    def __init__(self, SimModel: DiscreteModel,
                 Cell: CellProperties, dT: float,
                 shotstd=0.0, Istd=0.0, statenoise=True):

        sim_valid = SimModel.defined and (dT > 0.0)
        if sim_valid is False:
            if SimModel.defined is False:
                raise ValueError("Discrete Process Model not defined")
            else:
                raise ValueError("Time Step Not Valid")
        else:
            # we can validly initialise the simulator
            self.SimModel = SimModel
            self.Cell = Cell
            self.inputcnt = np.shape(SimModel.u0.vector)[0]
            self.order = np.shape(SimModel.x0.vector)[0]
            # all models must have concentration and ACD
            try:
                self.IndexcAl2O3 = SimModel.x0.dictionary['cAl2O3']
            except KeyError:
                raise ValueError("Model does not contain cAl2O3")
            try:
                self.IndexACD = SimModel.x0.dictionary['ACD']
            except KeyError:
                raise ValueError("Model does not contain ACD")
            try:
                self.IndexuAl2O3 = SimModel.x0.dictionary['uAl2O3']
            except KeyError:
                # if the model doesnt have undissolved,
                # ok to point it at dissolved
                self.IndexuAl2O3 = self.IndexcAl2O3
            try:
                self.IndexI = SimModel.u0.dictionary['I']
                self.InputIFlag = True
                self.Inom = SimModel.u0.variable('I')
            except KeyError:
                # the input vector doesnt contain total current so
                # assuming noise term is in a deltaI state
                self.InputIFlag = False
                try:
                    self.IndexI = SimModel.u0.dictionary['I0']
                    self.Inom = SimModel.u0.variable('I0')
                except KeyError:
                    raise ValueError("Model does not contain I or deltaI")
            try:
                self.IndexMove = SimModel.u0.dictionary['B']
            except KeyError:
                raise ValueError("Model does not contain cAl2O3")
            try:
                self.IndexFeed = SimModel.u0.dictionary['g']
            except KeyError:
                self.IndexFeed = SimModel.u0.dictionary['g0']
            # set up dimensions for noise vector and extract standard
            # deviations from the noise variance matrix Q
            self.Qsd = np.sqrt(SimModel.Q(dT))
            self.noisecnt = np.shape(self.Qsd)[0]
            self.sim_valid = sim_valid
            self.shotstd = shotstd
            self.Istd = Istd
            self.statenoise = statenoise

    def __valid_simrqst(self, dT, u0: VariableVector, x0: VariableVector):

        rqst_valid = self.sim_valid
        if rqst_valid is True:
            # a zero or negative timestep is invalid
            dT_valid = (dT > 0.0)
            # The variable vector definition must match that used to set model
            state_dict_valid = (x0.dictionary == self.SimModel.x0.dictionary)
            # The variable vector definition must match that used to set model
            input_dict_valid = (u0.dictionary == self.SimModel.u0.dictionary)
            # Sim initialisation only valid if all conditions met
            rqst_valid = (rqst_valid and dT_valid and state_dict_valid and
                          input_dict_valid)
        if rqst_valid is False:
            if self.sim_valid is False:
                raise ValueError("Discrete Process Model not defined")
            elif dT_valid is False:
                raise ValueError("Time Step Not Valid")
            else:
                raise ValueError("State or Input dictionary mismatch")
        else:
            # update the stored model to reflect the requested initial state
            self.SimModel.x0.setvector(x0.vector)
            self.SimModel.u0.setvector(u0.vector)
        return rqst_valid

    def __Inoise(self, stdev):

        if stdev is not None:
            self.Istd = stdev
        noise = np.zeros((self.inputcnt, 1))
        if (self.InputIFlag is True) and (self.Istd > 0.0):
            # Only impose noise on inputs if total current is input
            noise[self.IndexI, 0] = rd.gauss(0, self.Istd)
        return noise

    def __xnoise(self, noiseflag):

        if noiseflag is not None:
            self.statenoise = noiseflag
        noise = np.zeros((self.noisecnt, 1))
        for n in range(self.noisecnt):
            # calculate guassian for each element of noise vector
            if self.statenoise is True:
                noise[n, 0] = rd.gauss(0, self.Qsd[n, n])
            else:
                # leave the noise as zero if statenoise is OFF
                pass
        return noise

    def __beammove(self, movedist):

        move = np.zeros((self.inputcnt, 1))
        if movedist is not None:
            move[self.IndexMove, 0] = movedist
        return move

    def __input_default(self):

        # need to define default input with nominal current and
        # zero feed and beam move to use in case the initial
        # simulation has a feed or a move
        u_default = self.SimModel.u0.vector
        u_default[self.IndexI] = self.Inom
        u_default[self.IndexFeed] = 0.0
        u_default[self.IndexMove] = 0.0
        return u_default

    def __nomfeed(self, Inom, shotstd):

        if Inom is not None:
            # If nominal current has changed, update the attribute
            self.Inom = Inom
        if shotstd is not None:
            # If shot variation has changed, update the attribute
            self.shotstd = shotstd
        # establish nominal feed cycle time
        NomShot = self.Cell.shotSetPoint.value
        TrueShot = self.Cell.shotMass.value
        NomFeedCycle, validcycle = self.Cell.calc_feed_cycle(self.Inom)
        if validcycle is False:
            # shouldnt get this far, will have error earlier
            raise ValueError("Value out of range")
        return NomShot, TrueShot, NomFeedCycle

    def __calc_next_state(self, xi, xnoise, ui, dT):

        xip1 = (self.SimModel.F(dT) @ xi +
                self.SimModel.B(dT) @ ui +
                self.SimModel.G(dT) @ xnoise)
        # some state variables can't actually be < 0.0
        if xip1[self.IndexcAl2O3] < 0.0:
            xip1[self.IndexcAl2O3] = 0.0
        if xip1[self.IndexuAl2O3] < 0.0:
            xip1[self.IndexuAl2O3] = 0.0
        if xip1[self.IndexACD] < 0.0:
            xip1[self.IndexACD] = 0.0
        return xip1

    def __plotsim(self, plotflag, x, u, dT):

        if plotflag is True:
            xlist = list(self.SimModel.x0.dictionary.keys())
            plot_xsim(xlist, x, dT)
            ulist = list(self.SimModel.u0.dictionary.keys())
            plot_usim(ulist, u, dT)

    def sim_nofeed(self, dT: float, duration: float,
                   u0: VariableVector, x0: VariableVector, *,
                   Istd=None, Inom=None, plotresult=True, statenoise=None):

        rqst_valid = self.__valid_simrqst(dT, u0, x0)
        if rqst_valid is True:
            if Inom is not None:
                self.Inom = Inom
            # establish number of timesteps for duration of nofeeding
            rawsteps = duration/dT
            # need to add 1 step to account for step zero
            Tsteps = int(ma.ceil(rawsteps))+1
            # set up the input vector as a 3d array
            u = np.zeros((Tsteps, self.inputcnt, 1))
            x = np.zeros((Tsteps, self.order, 1))
            u_def = self.__input_default()
            u[0] = u_def
            x[0] = self.SimModel.x0.vector
            # simulation
            for i in range(Tsteps - 1):
                # generate the input noise for amperage
                unoise = self.__Inoise(Istd)
                u[i+1] = u_def + unoise
                # for each step, initialise noise vector
                xnoise = self.__xnoise(statenoise)
                # then calculate the next simualted state
                x[i+1] = self.__calc_next_state(x[i], xnoise, u[i], dT)
            # plot the results if requested
            self.__plotsim(plotresult, x, u, dT)
        return (x, u)

    # Simulation of in the inputs for nominal feeding
    def sim_nominalfeed(self, dT: float, NomDuration: float,
                        u0: VariableVector, x0: VariableVector, *,
                        shotstd=None, Istd=None, Inom=None, BeamMove=None,
                        plotresult=True, statenoise=None):

        rqst_valid = self.__valid_simrqst(dT, u0, x0)
        if rqst_valid is True:
            # resolve the shot mass and feed cycle time including any update
            # to the nominal current that might apply
            NomShot, TrueShot, NomFeedCycle = self.__nomfeed(Inom, shotstd)
            # calculate duration and number of tsteps required
            # round up to whole nubmer of nominal cycles
            NomCycles = int(ma.ceil(NomDuration/NomFeedCycle))
            duration = NomCycles*NomFeedCycle
            rawsteps = duration/dT
            # need to add 1 step to account for step zero
            Tsteps = int(ma.ceil(rawsteps))+1
            # set up the input vector as a 3d array
            unom = np.zeros((Tsteps, self.inputcnt, 1))
            u = np.zeros((Tsteps, self.inputcnt, 1))
            x = np.zeros((Tsteps, self.order, 1))
            u_def = self.__input_default()
            # do a small up movement at the start of each nominal feed period
            beammove = self.__beammove(BeamMove)
            unom[0] = u_def + beammove
            u[0] = u_def + beammove
            x[0] = self.SimModel.x0.vector
            feedcount = 0
            for i in range(Tsteps - 1):
                # need to use i+2, assume last overfeed happened at Tstep ~ -1
                remainder, rqst = ma.modf((i+2)*dT/NomFeedCycle)
                feedtrue = np.zeros((self.inputcnt, 1))
                feednom = np.zeros((self.inputcnt, 1))
                if rqst > feedcount:
                    # if a feed is due then create a feed event
                    feedcount = feedcount+1
                    if self.shotstd > 0.0:
                        nShot = rd.gauss(0, self.shotstd)
                    else:
                        nShot = 0.0
                    # note the "true" feed uses the actual shot mass with noise
                    # but the filter will use the assumed/nominal shot mass
                    feedtrue[self.IndexFeed, 0] = TrueShot + nShot
                    feednom[self.IndexFeed, 0] = NomShot
                # generate the input noise for amperage
                unoise = self.__Inoise(Istd)
                u[i+1] = u_def + unoise + feedtrue
                unom[i+1] = u_def + unoise + feednom
                # for each step, initialise noise vector
                xnoise = self.__xnoise(statenoise)
                # then calculate the next simualted state
                x[i+1] = self.__calc_next_state(x[i], xnoise, u[i], dT)
            # plot the results if requested
            self.__plotsim(plotresult, x, u, dT)
        return (x, u, unom)

    # Simulation of the inputs and states for overfeeding feeding
    def sim_overfeed(self, dT: float, OFMass: float, OFRate: float,
                     u0: VariableVector, x0: VariableVector, *,
                     shotstd=None, Istd=None, Inom=None,
                     plotresult=True, statenoise=None):

        rqst_valid = self.__valid_simrqst(dT, u0, x0)
        if rqst_valid is True:
            # resolve the shot mass and feed cycle time including any update
            # to the nominal current that might apply
            NomShot, TrueShot, NomFeedCycle = self.__nomfeed(Inom, shotstd)
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
            unom = np.zeros((Tsteps, self.inputcnt, 1))
            u = np.zeros((Tsteps, self.inputcnt, 1))
            x = np.zeros((Tsteps, self.order, 1))
            u_def = self.__input_default()
            unom[0] = u_def
            u[0] = u_def
            x[0] = self.SimModel.x0.vector
            # Simulation
            feedcount = -1
            for i in range(Tsteps - 1):
                # need to use i+2 assume that first feed happens at Tstep 1
                remainder, rqst = ma.modf((i+2)*dT/OFCycle)
                feedtrue = np.zeros((self.inputcnt, 1))
                feednom = np.zeros((self.inputcnt, 1))
                if ((rqst > feedcount) and
                    (feedcount < cycles)) or ((i+3) == Tsteps):
                    # if a feed is due then create a feed event
                    feedcount = feedcount+1
                    if self.shotstd > 0.0:
                        nShot = rd.gauss(0, self.shotstd)
                    else:
                        nShot = 0.0
                    # note the "true" feed uses the actual shot mass with noise
                    # but the filter will use the assumed/nominal shot mass
                    feedtrue[self.IndexFeed, 0] = TrueShot + nShot
                    feednom[self.IndexFeed, 0] = NomShot
                # generate the input noise for amperage
                unoise = self.__Inoise(Istd)
                u[i+1] = u_def + unoise + feedtrue
                unom[i+1] = u_def + unoise + feednom
                # for each step, initialise noise vector
                xnoise = self.__xnoise(statenoise)
                # then calculate the next simualted state
                x[i+1] = self.__calc_next_state(x[i], xnoise, u[i], dT)
            # plot the results if requested
            self.__plotsim(plotresult, x, u, dT)
        return (x, u, unom)

    # Simulation of the inputs and states for a full OF-Nom-UF cycle
    def sim_adaptfeedcycle(self, dT: float, OFMass: float, OFRate: float,
                           NomDuration: float, UFDuration: float,
                           u0: VariableVector, x0: VariableVector, *,
                           shotstd=None, Istd=None, Inom=None, BeamMove=None,
                           plotresult=True, statenoise=None, NumCycles=1):
        # for multiple cycles, need to progressively capture the initial state
        x0i = clone_vector(x0)
        for i in range(NumCycles):
            # no need to check validity as it is done in each sub-call
            # Overfeed
            xOF, uOF, unomOF = self.sim_overfeed(dT, OFMass, OFRate, u0, x0i,
                                                 shotstd=shotstd, Istd=Istd,
                                                 statenoise=statenoise,
                                                 plotresult=False)
            # all the return arrays must have same number of steps
            # the next x0 vectors will be the last obs of the previous
            x0nom = clone_vector(x0, newvector=xOF[-1])
            xNom, uNom, unomNom = self.sim_nominalfeed(dT, NomDuration, u0,
                                                       x0nom, shotstd=shotstd,
                                                       Istd=Istd,
                                                       BeamMove=BeamMove,
                                                       statenoise=statenoise,
                                                       plotresult=False)
            x0UF = clone_vector(x0, newvector=xNom[-1])
            xUF, uUF = self.sim_nofeed(dT, UFDuration, u0, x0UF,
                                       Istd=Istd, statenoise=statenoise,
                                       plotresult=False)
            # concatenate but drop the last observation on the join as it is
            # repeated as the first observation in the next
            ui = np.concatenate((uOF[:-1], uNom[:-1], uUF), axis=0)
            xi = np.concatenate((xOF[:-1], xNom[:-1], xUF), axis=0)
            unomi = np.concatenate((unomOF[:-1], unomNom[:-1], uUF), axis=0)
            x0i = clone_vector(x0, newvector=xi[-1])
            if (i == 0):
                u = ui
                x = xi
                unom = unomi
            else:
                u = np.concatenate((u[:-1], ui), axis=0)
                x = np.concatenate((x[:-1], xi), axis=0)
                unom = np.concatenate((unom[:-1], unomi), axis=0)
        # plot the results if requested
        self.__plotsim(plotresult, x, u, dT)
        return (x, u, unom)

    # simulation of z
    def sim_observation(self, dT, x, u, plotresult=True, simnoise=True):

        u0 = clone_vector(self.SimModel.u0, newvector=u[0])
        x0 = clone_vector(self.SimModel.x0, newvector=x[0])
        rqst_valid = self.__valid_simrqst(dT, u0, x0)
        if rqst_valid is True:
            # initialise the observation array to the correct size
            Tsteps = np.shape(x)[0]
            obscnt = np.shape(self.SimModel.z0.vector)[0]
            z = np.zeros((Tsteps, obscnt, 1))
            InvalidObs = 0
            # set up dimensions for noise vector and extract standard
            # deviations from the noise variance matrix R
            Rsd = np.sqrt(self.SimModel.R(dT))
            noisedim = np.shape(Rsd)[0]
            # simulation
            for i in range(Tsteps):
                # unpack each state simulation and input into variable vector
                xk = clone_vector(self.SimModel.x0, newvector=x[i])
                uk = clone_vector(self.SimModel.u0, newvector=u[i])
                cAl2O3, ACD, I0, deltaI = ukf_unpack(xk, uk)
                # calculate measurement noise
                znoise = np.zeros((noisedim, 1))
                for n in range(noisedim):
                    # only apply noise if requested (default is True)
                    if simnoise is True:
                        # calculate guassian for each element of noise vector
                        znoise[n, 0] = rd.gauss(0, Rsd[n, n])
                # non-linear function VIcell(x,u,...) as there is no H matrix
                Vcell, Icell, Vvalid = VIcell(cAl2O3, ACD, deltaI, I0,
                                              self.Cell)
                if (Vvalid is False):
                    ValidSim = False
                    if (InvalidObs == 0):
                        InvalidObs = i*dT
                # repack the simulated z in correct format
                zk = ukf_zpack(Vcell, Icell, self.SimModel.z0)
                # add the noise terms to reflect measurement error
                z[i] = zk.vector + znoise
            if (plotresult is True) or (ValidSim is False):
                zlist = list(zk.dictionary.keys())
                plot_zsim(zlist, z, dT, InvalidObs)
        return(z)

    def sim_transform(self, x: np.array, u: np.array, z: np.array,
                      Inom: float, ratefactor: float,
                      UKFModel: DiscreteModel):

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
            (self.order == 3) and (UKFModel.order <= 4)):
            # criteria to do transform is met
            for key in UKFModel.x0.dictionary:
                keyfound = True
                ukf_idx = UKFModel.x0.dictionary[key]
                try:
                    sim_idx = self.SimModel.x0.dictionary[key]
                except KeyError:
                    keyfound = False
                if keyfound is True:
                    # straight copy if variable matches
                    xt[:, ukf_idx, :] = x[:, sim_idx, :]
                elif key == 'deltaI':
                    # delta I is the variation from supplied nominal current
                    sim_idx = self.SimModel.u0.dictionary['I']
                    xt[:, ukf_idx, :] = (u[:, sim_idx, :] - Inom)
                else:
                    # delta g rate term assumed to be left
                    sim_idx = self.SimModel.x0.dictionary['uAl2O3']
                    xt[:, ukf_idx, :] = (x[:, sim_idx, :] / ratefactor)
            for key in UKFModel.u0.dictionary:
                keyfound = True
                ukf_idx = UKFModel.u0.dictionary[key]
                try:
                    sim_idx = self.SimModel.u0.dictionary[key]
                except KeyError:
                    keyfound = False
                if keyfound is True:
                    # straight copy if variable matches
                    ut[:, ukf_idx, :] = u[:, sim_idx, :]
                elif key == 'I0':
                    # insert nominal current for all time steps
                    ut[:, ukf_idx, :] = Inom
                elif key == 'g0':
                    # map the feed mass as nominal feed mass g0
                    sim_idx = self.SimModel.u0.dictionary['g']
                    ut[:, ukf_idx, :] = u[:, sim_idx, :]
                elif key == 'ACD':
                    # looks like we have ACD as input not state
                    sim_idx = self.SimModel.x0.dictionary['ACD']
                    ut[:, ukf_idx, :] = x[:, sim_idx, :]
                else:
                    validtfr = False
            for key in UKFModel.z0.dictionary:
                keyfound = True
                ukf_idx = UKFModel.z0.dictionary[key]
                try:
                    sim_idx = self.SimModel.z0.dictionary[key]
                except KeyError:
                    keyfound = False
                if keyfound is True:
                    # straight copy if variable matches (Vcell)
                    zt[:, ukf_idx, :] = z[:, sim_idx, :]
                else:
                    # measured I is mapped from the input
                    sim_idx = self.SimModel.u0.dictionary['I']
                    zt[:, ukf_idx, :] = u[:, sim_idx, :]
        else:
            validtfr = False
        return validtfr, xt, ut, zt
