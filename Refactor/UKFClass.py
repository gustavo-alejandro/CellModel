# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 17:02:48 2021

@author: Mark.Illingworth
"""

import numpy as np
import math as ma
from Refactor.Properties import ValidProperty, CellProperties
from Refactor.Model import VariableVector, ControlMatrix, clone_vector
from Refactor.CellVoltage import VIcell
from numpy.linalg import inv
from scipy.linalg import cholesky, LinAlgError
import matplotlib.pyplot as plt


# a fix because python does not do edge cases
def myminv(A):
    if np.size(A) == 1:
        return (1/A)
    else:
        return inv(A)


def ukf_gen_samples(xhat, P, gamma):

    # infer shape of state vector - order is number of states
    xdim = np.shape(xhat)
    pdim = np.shape(P)
    order = xdim[0]
    # dfine number of samples and set up arrays to return them
    nsamples = 2*order + 1
    xsample = np.zeros((nsamples, order, 1))
    weight = np.zeros((nsamples, 1, 1))
    # check dimensionality - if xhat and P must be 2-d
    valid = (len(xdim) == 2) and (len(pdim) == 2)
    # then check that xhat is 1 column and P is square of same size
    if valid is True:
        valid = valid and (xdim[1] == 1) and (pdim[0] == pdim[1] == order)
    # only if both checks pass do we attempt the calcs
    if valid is True:
        # evaluate sigmai from P matrix, note P must be positive definite
        P = P*(order + gamma)
        # note uses scipy cholesky which defaults to a upper triangular
        try:
            sigma = cholesky(P, lower=True)
        except LinAlgError:
            print("Not Positive Definite")
            valid = False
        if valid is True:
            # sample 0 is the input xhat mean, with weight as defined
            xsample[0] = xhat
            weight[0] = gamma/(order + gamma)
            for i in range(order):
                # extract each column from sigma
                sigmai = sigma[:, i:i+1]
                print(f'sigma[{i}]:')
                print(sigmai)
                # two samples for each column of sigma +/-
                xsample[2*i+1] = xhat + sigmai
                xsample[2*i+2] = xhat - sigmai
                # corresponding weights are the same
                weight[2*i+1] = 1/(2*(order+gamma))
                weight[2*i+2] = 1/(2*(order+gamma))
    return valid, xsample, weight


def ukf_unpack(xvect: VariableVector, uvect: VariableVector):

    # need to verfiy the contents of the vectors
    # all models should have cAl2O3 and ACD as variables of X
    cAl2O3 = xvect.variable('cAl2O3')
    try:
        ACD = xvect.variable('ACD')
    except KeyError:
        try:
            ACD = uvect.variable('ACD')
        except KeyError:
            ACD = 0.0
    try:
        deltaI = xvect.variable('deltaI')
    except KeyError:
        deltaI = 0.0
    try:
        Icell = uvect.variable('I')
    except KeyError:
        Icell = 0.0
    try:
        I0 = uvect.variable('I0')
    except KeyError:
        I0 = 0.0
    # if I0 is not in the vectors then we need to infer it
    if (I0 == 0.0):
        I0 = Icell - deltaI
    # note it is possible that both I0 and deltaI are set to zero
    # and the client will need to decife whether that is acceptable
    return cAl2O3, ACD, I0, deltaI


def ukf_zpack(volt: float, crrnt: float, ztemplate: VariableVector):

    # need to pack a z vector to match the template
    # start with a zero vector
    zpacked = clone_vector(ztemplate, zeros=True)
    try:
        zpacked.setvariable('Vcell', volt)
    except KeyError:
        pass
    try:
        zpacked.setvariable('Icell', crrnt)
    except KeyError:
        pass
    # note it is possible that both try statements fail
    # and the client will need to decife whether that is acceptable
    return zpacked


class DiscreteModel:

    def __init__(self, F: ControlMatrix, B: ControlMatrix,
                 G: ControlMatrix, Q: ControlMatrix,
                 R: ControlMatrix, z0: VariableVector,
                 u0: VariableVector, x0: VariableVector):

        # Check for consistent dimensionality first
        x_var_cnt = x0.vector.size
        u_var_cnt = u0.vector.size
        z_var_cnt = z0.vector.size
        # Note that we dont care about time step, only the size of the arrays
        F_row_cnt, F_col_cnt = np.shape(F.array_eval())
        B_row_cnt, B_col_cnt = np.shape(B.array_eval())
        G_row_cnt, G_col_cnt = np.shape(G.array_eval())
        Q_row_cnt, Q_col_cnt = np.shape(Q.array_eval())
        R_row_cnt, R_col_cnt = np.shape(R.array_eval())
        # F must be square matched with state vector (x) size
        F_valid = (F_row_cnt == F_col_cnt == x_var_cnt)
        # B must have row count matched with x size, col count matched with u
        B_valid = (B_row_cnt == x_var_cnt) and (B_col_cnt == u_var_cnt)
        # Q must be square and == state vector (x) size (shared not supported)
        Q_valid = (Q_row_cnt == Q_col_cnt == x_var_cnt)
        # R must be square and match the measure vector (z) size
        R_valid = (R_row_cnt == R_col_cnt == z_var_cnt)
        # G must have row count matched with x size, col count matched with Q
        G_valid = (G_row_cnt == x_var_cnt) and (G_col_cnt == Q_row_cnt)
        # variance matrices must be positive (sigma^2 elements)
        var_valid = not((Q.array_eval() < 0).any() or
                        (R.array_eval() < 0).any())
        init_valid = (F_valid and B_valid and Q_valid and R_valid
                      and G_valid and var_valid)
        if init_valid is False:
            # return a flag to indicate if the model was valid and initialised
            self.defined = False
            if F_valid is False:
                raise ValueError("F array not square or inconsistent with x")
            if B_valid is False:
                raise ValueError("B array inconsistent with x or u")
            if Q_valid is False:
                raise ValueError("Q array not square or inconsistent with x")
            if R_valid is False:
                raise ValueError("R array not square or inconsistent with z")
            if G_valid is False:
                raise ValueError("G array inconsistent with x or Q")
            if var_valid is False:
                raise ValueError("Q or R array contain negative variances")
        else:
            # lets initialise the Model with the control matrices private
            self._F = F    # discrete transformation matrix
            self._B = B    # discrete input/transformation matrix
            self._G = G    # discrete noise/transformation matrix
            self._Q = Q    # state noise terms (variances)
            self._R = R    # measurent noise (indep variances)
            self.x0 = clone_vector(x0)   # state vector definition
            self.u0 = clone_vector(u0)   # input vector definition
            self.z0 = clone_vector(z0)   # measurement vector definition
            self.order = x_var_cnt
            # return a flag to indicate if the model was valid and initialised
            self.defined = True

    def F(self, dT):

        # Public Method to evaluate the Transition Matrix F as a function of
        # the time step dT
        return self._F.array_eval(dT)

    def B(self, dT):

        # Public Method to evaluate the Input Matrix B as a function of
        # the time step dT
        return self._B.array_eval(dT)

    def G(self, dT):

        # Public Method to evaluate the Noise Trans Matrix G as a function of
        # the time step dT
        return self._G.array_eval(dT)

    def Q(self, dT=0.0):

        # Public Method to evaluate the State Noise Variance Matrix Q
        # less likely to have a dependency on timestep hence default is 0
        return self._Q.array_eval(dT)

    def R(self, dT=0.0):

        # Public Method to evaluate the Measurement Noise Variance Matrix R
        # less likely to have a dependency on timestep hence default is 0
        return self._R.array_eval(dT)


class UKF:

    def __init__(self, Model: DiscreteModel,
                 Cell: CellProperties, x0: VariableVector, dT: float, *,
                 xnorm=None, unorm=None, znorm=None):

        # the model must have been defined, if not then UKF not initialised
        UKF_valid = Model.defined
        if UKF_valid is True:
            # a zero or negative timestep is invalid
            dT_valid = (dT > 0.0)
            # The variable vector definition must match that used to set model
            state_dict_valid = (x0.dictionary == Model.x0.dictionary)
            # TODO proper handling of normalisation later
            if xnorm is None:
                xnorm = np.ones_like(x0.vector, dtype=float)
                xnorm_valid = True
            elif isinstance(xnorm, VariableVector) is False:
                xnorm_valid = False
            elif (Model.x0.dictionary == xnorm.dictionary):
                xnorm_valid = True
                xnorm = xnorm.vector
            else:
                xnorm_valid = False
            if unorm is None:
                unorm = np.ones_like(Model.u0.vector, dtype=float)
                unorm_valid = True
            elif isinstance(unorm, VariableVector) is False:
                unorm_valid = False
            elif (Model.u0.dictionary == unorm.dictionary):
                unorm_valid = True
                unorm = unorm.vector
            else:
                unorm_valid = False
            if znorm is None:
                znorm = np.ones_like(Model.z0.vector, dtype=float)
                znorm_valid = True
            elif isinstance(znorm, VariableVector) is False:
                znorm_valid = False
            elif (Model.z0.dictionary == znorm.dictionary):
                znorm_valid = True
                znorm = znorm.vector
            else:
                znorm_valid = False
            norm_valid = xnorm_valid and unorm_valid and znorm_valid
            # UKF initialisation only valid if all conditions met
            UKF_valid = (UKF_valid and dT_valid and state_dict_valid and
                         norm_valid)
        if UKF_valid is False:
            if Model.defined is False:
                raise ValueError("Discrete Process Model not defined")
            elif dT_valid is False:
                raise ValueError("Time Step Not Valid")
            elif state_dict_valid is False:
                raise ValueError("State dictionary mismatch")
            else:
                raise ValueError("Normalisation vector error/mismatch")

        else:
            # Once initialised, we force the use of the same model
            self.dT = None
            self.Model = Model
            self.Cell = Cell
            # TODO proper handling of normalisation later
            self.xnorm = xnorm
            self.unorm = unorm
            self.znorm = znorm
            # Flag attribute so UKF object tracks if it is doing
            # a prediction or an update. When the flag is True
            # if another prediction is called this signals that
            # the update/estimate is old so the prediction must
            # compound on the previous prediction
            self.predictionflag = False
            # default attribute for UKF sampling
            self.ukf_gamma = 1
            # initialise UKF Matrices based on initial time-step.
            updated = self.update_timestep(dT)
            if updated is True:
                # TODO proper handling of normalisation later
                # these are genuine initial conditions that get over-written
                x0n = x0.vector * self.xnorm
                self.xest = clone_vector(x0, newvector=x0n)
                # note these were originally 100 not 30
                self.Pest = 30*self.Qm
                self.xpred = clone_vector(x0, newvector=x0n)
                self.Ppred = 30*self.Qm
                self.uk = clone_vector(Model.u0)
                self.zk = clone_vector(Model.z0)
                self.determinant = np.linalg.det(self.Pest)
                self.condition_num = np.linalg.cond(self.Pest)
                self.initialised = True
            else:
                UKF_valid = False
                self.initialised = False
                raise ValueError("Time Step Update Failed on initialisation")

    def update_timestep(self, dT):

        # TODO proper handling of normalisation later
        # method used to verify if timestep has changed and update the
        # UKF model parameters appropriately
        if dT != self.dT:
            self.F = self.Model.F(dT) * self.xnorm / self.xnorm.T
            self.B = self.Model.B(dT) * self.xnorm / self.unorm.T
            self.G = self.Model.G(dT)
            Q = self.Model.Q() * self.xnorm.T * self.xnorm.T
            self.Qm = self.G @ Q @ self.G.T
            self.Rm = self.Model.R() * self.znorm.T * self.znorm.T
            self.dT = dT
            # print("F: ", self.F)
            # print("B: ", self.B)
            # print("G: ", self.G)
            # print("Qm: ", self.Qm)
            # print("Rm: ", self.Rm)
            updated = True
        else:
            updated = False
        return updated

    def predict_state(self, u: VariableVector, dT: float):

        # the UKF must have been initialised, if not then no prediction
        predict_valid = self.initialised
        if predict_valid is True:
            # a zero or negative timestep is invalid
            dT_valid = (dT > 0.0)
            # The variable vector definition must match that used to set model
            input_dict_valid = (u.dictionary == self.Model.u0.dictionary)
            # prediction is only valid if all conditions met
            predict_valid = predict_valid and dT_valid and input_dict_valid
        # only do a prediction if valid
        if predict_valid is True:
            # first check if that last action was a prediction
            # if it was then we need to compound the previous
            # prediction and use those as the local "estimates"
            if self.predictionflag is True:
                # compounded prediction
                xest = self.xpred.vector
                Pest = self.Ppred
            else:
                # new prediction
                xest = self.xest.vector
                Pest = self.Pest
            # store the latest input vector internally
            # TODO proper handling of normalisation later
            uk = u.vector * self.unorm
            self.uk.setvector(u.vector)
            # Update UKF Matrices if the next time-step has changed.
            updated = self.update_timestep(dT)
            if updated is False:
                # dont care about result in this instance as dT is valid
                pass
            # Calculate state prediction and prediction covariance
            xpred = self.F @ xest + self.B @ uk
            # update the stored variable vector instance
            self.xpred.setvector(xpred)
            self.Ppred = self.F @ Pest @ self.F.T + self.Qm
            self.predictionflag = True
            # print("Predict step", self.xpred.vector, self.Ppred)
        return predict_valid

    def ukf_gen_samples(self):

        # Model defines the order
        order = self.Model.order
        gamma = self.ukf_gamma
        # dfine number of samples and set up arrays to return them
        nsamples = 2*order + 1
        xsample = np.zeros((nsamples, order, 1))
        weight = np.zeros((nsamples, 1, 1))
        # evaluate sigmai from P matrix, note P must be positive definite
        P = self.Ppred*(order + gamma)
        # note uses scipy cholesky which defaults to a upper triangular
        try:
            sigma = cholesky(P, lower=True)
            posdef = True
        except LinAlgError:
            print("Not Positive Definite")
            posdef = False
        # can only evaluate the samples if sigma was evaluated from posdef P
        if posdef is True:
            # sample 0 is the input xhat mean, with weight as defined
            xsample[0] = self.xpred.vector
            weight[0] = gamma/(order + gamma)
            for i in range(order):
                # extract each column from sigma
                sigmai = sigma[:, i:i+1]
                # two samples for each column of sigma +/-
                xsample[2*i+1] = xsample[0] + sigmai
                xsample[2*i+2] = xsample[0] - sigmai
                # corresponding weights are the same
                weight[2*i+1] = 1/(2*(order+gamma))
                weight[2*i+2] = 1/(2*(order+gamma))
            # print("generated samples:", nsamples, xsample)
        return posdef, xsample, weight

    def ukf_eval_samples(self):

        # Model defines the order
        order = self.Model.order
        zsize = np.size(self.zk.vector)
        # dfine number of samples and set up arrays to return them
        nsamples = 2*order + 1
        zsample = np.zeros((nsamples, zsize, 1))
        # check that we have a valid set of predictions
        if hasattr(self, 'xpredsamples') is False:
            # there are no samples in the UKF object
            valid = False
            raise ValueError("UKF samples not initialised")
        elif np.all(self.xpredsamples == 0) is True:
            # there is nothing to predict
            valid = False
            raise ValueError("UKF samples are all zero")
        else:
            valid = True
            xsample = clone_vector(self.xpred)
            for i in range(nsamples):
                # unpack the state sample and the stored input vector
                xsample.setvector(self.xpredsamples[i] / self.xnorm)
                # print(" xsample", i, xsample.vector)
                cAl2O3, ACD, deltaI, I0 = ukf_unpack(xsample, self.uk)
                Vcell, Icell, Vvalid = VIcell(cAl2O3, ACD, deltaI,
                                              I0, self.Cell)
                zsamp = ukf_zpack(Vcell, Icell, self.zk)
                zsample[i] = zsamp.vector * self.znorm
                # print(" zsample", i, zsample[i])
        return valid, zsample

    def ukf_calc_zstats(self):

        # Model defines the order
        order = self.Model.order
        # dfine number of samples and set up arrays to return them
        nsamples = 2*order + 1
        xsize = np.shape(self.xpred.vector)[0]
        zsize = np.shape(self.zk.vector)[0]
        zbar = np.zeros((zsize, 1))
        Pzz = np.zeros((zsize, zsize))
        Pxz = np.zeros((xsize, zsize))
        # check that we have a valid set of predictions
        if hasattr(self, 'zpredsamples') is False:
            # there are no zsamples in the UKF object
            valid = False
            raise ValueError("UKF zsamples not initialised")
        elif np.all(self.zpredsamples == 0) is True:
            # there is nothing to predict
            valid = False
            raise ValueError("UKF zsamples are all zero")
        else:
            valid = True
            for i in range(nsamples):
                # Need zbar evaluated before we can calc Pzz and Pxz
                zbar = zbar + self.weights[i] * self.zpredsamples[i]
            for j in range(nsamples):
                # need raw delta for both z and x
                zdif = self.zpredsamples[j] - zbar
                xdif = self.xpredsamples[j] - self.xpredsamples[0]
                # now calculate the variances
                zz = zdif @ zdif.T
                xz = xdif @ zdif.T
                Pzz = Pzz + self.weights[j] * zz
                Pxz = Pxz + self.weights[j] * xz
        return valid, zbar, Pzz, Pxz

    def update_state(self, z: VariableVector):

        # the UKF must have been initialised, if not then no prediction
        # update cant be done if a prediction step hasn't been made
        # Note that the timestep was set with the previous prediction.
        update_valid = self.initialised

        # step 1 - confirm that the UKF is initialised with a valid
        # prediction available and a valid measurement vector.
        if update_valid is True:
            # update cant be done if a valid prediction step hasn't been made
            predict_valid = self.predictionflag
            # The variable vector definition must match that used to set model
            meas_dict_valid = (z.dictionary == self.Model.z0.dictionary)
            # prediction is only valid if all conditions met
            update_valid = (update_valid and predict_valid and meas_dict_valid)

        # Step 2 - generate state predictions for unscented variance calcs
        if update_valid is True:
            # store the latest measurement vector internally
            self.zk.setvector(z.vector)
            # need to generate the UKF samples first
            valid, xpredsamp, weight = self.ukf_gen_samples()
            if valid is False:
                update_valid = False
                # prediction failed, so zero the stored samples if they exist
                if hasattr(self, 'xpredsamples') is True:
                    self.xpredsamples *= 0
                raise ValueError("UKF samples unable to be calculated")

        # Step 3 - calculate expected measurements based on the state
        # predictions and inputs for unscented variance calcs
        if update_valid is True:
            # store the weiughts and samples internally
            self.weights = weight.copy()
            self.xpredsamples = xpredsamp.copy()
            # Evaluate the predicted measurements for each of the samples
            valid, zpredsamp = self.ukf_eval_samples()
            if valid is False:
                update_valid = False
                # prediction failed, so zero the stored samples if exist
                if hasattr(self, 'zpredsamples') is True:
                    self.zpredsamples *= 0
                raise ValueError("UKF samples unable to be summarised")

        # Step 4 - calculate predicted measurement statistics for the
        # unscented variance calcs
        if update_valid is True:
            # store results as UKF attribute
            self.zpredsamples = zpredsamp.copy()
            valid, zbar, Pzz, Pxz = self.ukf_calc_zstats()
            if valid is False:
                update_valid = False
                # prediction failed, so zero the stored samples if exist
                if hasattr(self, 'Pzz') is True:
                    self.Pzz *= 0
                if hasattr(self, 'Pxz') is True:
                    self.Pxz *= 0
                raise ValueError("UKF zstats unable to be calculated")

        # Step 5 - innovation and variance calculations for the UKF
        if update_valid is True:
            # store results as UKF attribute
            self.Pzz = Pzz
            self.Pxz = Pxz
            self.zbar = zbar
            # Innovation and Innovation Covariance
            innov = (z.vector * self.znorm) - zbar
            # print("zbar, innovation: ", zbar, innov)
            S = Pzz + self.Rm
            # Kgain calculation
            SI = myminv(S)
            Kgain = Pxz @ SI
            # and finally update estimates
            xest = self.xpred.vector + Kgain @ innov
            Pest = self.Ppred - Kgain @ S @ Kgain.T
            # include a step to make sure Pest is symmetrical
            # preserving the lower triangle (consistent with cholesky)
            Pest_sym = np.tril(Pest, 0) + np.tril(Pest, -1).T
            det = np.linalg.det(Pest_sym)
            cond = np.linalg.cond(Pest_sym)
            # store the results
            self.xest.setvector(xest)
            self.Pest = Pest_sym
            self.innov = innov
            self.innov_cov = S
            self.Kgain = Kgain
            self.determinant = det
            self.condition_num = cond
            self.predictionflag = False
        return update_valid


def ukf_execute(Filter: UKF, zt: np.array, ut: np.array,
                dT: float, plotresult=True, xtrue=None):

    nzt = np.shape(zt)[0]
    nsteps = np.shape(ut)[0]
    nobs = np.shape(zt)[1]
    ninputs = np.shape(ut)[1]
    zsize = np.shape(Filter.zk.vector)[0]
    usize = np.shape(Filter.uk.vector)[0]
    xsize = Filter.Model.order
    # to run the filter we need to have te same number of observations in
    # zt and ut, and the number of variables needs to match the initialised
    # filter
    if ((nzt == nsteps) and (Filter.initialised is True) and
        (nobs == zsize) and (ninputs == usize)):
        # initialise the result arrays to correct size
        xest = np.zeros((nsteps, xsize, 1))
        xpred = np.zeros((nsteps, xsize, 1))
        # Pzz = np.zeros((nsteps, zsize, zsize))
        innov = np.zeros((nsteps, zsize, 1))
        innov_cov = np.zeros((nsteps, zsize, zsize))
        Pest = np.zeros((nsteps, xsize, xsize))
        Ppred = np.zeros((nsteps, xsize, xsize))
        # an array to store condition and determinant
        Condition = np.zeros((nsteps, 2, 1))
        # capture the initialised values
        xest[0] = Filter.xest.vector / Filter.xnorm
        xpred[0] = Filter.xpred.vector / Filter.xnorm
        # innov[0] = Filter.innov
        innov_cov[0] = Filter.Rm
        Pest[0] = Filter.Pest
        Ppred[0] = Filter.Ppred
        InputVector = clone_vector(Filter.uk, zeros=True)
        ObsVector = clone_vector(Filter.zk, zeros=True)
        for i in range((nsteps-1)):
            valid_pred = False
            valid_update = False
            # prediction step
            InputVector.setvector(ut[i])
            valid_pred = Filter.predict_state(InputVector, dT)
            if valid_pred is True:
                # Update step
                ObsVector.setvector(zt[i+1])
                valid_update = Filter.update_state(ObsVector)
            else:
                print("Prediction step invalid at input", i)
            if valid_update is True:
                # capture the UKF states
                xest[i+1] = Filter.xest.vector / Filter.xnorm
                xpred[i+1] = Filter.xpred.vector / Filter.xnorm
                # Pzz[i+1] = Filter.Pzz
                innov[i+1] = Filter.innov
                innov_cov[i+1] = Filter.innov_cov
                Pest[i+1] = Filter.Pest
                Ppred[i+1] = Filter.Ppred
                Condition[i+1] = np.array([[Filter.condition_num],
                                          [Filter.determinant]])
            else:
                print("Update step invalid at observation", i+1)
            if (valid_pred is False) or (valid_update is False):
                # Stop executing the UKF at that point
                # The state of Filter is preserved and results to that
                # point returned
                if i == 0:
                    plotresult = False
                break
            # remainder, rqst = ma.modf(i/60)
            # if remainder == 0:
            #    tplot = i*dT
            #    print(f'Determinant at {tplot} sec: {Filter.determinant}')
            #    print(f'Cond Number at {tplot} sec: {Filter.condition_num}')
        if plotresult is True:
            zlist = list(ObsVector.dictionary.keys())
            xlist = list(Filter.xpred.dictionary.keys())
            ukf_plot_results(zt, zlist, innov, innov_cov,
                             xest, xpred, xlist, dT, Filter.znorm, xtrue)
            ukf_plot_condition(Condition, dT)
    return (xest, Pest, xpred, Ppred, innov, innov_cov, Condition)


# plot results
def ukf_plot_results(z, zlist, innov, innov_cov, xest, xpred,
                     xlist, dt, znorm, xtrue=None):

    tsteps, xcount, dontcare = np.shape(xest)
    zcount = np.shape(z)[1]
    t = np.linspace(0, dt*tsteps, tsteps)
    for xi in range(xcount):
        fig, axx = plt.subplots(num=xi, figsize=(8, 6), dpi=128, clear=True)
        if xtrue is not None:
            axx.plot(t, xtrue[:, xi, 0], label='Simulated')
        axx.plot(t, xest[:, xi, 0], label='Estimated')
        axx.plot(t, xpred[:, xi, 0], label='Predicted')
        axx.set_xlabel('time (sec)')
        axx.set_ylabel(xlist[xi])
        axx.set_title(f'UKF State Variable Outputs: {xlist[xi]}')
        axx.legend()
    # calculate the innovation standard deviation for ea
    # normalised variance term is scalar by definition
    S = np.sqrt(innov_cov)
    for si in range(zcount):
        Splus = S[:, si, si]
        Sminus = -S[:, si, si]
        fig, axs = plt.subplots(num=(si+xcount), figsize=(8, 6),
                                dpi=128, clear=True)
        axs.plot(t, innov[:, si, 0], label='Innovation')
        axs.plot(t, Splus, label='sqrt(Diag S)')
        axs.plot(t, Sminus, label='-sqrt(Diag S)')
        axs.set_xlabel('time (sec)')
        axs.set_ylabel(f'Innovation {zlist[si]}')
        axs.set_title(f'UKF Innovation for Measurement: {zlist[si]}')
        axs.legend()
    # calculate the innovation relative to z
    innov_cov_norm = np.zeros((tsteps, 1, 1), dtype=float)
    innov_cov_norm_avg = np.zeros((tsteps, 1, 1), dtype=float)
    for ni in range(tsteps):
        SI = myminv(innov_cov[ni])
        innov_cov_norm[ni] = innov[ni].T @ SI @ innov[ni]
        if ni == 0:
            innov_cov_norm_avg[ni] = innov_cov_norm[ni]
        else:
            window = min(120, ni+1)
            # innov_cov_norm_avg[ni] = (innov_cov_norm_avg[ni-1]*ni +
            #                          innov_cov_norm[ni])/(ni+1)
            windowsum = np.sum(innov_cov_norm[(ni-window+1):(ni+1)])
            innov_cov_norm_avg[ni] = windowsum/window
    innov = z - (innov / znorm)
    for zi in range(zcount):
        fig, axz = plt.subplots(num=(zi+xcount+zcount), figsize=(8, 6),
                                dpi=128, clear=True)
        axz.plot(t, innov[:, zi, 0], label='Prediction')
        axz.plot(t, z[:, zi, 0], label='Observation')
        axz.set_xlabel('time (sec)')
        axz.set_ylabel(zlist[zi])
        axz.set_title(f'UKF Prediction for Measurement: {zlist[zi]}')
        axz.legend()
    # normalised variance term is scalar by definition
    fig, axnc = plt.subplots(num=(zcount+xcount+zcount), figsize=(8, 6),
                             dpi=128, clear=True)
    axnc.axhline(y=2.60, color='g')
    axnc.plot(t, innov_cov_norm[:, 0, 0], label='Normalised Innov Covariance')
    axnc.plot(t, innov_cov_norm_avg[:, 0, 0], label='Running Average')
    axnc.axhline(y=1.48, color='g')
    axnc.set_xlabel('time (sec)')
    axnc.set_ylabel('normalised covariance')
    axnc.set_yscale('log')
    # note that these lines are for zdim 2 chi-sq test
    axnc.set_title("UKF Innovation Covariance Normalised")
    axnc.legend()


def ukf_plot_condition(Cond, dt):

    tsteps = np.shape(Cond)[0]
    t = np.linspace(0, dt*tsteps, tsteps)
    fig, axcon = plt.subplots(num=100, figsize=(8, 6),
                              dpi=128, clear=True)
    axcon.plot(t, Cond[:, 0, 0], label='Condition Number')
    axcon.set_xlabel('time (sec)')
    axcon.set_ylabel('Condition Number')
    axcon.set_title("UKF Condition Number on Pest")
    # fig, axcon = plt.subplots(num=101, figsize=(8, 6),
    #                           dpi=128, clear=True)
    # axcon.plot(t, Cond[:, 1, 0], label='Determinant')
    # axcon.set_xlabel('time (sec)')
    # axcon.set_ylabel('Determinant')
    # note that these lines are for zdim 2 chi-sq test
    # axcon.set_title("UKF Determinant on Pest")
