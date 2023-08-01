# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:45:34 2021

@author: Mark.Illingworth
"""

# Do imports first
import Global as Glob
import CellVoltage as Meas
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
from scipy.linalg import cholesky


# a fix because python does not do edge cases
def myminv(A):
    if np.size(A) == 1:
        return (1/A)
    else:
        return inv(A)


def ukf_samples(xhat, P, gamma, ut):

    # infer shape of state vector - order is number of states
    dims = np.shape(xhat)
    order = dims[0]
    # dfine number of samples and set up arrays to return them
    nsamples = 2*order + 1
    xsample = np.zeros((nsamples, order, 1))
    weight = np.zeros((nsamples, 1, 1))
    zsample = np.zeros((nsamples, 1, 1))
    # sample 0 is the input xhat mean, with weight as defined
    xsample[0] = xhat
    weight[0] = gamma/(order + gamma)
    # evaluate sigmai from P matrix, note P must be positive definite
    P = P*(order + gamma)
    sigma = cholesky(P, lower=True)
    # note uses scipy cholesky so defaults to a upper triangular
    # need to verify if it matters upper v lower
    # print(sigma)
    for i in range(order):
        # extract each column from sigma
        sigmai = sigma[:, i:i+1]
        # print(sigmai)
        # two samples for each column of sigma +/-
        xsample[2*i+1] = xhat + sigmai
        xsample[2*i+2] = xhat - sigmai
        # corresponding weights are the same
        weight[2*i+1] = 1/(2*(order+gamma))
        weight[2*i+2] = 1/(2*(order+gamma))
    # unpack the input vector for calculating all zsamples
    for i in range(nsamples):
        zsample[i], valid = Meas.Vcell(xsample[i], ut, Glob.BathConstVector)
    return xsample, weight, zsample


# Main Unscented Kalman Filter
# All in one file
def ukf(F, G, H, B, sigmaq, sigmar, zt, ut, x0, tstop):
    # first set up some space
    dims = np.shape(zt)
    tsteps = dims[0]
    zsize = dims[1]
    dims = np.shape(F)
    xsize = dims[0]
    xest = np.zeros((tsteps, xsize, 1))
    xpred = np.zeros((tsteps, xsize, 1))
    innov = np.zeros((tsteps, zsize, 1))
    Pest = np.zeros((tsteps, xsize, xsize))
    Ppred = np.zeros((tsteps, xsize, xsize))
    S = np.zeros((tsteps, zsize, zsize))
    #
    # some useful matrices - used multiple times
    FT = np.transpose(F)
    HT = np.transpose([H])  # why is this double bracket necessary?
    #
    # now initialise and account for sigmaq being a vector
    dims = np.shape(sigmaq)
    qdim = dims[0]
    Q = np.zeros((qdim, qdim))
    for n in range(qdim):
        # calculate  variance on assumption that all
        # cross correlation terms are zero
        Q[n, n] = sigmaq[n, 0]*sigmaq[n, 0]
    Rm = np.array([[sigmar * sigmar]])
    Qm = G @ Q @ np.transpose(G)
    # Rm = B @ R @ np.transpose(B)
    xest[0] = x0
    Pest[0] = 100*Qm
    xpred[0] = x0
    Ppred[0] = 100*Qm
    if tstop > 0:
        # debug stop and print results after tstop steps
        print("Q:\n", Q)
        print("Qm:\n", Qm)
        print("Initial X:\n", xest[0])
        print("Initial Pest:\n", Pest[0])
        nsteps = tstop
    else:
        nsteps = tsteps-1
    # for i in range(tsteps-1):
    for i in range(nsteps):
        # first do prediction and prediction covariance
        xpred[i+1] = F @ xest[i] + B @ ut[i]
        Ppred[i+1] = F @ Pest[i] @ FT + Qm
        # calculate samples, weights and obs around the prediction
        xsamples, weights, zsamples = ukf_samples(xpred[i+1], Ppred[i+1], 1, ut[i])
        dims = np.shape(zsamples)
        # Note all 3 returned arguments should have same number of samples
        # do we need to explicitly check?
        nsamples = dims[0]
        # print("Number of UKF Samples ", nsamples)
        # compute mean and the innovation
        zbar = 0
        for j in range(nsamples):
            zbar = zbar + weights[j]*zsamples[j]
        innov[i+1] = zt[i+1] - zbar
        # compute variance
        Pzz = np.zeros((zsize, zsize))
        for j in range(nsamples):
            Pzz = Pzz + weights[j]*(zsamples[j] - zbar)*(zsamples[j] - zbar)
        # compute cross correlation
        Pxz = np.zeros((xsize, zsize))
        for j in range(nsamples):
            # need to make the zsamples explicitly an array
            zdif = (zsamples[j] - zbar)
            # zdif = np.array([zdif])
            xdif = (xsamples[j] - xsamples[0])
            xz = xdif @ zdif
            Pxz = Pxz + weights[j] * xz
        # calculate innovation covariance
        S[i+1] = Pzz + Rm
        # S = np.array([S])
        SI = myminv(S[i+1])
        Kgain = Pxz @ SI
        # and finally update estimates
        xest[i+1] = xpred[i+1] + Kgain @ innov[i+1]
        Pest[i+1] = Ppred[i+1] - Kgain @ S[i+1] @ np.transpose(Kgain)
        det = np.linalg.det(Pest[i+1])
        cond = np.linalg.cond(Pest[i+1])
        if ((i+1) == nsteps):
            laststep = True
        else:
            laststep = False
        if ((tstop > 0) and laststep):
            print("==============================================")
            print("Number of Steps:\n", nsteps)
            print("xpred:\n", xpred[i+1])
            print("Ppred:\n", Ppred[i+1])
            print("xsamples:\n", xsamples)
            print("weights:\n", weights)
            print("zsamples:\n", zsamples)
            print("zbar:\n", zbar)
            print("innov:\n", innov[i+1])
            print("Pzz:\n", Pzz)
            print("Pxz:\n", Pxz)
            print("S:\n", S[i+1])
            print("SI:\n", SI)
            print("Kgain:\n", Kgain)
            print("xest:\n", xest[i+1])
            print("Pest:\n", Pest[i+1])
            print("Determinant Pest:\n", det)
            print("Condition Pest:\n", cond)
    # finally return results to be plotted
    return(xest, Pest, xpred, Ppred, innov, S)


# plot results
def plot_ukf(x, z, innov, xest, xpred, dt, S):

    dims = np.shape(x)
    tsteps = dims[0]
    t = np.linspace(0, dt*tsteps, tsteps)
    fig, axconc = plt.subplots()
    axconc.plot(t, x[:, 0, 0], label='Cd True')
    axconc.plot(t, xest[:, 0, 0], label='Cd Est')
    axconc.plot(t, xpred[:, 0, 0], label='Cd Pred')
    axconc.set_xlabel('time')
    axconc.set_ylabel('wt% alumina')
    axconc.set_title("UKF Dissolved Alumina")
    axconc.legend()
    fig, axsludge = plt.subplots()
    axsludge.plot(t, x[:, 1, 0], label='Cun True')
    axsludge.plot(t, xest[:, 1, 0], label='Cun Est')
    axsludge.plot(t, xpred[:, 1, 0], label='Cun Pred')
    axsludge.set_xlabel('time')
    axsludge.set_ylabel('wt% alumina')
    axsludge.set_title("UKF Undissolved Alumina")
    axsludge.legend()
    fig, axacd = plt.subplots()
    axacd.plot(t, x[:, 2, 0], label='D True')
    axacd.plot(t, xest[:, 2, 0], label='D Est')
    axacd.plot(t, xpred[:, 2, 0], label='D Pred')
    axacd.set_xlabel('time')
    axacd.set_ylabel('ACD (cm)')
    axacd.set_title("UKF ACD")
    axacd.legend()
    # calculate the innovation standard deviation
    S = np.sqrt(S)
    Sneg = -S
    fig, axinnov = plt.subplots()
    axinnov.plot(t, innov[:, 0, 0], label='Innovation')
    axinnov.plot(t, S[:, 0], label='sqrt(S)')
    axinnov.plot(t, Sneg[:, 0], label='-sqrt(S)')
    axinnov.set_xlabel('time')
    axinnov.set_ylabel('Differential Voltage (V)')
    axinnov.set_title("UKF Innovation")
    axinnov.legend()
    # calculate the innovation relative to z
    innov = z - innov
    fig, axzcomp = plt.subplots()
    axzcomp.plot(t, innov[:, 0, 0], label='Prediction')
    axzcomp.plot(t, z[:, 0, 0], label='Observation')
    axzcomp.set_xlabel('time')
    axzcomp.set_ylabel('Voltage (V)')
    axzcomp.set_title("UKF Prediction")
    axzcomp.legend()
