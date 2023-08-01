# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 13:34:03 2021

@author: Mark.Illingworth
"""

import numpy as np
from numpy.linalg import inv
from scipy.linalg import cholesky, LinAlgError
import pytest

from Refactor.Model import ControlMatrix
from Refactor.Model import VariableVector
from Refactor.UKFClass import DiscreteModel, UKF, ukf_unpack, ukf_zpack
from Refactor.Properties import ValidProperty, CellProperties
from Refactor.Model import VariableVector, ControlMatrix, clone_vector
from Refactor.CellVoltage import VIcell


def ukf_config_for_testing(dt):

    # define the initial vectors
    x0 = VariableVector(cAl2O3=3.5, uAl2O3=0.5, ACD=2.9193)
    u0 = VariableVector(I=126000.0, g=0.0, B=0.0)
    z0 = VariableVector(Vcell=4.465858)

    # Define a test cell and set the alumina conc
    cp = CellProperties()
    cp.cAl2O3.value = x0.variable('cAl2O3')

    # Define the model, multiple control matrices
    # F matrix - transformation
    F = ControlMatrix(3, 3, 1)
    test = F.set_array_layer(0,
                             np.array([[1.0, 0.0, 0.0],
                                       [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    assert test is True
    test = F.set_array_layer(1,
                             np.array([[0.0, cp.kdiss, 0.0],
                                       [0.0, -cp.kdiss, 0.0],
                                       [0.0, 0.0, 0.0]]))
    assert test is True
    # G matrix - transformation for noise terms
    G = ControlMatrix(3, 3, 1)
    test = G.set_array_layer(0,
                             np.array([[0.0, 0.0, 0.0],
                                       [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]))
    assert test is True
    test = G.set_array_layer(1,
                             np.array([[1.0, 0.0, 0.0],
                                       [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    assert test is True
    # B matrix - transformation for noise terms
    B = ControlMatrix(3, 3, 1)
    test = B.set_array_layer(0,
                             np.array([[0.0, 100*(1-cp.r)/cp.m, 0.0],
                                       [0.0, 100*cp.r/cp.m, 0.0],
                                       [0.0, 0.0, 1.0]]))
    assert test is True
    gamma, validg = cp.gamma
    assert validg is True
    test = B.set_array_layer(1,
                             np.array([[-gamma/cp.m, 0.0, 0.0],
                                       [0.0, 0.0, 0.0],
                                       [cp.alpha, 0.0, 0.0]]))
    assert test is True
    # Q Matrix from the UKF noise standard deviations
    q1dstd = 0.0002    # alumina concentration variation
    q1ustd = 0.0002    # undisolved alumina concentration variation
    q2std = 0.0000001  # acd varation needs to be much lower
    Q = ControlMatrix(3, 3, 0)
    test = Q.set_array_layer(0, np.array([[q1dstd**2, 0.0, 0.0],
                                          [0.0, q1ustd**2, 0.0],
                                          [0.0, 0.0, q2std**2]]))
    assert test is True
    # R Matrix from the UKF noise standard deviations
    rstd = 0.0010      # measurement variation (volts)
    R = ControlMatrix(1, 1, 0)
    test = R.set_array_layer(0, np.array([[rstd**2]]))
    assert test is True

    # OK now to define the complete model
    alumina3state = DiscreteModel(F, B, G, Q, R, z0, u0, x0)
    assert alumina3state.defined is True

    # Finally initialise the valid UKF
    UKF3state = UKF(alumina3state, cp, x0, dt)

    return UKF3state


testlist_ukfmethods = [
    (5, np.array([[126000.0], [0.0], [0.0]]), np.array([[4.47010343]]),
     np.array([[3.49917482], [0.495], [2.91927925]]),
     np.array([[3.4953567], [0.49496258], [2.91927928]]),
     np.array([[0.00226019]]), np.array([[1.02998578e-06]]),
     True),        # 5 sec step
    (10, np.array([[126000.0], [0.0], [0.0]]), np.array([[4.4727236]]),
     np.array([[3.49834964], [0.49], [2.9192585]]),
     np.array([[3.46802227], [0.4894117], [2.91925871]]),
     np.array([[0.00487121]]), np.array([[1.12058762e-06]]),
     True),        # 10 sec step
    ]


@pytest.mark.parametrize("dt," +
                         "uk," +
                         "zkp1," +
                         "xpred," +
                         "xest," +
                         "innov," +
                         "incov," +
                         "result", testlist_ukfmethods)
def test_ukf_methods(dt, uk, zkp1, xpred, xest, innov, incov, result):

    # generate the UKF instance
    TestUKF = ukf_config_for_testing(dt)
    # need the vectors to be in standard format
    ukvect = clone_vector(TestUKF.uk, newvector=uk)
    zkp1vect = clone_vector(TestUKF.zk, newvector=zkp1)

    # confirm it was initialised correctly
    assert TestUKF.initialised is result

    # prediction step using the current input vector
    TestUKF.predict_state(ukvect, dt)
    assert TestUKF.predictionflag is result
    assert np.array_equal(TestUKF.uk.vector, uk)
    print(TestUKF.xpred.vector, xpred)
    assert np.allclose(TestUKF.xpred.vector, xpred)

    # should be able to generate UKF xsamples now
    test, xsamparray, weightsarray = TestUKF.ukf_gen_samples()
    assert test is result
    assert np.allclose(xsamparray[0], xpred)
    assert np.shape(xsamparray)[0] == 7
    assert np.shape(weightsarray)[0] == 7

    # but cant evaluate them yet as they havent been stored in TestUKF
    with pytest.raises(ValueError) as error:
        test, zsamparray = TestUKF.ukf_eval_samples()
    assert str(error.value) == "UKF samples not initialised"

    # now we can do an estimation with the k+1 z vector
    TestUKF.update_state(zkp1vect)
    assert TestUKF.predictionflag is not result
    assert np.array_equal(TestUKF.zk.vector, zkp1)
    assert np.allclose(TestUKF.xest.vector, xest)
  
