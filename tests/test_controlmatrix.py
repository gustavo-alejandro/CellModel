# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 13:31:25 2021

@author: Mark.Illingworth
"""

import numpy as np
import pytest

from Refactor.Model import ControlMatrix
from Refactor.Model import VariableVector
from Refactor.UKFClass import DiscreteModel


testlist_initmatrix = [
    (0, 2, 2, np.arange(4).reshape((1, 2, 2)),
     True),        # zero order 2 x 2
    (1, 3, 2, np.arange(12),
     False),       # first order 3 x 2, but supplied matrix is wrong shape
    (2, 4, 4, np.arange(48).reshape((3, 4, 4)),
     True),        # second order 4 x 4
    (1, 3, 3, np.arange(18).reshape((2, 3, 3)),
     True),        # first order 3 x 3
    (0, 4, 3, np.arange(16).reshape((1, 4, 4)),
     False),       # zero order 4 x 3, but supplied matrix is wrong shape
    ]


@pytest.mark.parametrize("order," +
                         "rows," +
                         "columns," +
                         "testarray," +
                         "expectedresult", testlist_initmatrix)
def test_controlmatrix(order, rows, columns, testarray, expectedresult):

    CM = ControlMatrix(rows, columns, order)
    expected_shape = (order+1, rows, columns)
    expected_init = np.zeros(expected_shape, float)
    # confirm the created array is equal using NumPy array compare
    assert np.array_equal(CM.array_display, expected_init)
    assert CM.array_display.dtype == float
    # now attempt to set the array to an arbitrary value
    validresult = CM.set_array(testarray)
    assert validresult is expectedresult
    if validresult is True:
        # confirm that the matrix is updated
        assert np.array_equal(CM.array_display, testarray)
        assert CM.array_display.dtype == float
    else:
        # confirm that the matrix is still as initialised
        assert np.array_equal(CM.array_display, expected_init)


testlist_initlayer = [
    (0, 2, 2, 1, np.arange(4).reshape((2, 2)),
     False),        # zero order 2 x 2, but trying to change layer 1
    (1, 3, 2, 1, np.arange(6).reshape((3, 2)),
     True),         # first order 3 x 2, setting layer 1
    (2, 4, 4, 0, np.arange(12).reshape((3, 4)),
     False),        # second order 4 x 4, trying to change layer 0 wrong shape
    (1, 3, 3, 1, np.arange(18).reshape((2, 3, 3)),
     False),        # first order 3 x 3, trying to change layer but whole array
    (0, 4, 3, 0, np.arange(12).reshape((4, 3)),
     True),       # zero order 4 x 3, changing layer 0
    ]


@pytest.mark.parametrize("order," +
                         "rows," +
                         "columns," +
                         "layer," +
                         "testarray," +
                         "expectedresult", testlist_initlayer)
def test_controlmatrixlayer(order, rows, columns,
                            layer, testarray, expectedresult):

    CM = ControlMatrix(rows, columns, order)
    expected_shape = (order+1, rows, columns)
    expected_init = np.zeros(expected_shape, float)
    # confirm the created array is equal using NumPy array compare
    assert np.array_equal(CM.array_display, expected_init)
    assert CM.array_display.dtype == float
    # now attempt to set a layer to an arbitrary value
    validresult = CM.set_array_layer(layer, testarray)
    assert validresult is expectedresult
    if validresult is True:
        # confirm that the matrix is updated
        assert np.array_equal(CM.array_display[layer], testarray)
        assert CM.array_display.dtype == float
    else:
        # confirm that the matrix is still as initialised
        assert np.array_equal(CM.array_display, expected_init)


testlist_evaldT = [
    (np.arange(4).reshape((1, 2, 2)), 1,
     np.arange(4).reshape((2, 2))),  # zero order 2 x 2
    (np.arange(4).reshape((1, 2, 2)), 0,
     np.arange(4).reshape((2, 2))),  # zero order 2 x 2 zero dT
    (np.arange(4).reshape((1, 2, 2)), 3456.7899,
     np.arange(4).reshape((2, 2))),  # zero order 2 x 2 big float dT
    (np.arange(12).reshape((2, 3, 2)), 5,
     np.arange(30, 66, 6, dtype=float).reshape((3, 2))),  # first order 3 x 2
    (np.arange(48).reshape((3, 4, 4)), 1,
     np.arange(48, 96, 3, dtype=float).reshape((4, 4))),
    (np.arange(48).reshape((3, 4, 4)), 0,
     np.arange(16).reshape((4, 4))),
    (np.array([[[0, 1], [1, 0]], [[1, 0], [0, 1]]]), 10,
     np.array([[10, 1], [1, 10]])),
    (np.array([[[0, 3.2], [2.3, 0]], [[1, 0], [0, 1]]]), 3.14159,
     np.array([[3.14159, 3.2], [2.3, 3.14159]])),
    ]


@pytest.mark.parametrize("testarray," +
                         "dT," +
                         "expectedresult", testlist_evaldT)
def test_controleval(testarray, dT, expectedresult):
    order, rows, columns = np.shape(testarray)
    # zero order control matrix wrt dT still has a dimension oof 1 in that axis
    order -= 1
    CM = ControlMatrix(rows, columns, order)
    validresult = CM.set_array(testarray)
    # another test shouyld have failed if this part also does
    assert validresult is True
    assert np.array_equal(CM.array_display, testarray)
    assert CM.array_display.dtype == float
    resultdT = CM.array_eval(dT)
    # should get a result that is correct, float and 2D with shape (R, C)
    assert np.array_equal(resultdT, expectedresult)
    assert resultdT.dtype == float
    assert resultdT.shape == (rows, columns)


testlist_DiscreteModelTrue = [
    ((3, 3, 1), (3, 3, 1), (3, 3, 1), (3, 3, 0), (1, 1, 0),
     1, 3, 3, True),        # 1 measure, 3 states, 3 inputs
    ((4, 4, 1), (4, 3, 1), (4, 4, 2), (4, 4, 0), (2, 2, 0),
     2, 3, 4, True),        # 2 measure, 4 states, 3 inputs
    ((2, 2, 1), (2, 4, 1), (2, 2, 2), (2, 2, 0), (1, 1, 0),
     1, 4, 2, True),        # 1 measure, 2 states, 4 inputs
    ]


@pytest.mark.parametrize("Fm," +
                         "Bm," +
                         "Gm," +
                         "Qm," +
                         "Rm," +
                         "Zv," +
                         "Uv," +
                         "Xv," +
                         "expectedresult", testlist_DiscreteModelTrue)
def test_DiscreteModel(Fm, Bm, Gm, Qm, Rm, Zv, Uv, Xv, expectedresult):

    v1 = VariableVector(var1=0.0)
    v2 = VariableVector(var1=0.0, var2=0.0)
    v3 = VariableVector(var1=0.0, var2=0.0, var3=0.0)
    v4 = VariableVector(var1=0.0, var2=0.0, var3=0.0, var4=0.0)
    VecList = (v1, v2, v3, v4)
    F = ControlMatrix(Fm[0], Fm[1], Fm[2])
    B = ControlMatrix(Bm[0], Bm[1], Bm[2])
    G = ControlMatrix(Gm[0], Gm[1], Gm[2])
    Q = ControlMatrix(Qm[0], Qm[1], Qm[2])
    R = ControlMatrix(Rm[0], Rm[1], Rm[2])
    z0 = VecList[(Zv-1)]
    u0 = VecList[(Uv-1)]
    x0 = VecList[(Xv-1)]
    DM = DiscreteModel(F, B, G, Q, R, z0, u0, x0)
    assert DM.defined is expectedresult
    assert np.shape(DM.F(0)) == (Fm[0], Fm[1])
    assert np.shape(DM.B(0)) == (Bm[0], Bm[1])
    assert np.shape(DM.G(0)) == (Gm[0], Gm[1])
    assert np.shape(DM.Q(0)) == (Qm[0], Qm[1])
    assert np.shape(DM.R(0)) == (Rm[0], Rm[1])
    assert np.shape(DM.z0.vector) == (Zv, 1)
    assert np.shape(DM.u0.vector) == (Uv, 1)
    assert np.shape(DM.x0.vector) == (Xv, 1)


testlist_DiscreteModelFalse = [
    ((3, 2, 1), (3, 3, 1), (3, 3, 1), (3, 3, 0), (1, 1, 0),
     1, 3, 3, "F array not square or inconsistent with x"),    # F Not square
    ((3, 3, 1), (3, 3, 1), (3, 3, 1), (3, 3, 0), (1, 1, 0),
     1, 3, 4, "F array not square or inconsistent with x"),    # x size wrong
    ((3, 3, 1), (4, 3, 1), (3, 3, 1), (3, 3, 0), (1, 1, 0),
     1, 3, 3, "B array inconsistent with x or u"),    # B size wrong cf x
    ((3, 3, 1), (3, 3, 1), (3, 3, 1), (3, 3, 0), (1, 1, 0),
     1, 4, 3, "B array inconsistent with x or u"),    # u size wrong
    ((3, 3, 1), (3, 3, 1), (3, 3, 1), (3, 2, 0), (1, 1, 0),
     1, 3, 3, "Q array not square or inconsistent with x"),    # Q Not square
    ((3, 3, 1), (3, 3, 1), (3, 3, 1), (4, 4, 0), (1, 1, 0),
     1, 3, 3, "Q array not square or inconsistent with x"),    # Q too big
    ((3, 3, 1), (3, 3, 1), (3, 3, 1), (3, 3, 0), (2, 1, 0),
     1, 3, 3, "R array not square or inconsistent with z"),    # R Not square
    ((3, 3, 1), (3, 3, 1), (3, 3, 1), (3, 3, 0), (2, 2, 0),
     1, 3, 3, "R array not square or inconsistent with z"),    # R too big
    ((3, 3, 1), (3, 3, 1), (4, 3, 1), (3, 3, 0), (1, 1, 0),
     1, 3, 3, "G array inconsistent with x or Q"),     # G shape wrong cf x
    ((3, 3, 1), (3, 3, 1), (3, 4, 1), (3, 3, 0), (1, 1, 0),
     1, 3, 3, "G array inconsistent with x or Q"),     # G shape wrong cf Q
    ((3, 3, 1), (3, 3, 1), (3, 3, 1), (3, 3, 0), (1, 1, 0),
     1, 3, 3, "Q or R array contain negative variances"),  # negative variance
    ]


@pytest.mark.parametrize("Fm," +
                         "Bm," +
                         "Gm," +
                         "Qm," +
                         "Rm," +
                         "Zv," +
                         "Uv," +
                         "Xv," +
                         "expectedresult", testlist_DiscreteModelFalse)
def test_DiscreteModelFail(Fm, Bm, Gm, Qm, Rm, Zv, Uv, Xv, expectedresult):

    v1 = VariableVector(var1=0.0)
    v2 = VariableVector(var1=0.0, var2=0.0)
    v3 = VariableVector(var1=0.0, var2=0.0, var3=0.0)
    v4 = VariableVector(var1=0.0, var2=0.0, var3=0.0, var4=0.0)
    VecList = (v1, v2, v3, v4)
    F = ControlMatrix(Fm[0], Fm[1], Fm[2])
    B = ControlMatrix(Bm[0], Bm[1], Bm[2])
    G = ControlMatrix(Gm[0], Gm[1], Gm[2])
    Q = ControlMatrix(Qm[0], Qm[1], Qm[2])
    R = ControlMatrix(Rm[0], Rm[1], Rm[2])
    if expectedresult == "Q or R array contain negative variances":
        Q.set_array_layer(0, -1 * np.ones((Qm[0], Qm[1])))
    z0 = VecList[(Zv-1)]
    u0 = VecList[(Uv-1)]
    x0 = VecList[(Xv-1)]
    with pytest.raises(ValueError) as error:
        DM = DiscreteModel(F, B, G, Q, R, z0, u0, x0)
    assert str(error.value) == expectedresult
