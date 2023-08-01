# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 13:31:25 2021

@author: Mark.Illingworth
"""

import numpy as np
import pytest

from Refactor.Model import VariableVector

testlist_1var = [
    (4.532, 1),    # one float
    (-4.268, 1),   # one negative
    (42, 1),       # one integer
    ]


@pytest.mark.parametrize("var1," +
                         "expected_size", testlist_1var)
def test_1var_vector_create(var1, expected_size):

    testvector = VariableVector(variable1=var1)
    expected_vector = np.array([[var1]])
    expected_dictionary = {'variable1': 0}
    # confirm the vector is equal using NumPy array compare
    assert np.array_equal(testvector.vector, expected_vector)
    assert np.shape(testvector.vector) == (expected_size, 1)
    assert testvector.dictionary == expected_dictionary
    # confirm the variable elements in the vector are correctly set
    assert testvector.variable('variable1') == var1
    # confirm the variable elements in the vector are type float
    assert isinstance(testvector.variable('variable1'), float)


testlist_2var = [
    (4.532, 126000, 2),  # one float, one integer
    (4.268, 3.052, 2),   # two floats
    (42, 65, 2),         # two integers
    ]


@pytest.mark.parametrize("var1," +
                         "var2," +
                         "expected_size", testlist_2var)
def test_2var_vector_create(var1, var2, expected_size):

    testvector = VariableVector(variable1=var1, variable2=var2)
    expected_vector = np.array([[var1], [var2]])
    expected_dictionary = {'variable1': 0, 'variable2': 1}
    # confirm the vector is equal using NumPy array compare
    assert np.array_equal(testvector.vector, expected_vector)
    assert np.shape(testvector.vector) == (expected_size, 1)
    assert testvector.dictionary == expected_dictionary
    # confirm the variable elements in the vector are correctly set
    assert testvector.variable('variable1') == var1
    assert testvector.variable('variable2') == var2
    # confirm the variable elements in the vector are type float
    assert isinstance(testvector.variable('variable1'), float)
    assert isinstance(testvector.variable('variable2'), float)


testlist_3var = [
    (4.532, 126000, -1.7, 3),  # one float, one integer, on negative
    (4.268, 3.052, -6, 3),     # two floats, one negative int
    (42, 65, 1.6, 3),          # two integers, one float
    (-42, -65, -1.6, 3),       # two integers, one float, all negative
    ]


@pytest.mark.parametrize("var1," +
                         "var2," +
                         "var3," +
                         "expected_size", testlist_3var)
def test_3var_vector_create(var1, var2, var3, expected_size):

    testvector = VariableVector(variable1=var1, variable2=var2,
                                variable3=var3)
    expected_vector = np.array([[var1], [var2], [var3]])
    expected_dictionary = {'variable1': 0, 'variable2': 1, 'variable3': 2}
    # confirm the vector is equal using NumPy array compare
    assert np.array_equal(testvector.vector, expected_vector)
    assert np.shape(testvector.vector) == (expected_size, 1)
    assert testvector.dictionary == expected_dictionary
    # confirm the variable elements in the vector are correctly set
    assert testvector.variable('variable1') == var1
    assert testvector.variable('variable2') == var2
    assert testvector.variable('variable3') == var3
    # confirm the variable elements in the vector are type float
    assert isinstance(testvector.variable('variable1'), float)
    assert isinstance(testvector.variable('variable2'), float)
    assert isinstance(testvector.variable('variable3'), float)


testlist_vector_update = [
    (4.532, 126000, np.array([[4.1], [126500]]), True),  # single column
    (4.268, 3.052, np.array([4.1, 126500]), True),     # single row
    (4.268, 3.052, np.array([4.1, 126500, 1.1]), False),     # wrong size
    (4.268, 3.052, 4.126500, False),     # scalar instead of array
    ]


@pytest.mark.parametrize("var1," +
                         "var2," +
                         "update," +
                         "expectedresult", testlist_vector_update)
def test_vector_update(var1, var2, update, expectedresult):

    testvector = VariableVector(variable1=var1, variable2=var2)
    expected_vector = np.array([[var1], [var2]])
    expected_dictionary = {'variable1': 0, 'variable2': 1}
    # confirm the vector is equal using NumPy array compare
    assert np.array_equal(testvector.vector, expected_vector)
    assert testvector.dictionary == expected_dictionary
    # confirm the variable elements in the vector are correctly set
    assert testvector.variable('variable1') == var1
    assert testvector.variable('variable2') == var2
    # perform a vector update and test new values applied
    updateflag = testvector.setvector(update)
    assert updateflag == expectedresult
    if expectedresult is True:
        # confirm the variable elements in the vector are correctly set
        assert testvector.variable('variable1') == update.item(0)
        assert testvector.variable('variable2') == update.item(1)
    else:
        assert testvector.variable('variable1') == var1
        assert testvector.variable('variable2') == var2
    # confirm the variable elements in the vector are type float
    assert isinstance(testvector.variable('variable1'), float)
    assert isinstance(testvector.variable('variable2'), float)


testlist_variable_update = [
    (4.532, 126000, -1.7, 'variable1', 3, True),   # update the first variable
    (4.268, 3.052, -6, 'variable2', 1.678, True),  # update the second variable
    (42, 65, 1.6, 'variable3', 0, True),           # update the third variable
    (42, 65, 1.6, 'variable4', 5.8, False),        # update the third variable
    ]


@pytest.mark.parametrize("var1," +
                         "var2," +
                         "var3," +
                         "index," +
                         "newval," +
                         "expectedresult", testlist_variable_update)
def test_variable_update(var1, var2, var3, index, newval, expectedresult):

    testvector = VariableVector(variable1=var1, variable2=var2,
                                variable3=var3)
    expected_vector = np.array([[var1], [var2], [var3]])
    expected_dictionary = {'variable1': 0, 'variable2': 1, 'variable3': 2}
    # confirm the vector is equal using NumPy array compare
    assert np.array_equal(testvector.vector, expected_vector)
    assert testvector.dictionary == expected_dictionary
    # confirm the variable elements in the vector are correctly set
    assert testvector.variable('variable1') == var1
    assert testvector.variable('variable2') == var2
    assert testvector.variable('variable3') == var3
    # perform a variable update and test new value applied
    updateflag = testvector.setvariable(index, newval)
    assert updateflag == expectedresult
    if expectedresult is True:
        # confirm the new variable is updated
        assert testvector.variable(index) == newval
    else:
        # confirm the variable elements in the vector are unchanged
        assert testvector.variable('variable1') == var1
        assert testvector.variable('variable2') == var2
        assert testvector.variable('variable3') == var3
    # confirm the variable elements in the vector are type float
    assert isinstance(testvector.variable('variable1'), float)
    assert isinstance(testvector.variable('variable2'), float)
    assert isinstance(testvector.variable('variable3'), float)
