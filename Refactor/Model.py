# -*- coding: utf-8 -*-
"""
Created on Tue 29 Jun 2021 08:00

@author: Mark.Illingworth
"""

import numpy as np
from copy import deepcopy


class VariableVector:

    # A class that defines a vector in a KF system
    # vector can be used for state, input or measurements generically

    def __init__(self, **keywords: float):
        vectordict = {}
        vector = np.array([], dtype=float)
        index = -1
        for key in keywords:
            index += 1
            vectordict[key] = index
            appended = np.array([keywords[key]])
            vector = np.append(vector, appended)
        vector = vector[:, np.newaxis]
        self._vector = vector
        self._dictionary = vectordict

    @property
    def vector(self):

        # public function to return the vector
        return self._vector

    @property
    def dictionary(self):

        # public function to return the dictionary of the vector
        return self._dictionary

    def variable(self, itemstr: str):

        # public function to return the named item according to dictionary
        # returnvar = self.vector[self.dictionary[itemstr]].item()
        return self.vector[self.dictionary[itemstr]].item()

    def setvariable(self, itemstr: str, newValue: float):

        # public function to set an individual variable by name
        updatedflag = True
        if itemstr in self._dictionary:
            indx = self.dictionary[itemstr]
            self._vector[indx] = np.array([newValue])
        else:
            updatedflag = False
        return updatedflag

    def setvector(self, newVect: np.array):

        # public function to set the whole vector
        # checks for the size being the same but copes with either a row
        # or a column.
        size = np.size(newVect)
        shape = np.shape(newVect)
        updatedflag = True
        if (np.size(self.vector) == size) and (shape.count(size) >= 1):
            # Note that a 1x1 vector will have two axes with size 1
            np.copyto(self._vector, newVect.reshape((size, 1)), casting='safe')
        else:
            # raise ValueError("Incorrect Vector Size")
            updatedflag = False
        return updatedflag


class ControlMatrix:

    # A class that defines a 2D matrix in a KF linear discrete system
    # with dependency on a timestep delta T
    # order 0 is no dependency on delta T
    # order 1 is linear dependency on delta T
    # order 2 is quadratic etc

    def __init__(self, rows, columns, order=1):

        self._array = np.zeros((order+1, rows, columns), float)

    @property
    def array_display(self):

        # public function to return the vector
        return self._array.astype(float)

    def array_eval(self, dT=0):

        result = np.zeros_like(self._array[0], dtype=float)
        dims = np.shape(self._array)
        for index in range(dims[0]):
            result = result + self._array[index]*(dT**index)
        return result

    def set_array_layer(self, index: int, inputarray: np.array):

        # note that this accepts any array but will only update
        # if the shape is correct
        inputshape = np.shape(inputarray)
        updatedflag = True
        if (index < 0) or (index >= self.array_display.shape[0]):
            print("Wrong Order")
            updatedflag = False
        else:
            selfshape = np.shape(self._array[index])
            if (inputshape == selfshape):
                self._array[index] = inputarray.astype(float)
            else:
                print("Wrong Shape")
                updatedflag = False
        return updatedflag

    def set_array(self, inputarray: np.array):

        # note that this accepts any array but will only update
        # if the shape is correct with 3 axes
        inputshape = inputarray.shape
        updatedflag = True
        selfshape = self.array_display.shape
        if (inputshape == selfshape):
            self._array = inputarray.astype(float)
        else:
            print("Wrong Shape")
            updatedflag = False
        return updatedflag


def clone_vector(varvector: VariableVector, *, zeros=False, newvector=None):

    # zeros parameter has highest priority, but if its false
    # then newvector parameter is checked
    newvarvector = deepcopy(varvector)
    cloneshape = np.shape(newvarvector.vector)
    if zeros is True:
        newshape = cloneshape
        newvector = np.zeros(newshape)
    elif (newvector is None):
        newshape = None
    else:
        newshape = np.shape(newvector)
    if cloneshape == newshape:
        newvarvector.setvector(newvector)
    return newvarvector
