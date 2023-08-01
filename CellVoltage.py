# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:38:06 2021

@author: Mark.Illingworth
"""
import Global as Glob
import BathProperties as Bath
import ElectroChemistry as Electro


def Vcell(StateVector, InputVector, ConstVector):
    """
    Vcell = Erev + Eca + Esa + Ecc + I_line*Rpath + I_line*Glob.Rext
    Erev: ReversiblePotential
    Eca: AnodeConcOverVolt
    Esa: AnodeSurfOverVolt
    Ecc: CathodeConcOverVolt
    I_line:
    Rpath: external resistance ohms
    Rext:
    :param StateVector:
    :param InputVector:
    :param ConstVector:
    :return:
    """

    # Non-linear function to define H
    # Note that while Cun and Cd are related in the state equations
    # Cun is not used to calculate Vcell

    # Unpack State, Input and Constant Vectors - need to make this more 
    # general in future
    Cd = StateVector[0]
    Cun = StateVector[1]
    D = StateVector[2]
    I_line = InputVector[0]
    g = InputVector[1]
    BM = InputVector[2]
    C_AlF3 = ConstVector[0]
    C_CaF2 = ConstVector[1]
    C_MgF2 = ConstVector[2]
    C_LiF = ConstVector[3]
    T_bath = ConstVector[4]
    # Valid Flag will be set to false if any of the properties or
    # electrochemical functions returns an invalid argument response
    ValidFlag = True

    # Calculate bath and bubble characteristics currently assuming
    # constant bath chemistry and temperature
    (BR, valid) = Bath.Ratio(Cd, C_AlF3, C_CaF2, C_LiF, C_MgF2)
    if valid is False:
        ValidFlag = False
    (kbath, valid) = Bath.Conductivity(Cd, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_bath)
    if valid is False:
        ValidFlag = False

    # calculate bubble-related factors
    dB = Electro.BubbleThickness(I_line)
    (BCover, valid) = Electro.BubbleCoverage(Cd,  BR, I_line)
    if valid is False:
        ValidFlag = False

    # Calculate Electrochemical potential components
    (Erev, valid) = Electro.ReversiblePotential(Cd, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_bath)
    if valid is False:
        ValidFlag = False
    (Eca, valid) = Electro.AnodeConcOverVolt(Cd, T_bath, I_line)
    if valid is False:
        ValidFlag = False
    (Esa, valid) = Electro.AnodeSurfOverVolt(Cd, T_bath, I_line)
    if valid is False:
        ValidFlag = False
    (Ecc, valid) = Electro.CathodeConcOverVolt(BR, T_bath, I_line)
    if valid is False:
        ValidFlag = False

    # Calculate Ohmic reisstance components
    Rbath = Electro.BathRes(D, dB, kbath)
    Rbub = Electro.BubbleRes(BCover, dB, kbath)
    Rpath = Glob.Ran + Rbath + Rbub + Glob.Rca

    # Overall cell voltage as expressed in modified equation A1

    Vcell = Erev + Eca + Esa + Ecc + I_line*Rpath + I_line*Glob.Rext

    return (Vcell, ValidFlag)
