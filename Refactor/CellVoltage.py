# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:38:06 2021

@author: Mark.Illingworth
"""
from Refactor.Properties import CellProperties


def VIcell(cAl2O3: float, ACD: float, deltaI: float, I0: float,
           Cell: CellProperties):

    # Non-linear function to define H
    # Note that while Cun and Cd are related in the state equations
    # Cun is not used to calculate Vcell

    # Valid Flag will be set to false if any of the properties or
    # electrochemical functions returns an invalid argument response
    ValidVolt = True

    # Update the Alumina concentration property with the requested state value
    Cell.cAl2O3.value = cAl2O3

    # Calculate Actual/Total Current
    Icell = I0 + deltaI

    # Calculate Electrochemical potential components
    (Erev, valid) = Cell.ReversiblePotential
    if valid is False:
        ValidVolt = False
        raise ValueError("Out of Range Reversible Potential")

    (Eca, valid) = Cell.AnodeConcOverVolt(Icell)
    if valid is False:
        ValidVolt = False
        raise ValueError("Out of Range Anode Concentration Overvoltage")

    (Esa, valid) = Cell.AnodeSurfOverVolt(Icell)
    if valid is False:
        ValidVolt = False
        raise ValueError("Out of Range Anode Surface Overvoltage")

    (Ecc, valid) = Cell.CathodeConcOverVolt(Icell)
    if valid is False:
        ValidVolt = False
        raise ValueError("Out of Range Cathode Concentration Overvoltage")

    # Calculate Ohmic resistance components
    Rbath, valid = Cell.BathRes(ACD)
    if valid is False:
        ValidVolt = False
        raise ValueError("Out of Range Bath Resistance")

    Rbub, valid = Cell.BubbleRes(Icell)
    if valid is False:
        ValidVolt = False
        raise ValueError("Out of Range Bath Resistance")

    # Calculate total cell resistance
    Rpath = Cell.Ran + Rbath + Rbub + Cell.Rca

    # Overall cell voltage as expressed in modified equation A1
    Vcell = Erev + Eca + Esa + Ecc + Icell*(Rpath + Cell.Rext)

    # return (Vcell, Icell, ValidVolt)
    return (Vcell, Icell, Erev, Eca, Esa, Ecc, Rpath, ValidVolt)
