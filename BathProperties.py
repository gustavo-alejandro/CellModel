# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 09:21:14 2021

@author: Mark.Illingworth

# Collection of functions relating to physical bath properties
# Functions were extracted/adapted from two sources:
# - RTA PacOps ACD Calculator
# - Published paper 'Estimation of Spatial Alumina Concentration in an
#   Aluminum Reduction Cell Using a Multilevel State Observer,
#   Yuchen Yao et al, American Institute of Chemical Engineers
#   AIChE J, 63: 2806–2818, 2017'
# Calcualted outcomes should be consistent with the RTA PacOps current
# standards.
"""
import math as ma
import Global as Glob


def Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2):

    # All inputs in wt%
    # Function to evaluate Bath Ratio via eq A20a
    # Note that while all arguments need to be valid, at the moment
    # only checking alumina as it is simulated/filtered
    ValidFlag = True
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, Glob.Al2O3min)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    Base = 100.0 - (C_Al2O3 + C_CaF2 + C_LiF + C_MgF2)
    AlF3Ratio = C_AlF3 / Base
    Ratio = (1 - AlF3Ratio) / (2/3 + AlF3Ratio)
    return (Ratio, ValidFlag)


def AluminaSaturation(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath):
    """
    Skybakmoen, E., Solheim, A. & Sterten, Å. Alumina solubility in molten salt systems of interest for aluminum
    electrolysis and related phase diagram data. Metall Mater Trans B 28, 81–86 (1997).
    https://doi.org/10.1007/s11663-997-0129-9
    All inputs in wt%, temperature in C
    :param C_Al2O3:
    :param C_AlF3:
    :param C_CaF2:
    :param C_LiF:
    :param C_MgF2:
    :param T_Bath:
    :return:
    """
    # All inputs in wt%, temperature in C
    # Function to evaluate Alumina Saturation Solubility via eqa A5a
    # Note that while all arguments need to be valid, at the moment
    # only checking alumina as it is simulated/filtered
    ValidFlag = True
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, Glob.Al2O3min)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    a = 11.9 - 0.062*C_AlF3 - 0.0031*C_AlF3**2 - 0.2*C_CaF2
    a = a - 0.5*C_LiF - 0.3*C_MgF2 + 42.0*C_AlF3*C_LiF/(2000 + C_AlF3*C_LiF)
    B = 4.8 - 0.048*C_AlF3 + 2.2*(C_LiF**1.5)/(10 + C_LiF + 0.001*C_AlF3**3)
    sat = a*(T_Bath/1000)**B
    return (sat, ValidFlag)


def AluminaActivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath):
    """
    From paper titled 'Haupin, W. Interpreting the components of cell voltage'
    Eq. 10
    https://doi.org/10.1007/978-3-319-48156-2_21
    All inputs in wt%, temperature in C
    :param C_Al2O3:
    :param C_AlF3:
    :param C_CaF2:
    :param C_LiF:
    :param C_MgF2:
    :param T_Bath:
    :return:
    """
    # All inputs in wt%, temperature in C
    # Function to evaluate Alumina Activity via eq A4, A5a
    # Note that while all arguments need to be valid, at the moment
    # only checking alumina as it is simulated/filtered
    ValidFlag = True
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, Glob.Al2O3min)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    Sat, valid = AluminaSaturation(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)
    if valid is False:
        ValidFlag = False
    ROS = C_Al2O3/Sat
    act = -0.03791*ROS + 2.364*(ROS**2) - 2.194*(ROS**3) + 0.8686*(ROS**4)
    return (act, ValidFlag)


def Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath):
    """
    There is a reference empirical cryolite bath conductivity eq. in ref
    https://doi.org/10.1007/BF02915051
    All inputs in wt%, temperature in C. Conductivity in S/cm
    :param C_Al2O3:
    :param C_AlF3:
    :param C_CaF2:
    :param C_LiF:
    :param C_MgF2:
    :param T_Bath:
    :return:
    """
    # All inputs in wt%, temperature in C
    # Function to evaluate Bath Conductivity via eq A17a,b,c
    # Note that while all arguments need to be valid, at the moment
    # only checking alumina as it is simulated/filtered
    ValidFlag = True
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, Glob.Al2O3min)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    C_NaF = 0.6*(100 - C_AlF3 - C_Al2O3 - C_CaF2 - C_LiF - C_MgF2)
    BWR = (C_NaF/41.99 + C_LiF/25.94 - C_MgF2/62.31)
    BWR = BWR/2/((C_AlF3 + 2*C_NaF/3)/83.98)
    condexp = 1.9362 + 0.3092*BWR - 0.004132*C_CaF2 - 0.01755*C_Al2O3
    condexp = condexp + 0.008123*C_LiF - 0.00398*C_MgF2
    condexp = condexp - 1751.1/(T_Bath + 273.15)
    conductivity = ma.exp(condexp)
    return (conductivity, ValidFlag)
