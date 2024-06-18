# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 11:48:07 2021

@author: Mark.Illingworth

# Collection of functions relating to Electrical Potential
# Functions that were extracted/adapted from two sources:
# - RTA PacOps ACD Calculator
# - Published paper 'Estimation of Spatial Alumina Concentration in an
#   Aluminum Reduction Cell Using a Multilevel State Observer,
#   Yuchen Yao et al, American Institute of Chemical Engineers
#   AIChE J, 63: 2806–2818, 2017'
# Calculated outcomes should be consistent with the RTA PacOps current
# standards.
"""

import math as ma
import Global as Glob
import BathProperties as Bath


def ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath):
    """
    From paper titled 'Haupin, W. Interpreting the components of cell voltage'
    Eq. 5
    https://doi.org/10.1007/978-3-319-48156-2_21
    :param C_Al2O3:
    :param C_AlF3:
    :param C_CaF2:
    :param C_LiF:
    :param C_MgF2:
    :param T_Bath:
    :return:
    """
    # All inputs in wt%, temperature in C
    # Function to evaluate Reversible Potential Erev via eq A2a,A3a
    ValidFlag = True
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, Glob.Al2O3min)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    E0 = 1.896 - 0.000572*(T_Bath + 273.15)
    act, valid = Bath.AluminaActivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)
    if valid is False:
        ValidFlag = False
    logact = ma.log((1/act**2))
    Erev = E0 + Glob.Rgas*(T_Bath + 273.15)*logact/(12*Glob.Faraday)
    return (Erev, ValidFlag)
    #return Erev


def AnodeCurrentDensity(I_line):
    # simple function to evaluate anode current density
    Icell = I_line/Glob.f/Glob.AAnode
    return Icell


def CriticalCurrent(C_Al2O3, T_Bath, I_line):
    """
    Equation from GRJOTHEIM, Kai; WELCH, Barry J. Aluminium Smelter Technology--a Pure and Applied Approach.
    Aluminium-Verlag, P. O. Box 1207, Konigsallee 30, D 4000 Dusseldorf 1, FRG, 1988., 1988.
    Chapter 5, Equation 13, page 131.
    :param C_Al2O3:
    :param T_Bath:
    :param I_line:
    :return:
    """
    # All inputs in wt%, temperature in C, current in amps
    # Function to evaluate Critical Current Density icell via eq A7a
    ValidFlag = True
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, Glob.Al2O3min)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    tmp1 = ((C_Al2O3/100)**0.5-0.04)*(Glob.A**-0.1)
    icrit = 0.0001*((550000+1800*((T_Bath + 273.15)-1323))*tmp1)
    return (icrit, ValidFlag)


def AnodeConcOverVolt(C_Al2O3, T_Bath, I_line):
    """
    From paper titled 'Haupin, W. Interpreting the components of cell voltage'
    Eq. 23
    https://doi.org/10.1007/978-3-319-48156-2_21
    :param C_Al2O3:
    :param T_Bath:
    :param I_line:
    :return:
    """
    # All inputs in wt%, temperature in C, current in amps
    # Function to evaluate Anode Conc Overvoltage via eq A6a, A8a
    ValidFlag = True
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, Glob.Al2O3min)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    icell = AnodeCurrentDensity(I_line)
    ic, valid = CriticalCurrent(C_Al2O3, T_Bath, I_line)
    if valid is False:
        ValidFlag = False
    if icell >= (ic - 0.01):
        print(icell, ic, C_Al2O3, T_Bath, I_line)
        # as ic gets closer to icell, the ratio gets infinitely large
        # This arbitrary value results in AE voltage of ~ 30V
        log_ic = 500
    else:
        log_ic = ma.log(ic/(ic - icell))
    Eca = (T_Bath + 273.15)*log_ic/23209 # Value 23210 is R/(2F)
    return (Eca, ValidFlag)


def AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line):
    """
    From paper titled 'Haupin, W. Interpreting the components of cell voltage'
    Eq. 26
    https://doi.org/10.1007/978-3-319-48156-2_21
    :param C_Al2O3:
    :param T_Bath:
    :param I_line:
    :return:
    """
    # All inputs in wt%, temperature in C, current in amps
    # Function to evaluate Anode Surface Overvoltage via eq A11a, A12a
    ValidFlag = True
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, Glob.Al2O3min)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    icell = AnodeCurrentDensity(I_line)
    ir = 0.0029*(C_Al2O3**0.56)
    log_ir = ma.log(icell/ir)
    Esa = (T_Bath + 273.15)*log_ir/12533
    return (Esa, ValidFlag)


def CathodeConcOverVolt(Ratio, T_Bath, I_line):
    """
    From paper titled 'Haupin, W. Interpreting the components of cell voltage'
    Eq. 28
    https://doi.org/10.1007/978-3-319-48156-2_21
    :param Ratio:
    :param T_Bath:
    :param I_line:
    :return:
    """
    # All inputs in wt%, temperature in C, current in amps
    # Function to evaluate Cathode Conc Overvoltage via eq A13a,b
    ica = I_line/Glob.Aca
    if ica <= 0:
        log_ica = 0
        ValidFlag = False
    else:
        log_ica = ma.log(3.891050584*ica)
        ValidFlag = True
    Ecc = 0.00005744815304*(T_Bath + 273.15)*(1.375 - 0.25*Ratio)*log_ica
    return (Ecc, ValidFlag)


def BubbleThickness(I_line):
    """
    From paper titled 'Haupin, W. Interpreting the components of cell voltage'
    Eq. 30
    https://doi.org/10.1007/978-3-319-48156-2_21
    :param I_line:
    :return:
    """
    # All inputs in wt%, temperature in C, current in amps
    # Function to evaluate Bubble Thickness via eq A16 in cm
    icell = AnodeCurrentDensity(I_line)
    #icell = I_line / Glob.f / Glob.AAnode
    #f = 1.1382  # fanning factor for effective anode area
    #AAnode = 2 * 104652  # geometric anode area cm2
    #icell = I_line / f / AAnode
    db = (0.5517 + icell)/(1 + 2.167*icell)
    return db


def BathRes(ACD, db, k_bath):

    # inputs in units from previous functions (cm based)
    # Function to evaluate Bath Resistance via eq A15
    # Note RTA doesn't subtract the bubble thickness layer form this
    # calc, so will exclude it for now
    f = 1.1382e0  # fanning factor for effective anode area
    AAnode = 2 * 104652e0  # geometric anode area cm2
    R_Bath = (ACD - db) / (k_bath * f * AAnode)
    #R_Bath = (ACD - db) / (k_bath * Glob.f * Glob.AAnode)
    #R_Bath = (ACD - 0*db)/(k_bath*Glob.f*Glob.AAnode)
    return R_Bath


def BubbleCoverage(C_Al2O3, Ratio, I_line):
    """
    From paper titled 'Haupin, W. Interpreting the components of cell voltage'
    Eq. 31
    https://doi.org/10.1007/978-3-319-48156-2_21
    :param C_Al2O3:
    :param Ratio:
    :param I_line:
    :return:
    """
    # All inputs in wt%, temperature in C, current in amps
    # Function to evaluate Bubble Coverage via eq A18a
    AEOR = Glob.AeOr
    ValidFlag = True
    # bubble coverage tends towards 100% as C_Al2O3 approaches AEOR
    # limit the impact to a practical value that will apply for
    # all alumina content less than that
    adjusted, valid = Glob.InRange(C_Al2O3, Glob.Al2O3max, AEOR+0.01)
    if valid is False:
        C_Al2O3 = adjusted
        ValidFlag = False
    icell = AnodeCurrentDensity(I_line)
    BRc = (0.4322 - 0.3781*Ratio)/(1 + 1.637*Ratio)
    Aluminac = (0.431 - 0.1437*(C_Al2O3 - AEOR))/(1 + 7.353*(C_Al2O3 - AEOR))
    Coverage = 0.509 + 0.1823*icell - 0.1723*(icell**2) + 0.05504*(icell**3)
    Coverage = Coverage + BRc + Aluminac
    Coverage = 0.9*Coverage
    return (Coverage, ValidFlag)


def BubbleRes(Coverage, db, k_bath):
    """
    From paper titled 'Haupin, W. Interpreting the components of cell voltage'
    Eq. 32
    https://doi.org/10.1007/978-3-319-48156-2_21
    :param Coverage:
    :param db:
    :param k_bath:
    :return:R)_
    """
    # inputs in units from previous functions (cm based)
    # Function to evaluate Bubble Resistance via eq A19
    f = 1.1382e0  # fanning factor for effective anode area
    AAnode = 2 * 104652e0  # geometric anode area cm2
    #R_Bub = (db*Coverage)/(k_bath*Glob.f*Glob.AAnode*(1-Coverage))
    R_Bub = (db * Coverage) / (k_bath * f * AAnode * (1 - Coverage))
    return R_Bub
