# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:07:07 2021

@author: Mark.Illingworth
"""

# Do imports first
import numpy as np

# constants for state amd model equations
CE = 0.945                 # current efficiency (fraction)
z = 6                      # electrons per mol Al2O3
Faraday = 96485.3329       # Farady constant s A / mol
Rgas = 8.3145              # Gas constant R
m = 1800                   # Assumed Bath mass (kg)
MAlumina = 101.9582        # Molar mass Alumina (g/mol)
r = 0.6                    # assumed ratio of alumina fed that doesnt dissolve
# note the paper had kdiss = 0.005, but this seemed too fast for sludge
kdiss = 0.002              # assumed slow rate for undissolved alumina sec-1
dcon = 1.03/(86400*0.5)    # assumed anode consumption cm/sec (10.3mm/12hr)
daccum = 3.63/(86400*1.5)  # assumed metal accumulation cm/sec (36.3mm/36hr)
# print(dcon, daccum)
dcon = dcon/126500         # convert anode consumption to cm/sec/A
daccum = daccum/126500     # convert metal accumulation to cm/sec/A
# print(dcon, daccum)

# Bath properties currently assumed, could be measurements
cAlF3 = 12.0             # assumed excess AlF3 in bath (wt%)
cCaF2 = 4.1              # assumed CaF2 in bath (wt%)
cMgF2 = 0.25             # assumed MgF2 in bath (wt%)
cLiF = 0.0               # assumed LiF in bath (wt%)
bTemp = 963.0            # assumed bath temperature (C)

# Cell design related constants
Aca = 115348             # geometric cathode area cm2##original number:115348
nAnode = 20#20              # number of anodes in cell
AAnode = 2*104652          # geometric anode area cm2
A = AAnode/nAnode        # geometric area of single anode cm2
f = 1.1382               # fanning factor for effective anode area
Ran = 3.02/1000000       # anode resistance ohms
Rca = 2.15/1000000       # cathode resistance ohms
Rext = 2.01/1000000      # external resistance ohms
AeOr = 1.80              # assumed alumina content at AE (wt%)

# Alumina shot mass and feed cycle time calculation
MAl = 26.98              # Molar mass Aluminium metal (g/mol)
AlDump = 2.9*2           # kg Alumina, Note BBA feeds both ends concurrently


def AluminaConsumption(Crrnt, CE):

    # Calcualte consumption rate kg/sec using global constants
    consumption = MAlumina/(z*Faraday)*Crrnt*CE/1000
    return consumption


def FeedCycleCalc(Crrnt, CE, Shotmass):

    # Calculate feed cycle time in seconds using global constants
    Consumption = AluminaConsumption(Crrnt, CE)
    FeedCycle = Shotmass/Consumption
    return FeedCycle


# Constants for acceptable ranges not all yet used
Al2O3max = 10.0
# practical minimum for gross bath alumina content would be at or near
# onset of AE, but allow it to go lower though not zero
Al2O3min = 1.0
AlF3max = 25.0
AlF3min = -10.0
CaF2max = 10.0
CaF2min = 0.0
LiFmax = 10.0
LiFmin = 0.0
MgF2max = 10.0
MgF2min = 0.0
BTempmax = 1200
BTempmin = 850

# Global vector of constants to call Vcell
BathConstVector = np.array([[cAlF3], [cCaF2], [cMgF2], [cLiF], [bTemp]])


def InRange(value, maximum, minimum):

    # Simple function to check a value is in range
    Valid = True
    ReturnVal = value
    if (value < minimum):
        Valid = False
        ReturnVal = minimum
    if (value > maximum):
        ReturnVal = maximum
        Valid = False
    return (ReturnVal, Valid)
