import tkinter as tk
from tkinter import ttk
import Global as Glob
import BathProperties as Bath
import ElectroChemistry as Electro
import CellVoltage as Volt
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

StateVector = [4.2, 0, 2.9]
InputVector = [126000, 0, 0]
ConstVector = [10.3, 7, .3, 0, 963]

C_AlF3 = ConstVector[0]
C_CaF2 = ConstVector[1]
C_LiF = ConstVector[3]
C_MgF2 = ConstVector[2]
T_Bath = ConstVector[4]
I_line = InputVector[0]
ACD = StateVector[2]

# Initialize empty lists to store the calculated values
Erev_values = []
Eca_values = []
Esa_values = []
Ecc_values = []
Vbath_values = []
Vbub_values = []
V_an_values = []
V_ca_values = []
V_ext_values = []
V_cell_values = []
# Define a range of C_Al2O3 values from 2 to 8
time = np.linspace(0, 1400, 701)
C_Al2O3_values = []

for t in time:
    C_Al2O3 = 7.5977e-18 * pow(t, 6) - 3.8081e-14 * pow(t, 5) + 7.1245e-11 * pow(t, 4) \
              - 6.0775e-08 * pow(t, 3) + 2.2483e-05 * pow(t, 2) - 1.2591e-03 * pow(t, 1) + 7.2423e-01
    C_Al2O3_values.append(C_Al2O3)
    BR = Bath.Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2)[0]
    Erev = Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]
    Erev_values.append(Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0])

    Esa = Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0]
    Esa_values.append(Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0])
    BCover = Electro.BubbleCoverage(C_Al2O3, BR, I_line)[0]

    Vbub = I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line),
                                      Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[
                                          0])
    Vbub_values.append(I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line),
                                                  Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[
                                                      0]))

    Eca = Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0]

    Eca_values.append(Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0])

    Ecc = Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0]
    Ecc_values.append(Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0])

    dB = Electro.BubbleThickness(I_line)
    kbath = Bath.Conductivity(C_Al2O3, ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2],
                              ConstVector[4])[0]
    Rbath = Electro.BathRes(ACD, dB, kbath)

    Vbath = I_line * Rbath
    Vbath_values.append(I_line * Rbath)

    V_ca = I_line * Glob.Ran
    V_an = I_line * Glob.Rca
    V_ext = I_line * Glob.Rext

    V_ca_values.append(I_line * Glob.Ran)
    V_an_values.append(I_line * Glob.Rca)
    V_ext_values.append(I_line * Glob.Rext)

    V_cell_values.append(Erev + Esa + Vbub + Eca + Ecc + Vbath + V_ca + V_an + V_ext)


