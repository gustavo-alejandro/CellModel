from IPython.terminal.pt_inputhooks import qt

import Global as Glob
import BathProperties as Bath
import ElectroChemistry as Electro
import CellVoltage as Volt
import pandas as pd
##
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib
# # #matplotlib.use('qt')
# matplotlib.use('Qt5Agg')
from mpl_toolkits.mplot3d import Axes3D
# Set the backend to %matplotlib qt
import matplotlib
matplotlib.use('TkAgg')

C_Al2O3 = 4.2
C_AlF3 =10.3
C_CaF2 =7
C_LiF = 0
C_MgF2 = 0.3
T_Bath = 963
I_line = 126000
ACD = 2.9

StateVector = [C_Al2O3, 0, ACD]
InputVector = [I_line, 0, 0]
ConstVector = [C_AlF3, C_CaF2, C_MgF2, C_LiF, T_Bath]

Vcell_total = Volt.Vcell(StateVector, InputVector, ConstVector)[0]

# Define a range of C_Al2O3 values from 2 to 8
C_Al2O3_values = np.linspace(1.9, 8, 100)
# Initialize empty lists to store the calculated values
Erev_values = []
Esa_values = []
Vbub_values = []
Eca_values = []
Ecc_values = []
V_cell_values = []


# Calculate the values for each C_Al2O3
for C_Al2O3 in C_Al2O3_values:
    BR = Bath.Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2)[0]
    Erev = Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]
    Erev_values.append(Erev)
    Esa = Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0]
    Esa_values.append(Esa)
    BCover = Electro.BubbleCoverage(C_Al2O3, BR, I_line)[0]
    Vbub = I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line), Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0])
    Vbub_values.append(Vbub)
    Eca = Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0]
    Eca_values.append(Eca)
    Ecc = Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0]
    Ecc_values.append(Ecc)
    Vcell = Erev+Esa+Eca+Ecc+Vbub
    V_cell_values.append(Vcell)


# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(C_Al2O3_values, Erev_values, label='Erev')
plt.plot(C_Al2O3_values, Esa_values, label='Esa')
plt.plot(C_Al2O3_values, Vbub_values, label='Vbub')
plt.plot(C_Al2O3_values, Eca_values, label='Eca')
plt.plot(C_Al2O3_values, Ecc_values, label='Ecc')
plt.plot(C_Al2O3_values, V_cell_values, label='Vcell')
plt.xlabel('C_Al2O3')
plt.ylabel('Voltage (V)')
plt.legend()
plt.title('Effect of C_Al2O3 on Erev, Esa, Vbub, Eca, and Ecc')
plt.grid(True)
plt.show()
