from IPython.terminal.pt_inputhooks import qt

import Global as Glob
import BathProperties as Bath
import ElectroChemistry as Electro
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
"""
Vcell = Erev + Eca + Esa + Ecc + I_line*Rpath + I_line*Glob.Rext
"""

C_Al2O3 = 4.2
C_AlF3 =10.3
C_CaF2 =7
C_LiF = 0
C_MgF2 = 0.3
T_Bath = 963
I_line = 126000
ACD = 2.9



BR = Bath.Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2)[0]
Erev = Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]
Eca = Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0]
Esa = Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0]
Ecc = Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0]

db = Electro.BubbleThickness(I_line)
kbath = Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]
BCover = Electro.BubbleCoverage(C_Al2O3, BR, I_line)[0]

# Calculate Ohmic reisstance components
Rbath = Electro.BathRes(ACD, db, kbath)
Rbub = Electro.BubbleRes(BCover, db, kbath)
Rpath = Glob.Ran + Rbath + Rbub + Glob.Rca
Van = I_line*Glob.Ran
Vbath = I_line*Rbath
Vbub = I_line*Rbub
Vca = I_line*Glob.Rca
Vdroppath = I_line*Rpath
Vdropext = I_line*Glob.Rext
# ibottomfraction = 0.8275
# ibottom = I_line*ibottomfraction/Glob.nAnode
# ibottomcd = ibottom/Glob.A
# Vbath = ibottomcd*(ACD-db)*(1/kbath)
# Vdroppath = Van+Vbath+Vbub+Vca
Rpath = Glob.Ran + Rbath + Rbub + Glob.Rca

Vcell = Erev + Eca + Esa + Ecc + Vdroppath + I_line*Glob.Rext
Volt_table = pd.DataFrame([{'Erev': Erev, 'Eca': Eca, 'Esa': Esa, 'Ecc': Ecc, 'Vdroppath': Vdroppath, 'Vdropext': Vdropext}])
print(f"Cell voltage: {Volt_table.sum(axis=1)[0]}")


# Define a range of C_Al2O3 values from 2 to 8
C_Al2O3_values = np.linspace(1.9, 8, 100)

# Initialize empty lists to store the calculated values
Erev_values = []
Esa_values = []
Vbub_values = []
Eca_values = []
Ecc_values = []
#V_cell_values = []


# Calculate the values for each C_Al2O3
for C_Al2O3 in C_Al2O3_values:
    BR = Bath.Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2)[0]
    Erev_values.append(Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0])
    Esa_values.append(Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0])
    BCover = Electro.BubbleCoverage(C_Al2O3, BR, I_line)[0]
    Vbub_values.append(I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line), Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]))
    Eca_values.append(Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0])
    Ecc_values.append(Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0])

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(C_Al2O3_values, Erev_values, label='Erev')
plt.plot(C_Al2O3_values, Esa_values, label='Esa')
plt.plot(C_Al2O3_values, Vbub_values, label='Vbub')
plt.plot(C_Al2O3_values, Eca_values, label='Eca')
plt.plot(C_Al2O3_values, Ecc_values, label='Ecc')
plt.xlabel('C_Al2O3')
plt.ylabel('Voltage (V)')
plt.legend()
plt.title('Effect of C_Al2O3 on Erev, Esa, Vbub, Eca, and Ecc')
plt.grid(True)
plt.show()


#
#
# # Define a range of C_Al2O3 and C_AlF3 values from 2 to 8
# C_Al2O3_values = np.linspace(2, 8, 30)
# C_AlF3_values = np.linspace(2, 8, 30)
#
# # Create a meshgrid of C_Al2O3 and C_AlF3 values
# C_Al2O3_grid, C_AlF3_grid = np.meshgrid(C_Al2O3_values, C_AlF3_values)
#
# # Flatten the meshgrid arrays to 1D arrays
# C_Al2O3_flat = C_Al2O3_grid.flatten()
# C_AlF3_flat = C_AlF3_grid.flatten()
#
# # Initialize empty arrays to store the calculated values
# Erev_values = np.zeros_like(C_Al2O3_flat)
# Esa_values = np.zeros_like(C_Al2O3_flat)
# Vbub_values = np.zeros_like(C_Al2O3_flat)
# Eca_values = np.zeros_like(C_Al2O3_flat)
# Ecc_values = np.zeros_like(C_Al2O3_flat)
#
# # Calculate the values for each combination of C_Al2O3 and C_AlF3
# for i in range(len(C_Al2O3_flat)):
#     C_Al2O3 = C_Al2O3_flat[i]
#     C_AlF3 = C_AlF3_flat[i]
#     BR = Bath.Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2)[0]
#     Erev_values[i] = Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]
#     Esa_values[i] = Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0]
#     BCover = Electro.BubbleCoverage(C_Al2O3, BR, I_line)[0]
#     Vbub_values[i] = I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line), Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0])
#     Eca_values[i] = Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0]
#     Ecc_values[i] = Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0]
#
# # Reshape the 1D arrays back to the original grid shape
# Erev_values = Erev_values.reshape(C_Al2O3_grid.shape)
# Esa_values = Esa_values.reshape(C_Al2O3_grid.shape)
# Vbub_values = Vbub_values.reshape(C_Al2O3_grid.shape)
# Eca_values = Eca_values.reshape(C_Al2O3_grid.shape)
# Ecc_values = Ecc_values.reshape(C_Al2O3_grid.shape)
#
# # Plot the 3D surfaces
# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111, projection='3d')
#
# # Plot the surfaces
# ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Erev_values, label='Erev')
# ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Esa_values, label='Esa')
# ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Vbub_values, label='Vbub')
# ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Eca_values, label='Eca')
# ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Ecc_values, label='Ecc')
#
# # Set axis labels
# ax.set_xlabel('C_Al2O3')
# ax.set_ylabel('C_AlF3')
# ax.set_zlabel('Voltage (V)')
#
# # Add a legend
# #ax.legend()
#
# # Set plot title
# plt.title('Effect of C_Al2O3 and C_AlF3 on Erev, Esa, Vbub, Eca, and Ecc')
#
# # Show the plot
# #plt.show()
# plt.show()