
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



