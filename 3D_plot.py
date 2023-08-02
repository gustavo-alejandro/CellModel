
# Define a range of C_Al2O3 and C_AlF3 values from 2 to 8
C_Al2O3_values = np.linspace(2, 8, 30)
C_AlF3_values = np.linspace(2, 8, 30)

# Create a meshgrid of C_Al2O3 and C_AlF3 values
C_Al2O3_grid, C_AlF3_grid = np.meshgrid(C_Al2O3_values, C_AlF3_values)

# Flatten the meshgrid arrays to 1D arrays
C_Al2O3_flat = C_Al2O3_grid.flatten()
C_AlF3_flat = C_AlF3_grid.flatten()

# Initialize empty arrays to store the calculated values
Erev_values = np.zeros_like(C_Al2O3_flat)
Esa_values = np.zeros_like(C_Al2O3_flat)
Vbub_values = np.zeros_like(C_Al2O3_flat)
Eca_values = np.zeros_like(C_Al2O3_flat)
Ecc_values = np.zeros_like(C_Al2O3_flat)

# Calculate the values for each combination of C_Al2O3 and C_AlF3
for i in range(len(C_Al2O3_flat)):
    C_Al2O3 = C_Al2O3_flat[i]
    C_AlF3 = C_AlF3_flat[i]
    BR = Bath.Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2)[0]
    Erev_values[i] = Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]
    Esa_values[i] = Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0]
    BCover = Electro.BubbleCoverage(C_Al2O3, BR, I_line)[0]
    Vbub_values[i] = I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line), Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0])
    Eca_values[i] = Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0]
    Ecc_values[i] = Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0]

# Reshape the 1D arrays back to the original grid shape
Erev_values = Erev_values.reshape(C_Al2O3_grid.shape)
Esa_values = Esa_values.reshape(C_Al2O3_grid.shape)
Vbub_values = Vbub_values.reshape(C_Al2O3_grid.shape)
Eca_values = Eca_values.reshape(C_Al2O3_grid.shape)
Ecc_values = Ecc_values.reshape(C_Al2O3_grid.shape)

# Plot the 3D surfaces
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the surfaces
ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Erev_values, label='Erev')
ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Esa_values, label='Esa')
ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Vbub_values, label='Vbub')
ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Eca_values, label='Eca')
ax.plot_surface(C_Al2O3_grid, C_AlF3_grid, Ecc_values, label='Ecc')

# Set axis labels
ax.set_xlabel('C_Al2O3')
ax.set_ylabel('C_AlF3')
ax.set_zlabel('Voltage (V)')

# Add a legend
#ax.legend()

# Set plot title
plt.title('Effect of C_Al2O3 and C_AlF3 on Erev, Esa, Vbub, Eca, and Ecc')

# Show the plot
plt.show()