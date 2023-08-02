import tkinter as tk
from tkinter import ttk
import BathProperties as Bath
import ElectroChemistry as Electro
import CellVoltage as Volt

StateVector = [4.2, 0, 2.9]
InputVector = [126000, 0, 0]
ConstVector = [10.3, 7, .3, 0, 963]

##Create GUI window
root = tk.Tk()
root.title("Aluminium Cell Model")
plot_frame = ttk.Frame(root)
plot_frame.grid(row=1, column=7, rowspan=6, padx=10, pady=10)

# Variables to store the updated values
C_Al2O3 = tk.DoubleVar(value=StateVector[0])
C_AlF3 = tk.DoubleVar(value=ConstVector[0])
C_CaF2 = tk.DoubleVar(value=ConstVector[1])
C_LiF = tk.DoubleVar(value=ConstVector[3])
C_MgF2 = tk.DoubleVar(value=ConstVector[2])
T_Bath = tk.DoubleVar(value=ConstVector[4])

I_line = tk.DoubleVar(value=InputVector[0])
ACD = tk.DoubleVar(value=StateVector[2])

##Insert aluminium cell image
image_path = "C:/Users/n11675250/OneDrive - Queensland University of Technology/Aluminium Cell/E/hhprg/cellschema/images/left_riser_gesamt1.png"
image = tk.PhotoImage(file=image_path)
image = image.subsample(2,2)
image_label = ttk.Label(root, image=image, width=100)
image_label.grid(row=0, column=0, columnspan=3, padx=10, pady=10)

#Sliders
bath_comp_frame = ttk.LabelFrame(root, text="Bath comp wt%")
bath_comp_frame.grid(row=0, column=5, columnspan=1, padx=5, pady=5, sticky="nw")
slider_labels = ['C_Al2O3', 'C_AlF3', 'C_CaF2', 'C_MgF2', 'C_LiF', 'T_Bath']
slider_vars = [C_Al2O3, C_AlF3, C_CaF2, C_MgF2, C_LiF, T_Bath]
slider_ranges = [(0, 11), (0, 11), (0, 11), (0, 11), (0, 11), (800, 1000)]

for idx, label in enumerate(slider_labels):
    ttk.Label(bath_comp_frame, text=label, width=10).grid(row=idx, column=0, padx=5, pady=5)
    slider = ttk.Scale(
    bath_comp_frame, variable=slider_vars[idx], from_=slider_ranges[idx][0], to=slider_ranges[idx][1], length=50, orient="horizontal")
    slider.grid(row=idx, column=1, padx=5, pady=5)
    slider.set(slider_vars[idx].get())  # Set the slider value to the initial attribute value
    #slider.bind("<ButtonRelease-1>", on_slider_release())  # Bind the slider release event
    ttk.Label(bath_comp_frame, textvariable=slider_vars[idx], width=4).grid(row=idx, column=2, padx=5, pady=5)


# BathProperties calcs

cryolite_ratio_frame = ttk.LabelFrame(root, text="Bath ratio")
cryolite_ratio_frame.grid(row=1, column=5, padx=10, pady=10, sticky="w")
cryolite_ratio_label = ttk.Label(cryolite_ratio_frame, text="Bath ratio: ")
cryolite_ratio_label.grid(row=0, column=0, padx=5, pady=5)

alumina_sat_frame = ttk.LabelFrame(root, text="Alumina saturation")
alumina_sat_frame.grid(row=2, column=5, padx=10, pady=10, sticky="w")
alumina_sat_label = ttk.Label(alumina_sat_frame, text="Sat: ")
alumina_sat_label.grid(row=0, column=0, padx=5, pady=5)

#Cell input frame GUI

cell_frame = ttk.LabelFrame(root, text="Cell input data")
cell_frame.grid(row=1, column=0, padx=5, pady=10, sticky="nw")

ttk.Label(cell_frame, text="Current", width=10).grid(row=0, column=0, padx=5, pady=5)
cell_current_entry = ttk.Entry(cell_frame, textvariable=I_line, width=6)
cell_current_entry.grid(row=0, column=1, padx=1, pady=1)

ttk.Label(cell_frame, text="ACD", width=10).grid(row=1, column=0, padx=5, pady=5)
ACD_entry = ttk.Entry(cell_frame, textvariable=ACD, width=6)
ACD_entry.grid(row=1, column=1, padx=1, pady=1)

#Output
V_cell_frame = ttk.LabelFrame(root, text="V_cell")
V_cell_frame.grid(row=2, column=0, padx=10, pady=10, sticky="w")
V_cell_label = ttk.Label(V_cell_frame, text="V_cell: ")
V_cell_label.grid(row=0, column=0, padx=5, pady=5)


def update_values():
    # Update StateVector and ConstVector with the current slider values
    StateVector[0] = C_Al2O3.get()
    ConstVector[0] = C_AlF3.get()
    ConstVector[1] = C_CaF2.get()
    ConstVector[2] = C_MgF2.get()
    ConstVector[3] = C_LiF.get()
    ConstVector[4] = T_Bath.get()
    InputVector[0] = I_line.get()
    StateVector[2] = ACD.get()

    Vcell_total = Volt.Vcell(StateVector, InputVector, ConstVector)[0]
    print(Vcell_total)

    #Bath ratio calc to display in GUI
    bath_ratio = Bath.Ratio(StateVector[0], ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2])[0]
    cryolite_ratio_label.config(text=f"bath_ratio: {bath_ratio:.4f}")

    #Alumina saturation
    alumina_sat = Bath.AluminaSaturation(StateVector[0], ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2], ConstVector[4])[0]
    alumina_sat_label.config(text=f"Alumina sat: {alumina_sat:.4f}")

    V_cell_label.config(text=f"V_cell: {Vcell_total:.4f}")

update_button = ttk.Button(root, text="Update Values", command=update_values)
update_button.grid(row=len(slider_labels), column=0, columnspan=2, pady=10)

root.mainloop()
