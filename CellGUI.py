import Global as Glob
import BathProperties as Bath
import ElectroChemistry as Electro
import CellVoltage as Volt
import tkinter as tk
from tkinter import ttk
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

##Create GUI window
root = tk.Tk()
root.title("Aluminium Cell Model")
plot_frame = ttk.Frame(root)
plot_frame.grid(row=1, column=7, rowspan=6, padx=10, pady=10)

###Inputs from GUI
def setvars(C_Al2O3_val, C_AlF3_val, C_CaF2_val, C_LiF_val, C_MgF2_val, T_Bath_val, I_line_val, ACD_val):
    C_Al2O3 = tk.DoubleVar(value=C_Al2O3_val)
    C_AlF3 = tk.DoubleVar(value=C_AlF3_val)
    C_CaF2 = tk.DoubleVar(value=C_CaF2_val)
    C_LiF =  tk.DoubleVar(value=C_LiF_val)
    C_MgF2 = tk.DoubleVar(value=C_MgF2_val)
    T_Bath = tk.DoubleVar(value=T_Bath_val)
    I_line = tk.DoubleVar(value=I_line_val)
    ACD = tk.DoubleVar(value=ACD_val)

    StateVector = [C_Al2O3, 0, ACD]
    InputVector = [I_line, 0, 0]
    ConstVector = [C_AlF3, C_CaF2, C_MgF2, C_LiF, T_Bath]
    return StateVector, InputVector, ConstVector

StateVector, InputVector, ConstVector = setvars(C_Al2O3_val=4.2, C_AlF3_val=10.3, C_CaF2_val=7, C_LiF_val=0, C_MgF2_val=.3, T_Bath_val=963, I_line_val=126000, ACD_val=2.9)


#    StateVector[0] = StateVector[0].get()

##Insert aluminium cell image
image_path = "C:/Users/n11675250/OneDrive - Queensland University of Technology/Aluminium Cell/E/hhprg/cellschema/images/left_riser_gesamt1.png"
image = tk.PhotoImage(file=image_path)
image = image.subsample(2,2)
image_label = ttk.Label(root, image=image, width=100)
image_label.grid(row=0, column=0, columnspan=3, padx=10, pady=10)

###Insert sliders block

attributes_frame = ttk.LabelFrame(root, text="Bath comp wt%")
attributes_frame.grid(row=0, column=5, columnspan=1, padx=5, pady=5, sticky="nw")
attributes_labels = ['C_Al2O3', 'C_AlF3', 'C_CaF2', 'C_MgF2', 'C_LiF', 'Bath Temperature (K)']
attributes_vars = [StateVector[0], ConstVector[0], ConstVector[1], ConstVector[2], ConstVector[3], ConstVector[4]]
slider_ranges = [(0, 11), (0, 11), (0, 11), (0, 11), (0, 11), (800, 1000)]
for idx, label in enumerate(attributes_labels):
    ttk.Label(attributes_frame, text=label, width=10).grid(row=idx, column=0, padx=5, pady=5)
    slider = ttk.Scale(
    attributes_frame, variable=attributes_vars[idx], from_=slider_ranges[idx][0], to=slider_ranges[idx][1], length=50, orient="horizontal")
    slider.grid(row=idx, column=1, padx=5, pady=5)
    slider.set(attributes_vars[idx].get())  # Set the slider value to the initial attribute value
    #slider.bind("<ButtonRelease-1>", on_slider_release())  # Bind the slider release event
    ttk.Label(attributes_frame, textvariable=attributes_vars[idx], width=4).grid(row=idx, column=2, padx=5, pady=5)


# def on_slider_release():
#     update_results()
#
#


#List of GUI Elements
#Cell input frame GUI

cell_frame = ttk.LabelFrame(root, text="Cell input data")
cell_frame.grid(row=1, column=0, padx=5, pady=10, sticky="nw")

ttk.Label(cell_frame, text="Current", width=10).grid(row=0, column=0, padx=5, pady=5)
cell_current_entry = ttk.Entry(cell_frame, textvariable=InputVector[0], width=6)
cell_current_entry.grid(row=0, column=1, padx=1, pady=1)

ttk.Label(cell_frame, text="ACD", width=10).grid(row=1, column=0, padx=5, pady=5)
ACD_entry = ttk.Entry(cell_frame, textvariable=StateVector[2], width=6)
ACD_entry.grid(row=1, column=1, padx=1, pady=1)
###

###BathProperties elements

resistivity_frame = ttk.LabelFrame(root, text="Bath Resistivity")
resistivity_frame.grid(row=4, column=5, padx=10, pady=10, sticky="w")

resistivity_label = ttk.Label(resistivity_frame, text="Bath Resistivity: ")
resistivity_label.grid(row=0, column=0, padx=5, pady=5)

Equil_potential_frame = ttk.LabelFrame(root, text="Equil_potential")
Equil_potential_frame.grid(row=1, column=5, padx=10, pady=10, sticky="nw")

Equil_potential_label = ttk.Label(Equil_potential_frame, text="Equil_potential: ")
Equil_potential_label.grid(row=3, column=5, padx=5, pady=5)

bath_ratio_frame = ttk.LabelFrame(root, text="bath ratio")
bath_ratio_frame.grid(row=2, column=5, padx=10, pady=10, sticky="nw")

bath_ratio_label = ttk.Label(bath_ratio_frame, text="bath_ratio: ")
bath_ratio_label.grid(row=0, column=0, padx=5, pady=5)

# def update_results():
#     #StateVector[0].set(attributes_vars[0])
#     C_Al2O3 = tk.DoubleVar(value=StateVector[0].get())
#     print(StateVector[0].get())
#     resistivity_label.config(text=f"Bath Resistivity: {StateVector[0].get():.4f}")
#
# ttk.Button(root, text="Calculate", command=update_results()).grid(row=5, column=0, columnspan=3, padx=10, pady=10)



root.mainloop()