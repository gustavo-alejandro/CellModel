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

StateVector = [2.35, 0, 2.9]
InputVector = [126000, 0, 0]
ConstVector = [10.3, 7, .3, 0, 963]

##Create GUI window
root = tk.Tk()
root.title("Aluminium Cell Model")
#root.geometry('800x800')
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
image = image.subsample(1,1)
image_label = ttk.Label(root, image=image, width=100)
image_label.grid(row=0, column=2, columnspan=3, padx=10, pady=10)

#Sliders
bath_col=1
bath_comp_frame = ttk.LabelFrame(root, text="Bath comp wt%")
bath_comp_frame.grid(row=0, column=bath_col, columnspan=1, padx=5, pady=5, sticky="w")
slider_labels = ['C_Al2O3', 'C_AlF3', 'C_CaF2', 'C_MgF2', 'C_LiF', 'T_Bath']
slider_vars = [C_Al2O3, C_AlF3, C_CaF2, C_MgF2, C_LiF, T_Bath]
slider_ranges = [(0, 50), (0, 11), (0, 11), (0, 11), (0, 11), (800, 1000)]

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
cryolite_ratio_frame.grid(row=1, column=bath_col, padx=10, pady=10, sticky="w")
cryolite_ratio_label = ttk.Label(cryolite_ratio_frame, text="Bath ratio: ")
cryolite_ratio_label.grid(row=0, column=0, padx=5, pady=5)

alumina_sat_frame = ttk.LabelFrame(root, text="Alumina saturation")
alumina_sat_frame.grid(row=2, column=bath_col, padx=10, pady=10, sticky="w")
alumina_sat_label = ttk.Label(alumina_sat_frame, text="Sat: ")
alumina_sat_label.grid(row=0, column=0, padx=5, pady=5)

alumina_activ_frame = ttk.LabelFrame(root, text="Alumina activity")
alumina_activ_frame.grid(row=3, column=bath_col, padx=10, pady=10, sticky="w")
alumina_activ_label = ttk.Label(alumina_activ_frame, text="activ: ")
alumina_activ_label.grid(row=0, column=0, padx=5, pady=5)

bath_cond_frame = ttk.LabelFrame(root, text="bath conductivity")
bath_cond_frame.grid(row=4, column=bath_col, padx=10, pady=10, sticky="w")
bath_cond_label = ttk.Label(bath_cond_frame, text="bath cond: ")
bath_cond_label.grid(row=0, column=0, padx=5, pady=5)

#Cell input frame GUI

cell_frame = ttk.LabelFrame(root, text="Cell input data")
cell_frame.grid(row=0, column=0, padx=5, pady=10, sticky="w")

ttk.Label(cell_frame, text="Current [A]", width=10).grid(row=0, column=0, padx=5, pady=5)
cell_current_entry = ttk.Entry(cell_frame, textvariable=I_line, width=6)
cell_current_entry.grid(row=0, column=1, padx=1, pady=1)

ttk.Label(cell_frame, text="ACD [cm]", width=10).grid(row=1, column=0, padx=5, pady=5)
ACD_entry = ttk.Entry(cell_frame, textvariable=ACD, width=6)
ACD_entry.grid(row=1, column=1, padx=1, pady=1)

#Output
V_cell_output_col = 6
V_cell_frame = ttk.LabelFrame(root, text="V_cell [V]")
V_cell_frame.grid(row=0, column=V_cell_output_col, padx=10, pady=10, sticky="w")

V_cell_label = ttk.Label(V_cell_frame, text="V_cell: ")
V_cell_label.grid(row=0, column=V_cell_output_col, padx=5, pady=5)

E_rev_label = ttk.Label(V_cell_frame, text="E_rev: ")
E_rev_label.grid(row=1, column=V_cell_output_col, padx=5, pady=5)

E_ca_label = ttk.Label(V_cell_frame, text="E_ca: ")
E_ca_label.grid(row=2, column=V_cell_output_col, padx=5, pady=5)

E_sa_label = ttk.Label(V_cell_frame, text="E_sa: ")
E_sa_label.grid(row=3, column=V_cell_output_col, padx=5, pady=5)

E_cc_label = ttk.Label(V_cell_frame, text="E_cc: ")
E_cc_label.grid(row=4, column=V_cell_output_col, padx=5, pady=5)

E_cc_label = ttk.Label(V_cell_frame, text="E_cc: ")
E_cc_label.grid(row=4, column=V_cell_output_col, padx=5, pady=5)

Vbath_label = ttk.Label(V_cell_frame, text="Vbath: ")
Vbath_label.grid(row=5, column=V_cell_output_col, padx=5, pady=5)

Vbub_label = ttk.Label(V_cell_frame, text="Vbub: ")
Vbub_label.grid(row=6, column=V_cell_output_col, padx=5, pady=5)

V_an_label = ttk.Label(V_cell_frame, text="V_an: ")
V_an_label.grid(row=7, column=V_cell_output_col, padx=5, pady=5)

V_ca_label = ttk.Label(V_cell_frame, text="V_ca: ")
V_ca_label.grid(row=8, column=V_cell_output_col, padx=5, pady=5)

V_ext_label = ttk.Label(V_cell_frame, text="V_ext: ")
V_ext_label.grid(row=9, column=V_cell_output_col, padx=5, pady=5)


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

    #Bath ratio calc to display in GUI
    bath_ratio = Bath.Ratio(StateVector[0], ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2])[0]
    cryolite_ratio_label.config(text=f"bath_ratio: {bath_ratio:.4f}")

    #Alumina saturation
    alumina_sat = Bath.AluminaSaturation(StateVector[0], ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2], ConstVector[4])[0]
    alumina_sat_label.config(text=f"Alumina sat: {alumina_sat:.4f}")

    #Alumina activity
    alumina_actv = Bath.AluminaActivity(StateVector[0], ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2], ConstVector[4])[0]
    alumina_activ_label.config(text=f"Alumina actv: {alumina_actv:.4f}")

    #Bath conductivity
    bath_cond = Bath.Conductivity(StateVector[0], ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2],
                                        ConstVector[4])[0]
    bath_cond_label.config(text=f"Bath cond: {bath_cond:.4f}")

    #Vcell total volts
    Vcell_total = Volt.Vcell(StateVector, InputVector, ConstVector)[0]
    V_cell_label.config(text=f"V_cell: {Vcell_total:.4f}")
    print(Vcell_total)

    #Erev volts
    E_rev = Electro.ReversiblePotential(StateVector[0], ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2],
                                        ConstVector[4])[0]
    E_rev_label.config(text=f"E_rev: {E_rev:.4f}")

    E_ca = Electro.AnodeConcOverVolt(StateVector[0], ConstVector[4], InputVector[0])[0]
    E_ca_label.config(text=f"E_ca: {E_ca:.4f}")

    E_sa = Electro.AnodeSurfOverVolt(StateVector[0], ConstVector[4], InputVector[0])[0]
    E_sa_label.config(text=f"E_sa: {E_sa:.4f}")

    E_cc = Electro.CathodeConcOverVolt(bath_ratio, ConstVector[4], InputVector[0])[0]
    E_cc_label.config(text=f"E_cc: {E_cc:.4f}")

    dB = Electro.BubbleThickness(InputVector[0])
    kbath = Bath.Conductivity(StateVector[0], ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2], ConstVector[4])[0]
    Rbath = Electro.BathRes(StateVector[2], dB, kbath)
    Vbath = InputVector[0]*Rbath
    Vbath_label.config(text=f"Vbath: {Vbath:.4f}")

    BubbleCoverage = Electro.BubbleCoverage(StateVector[0],bath_ratio,InputVector[0])[0]
    BubbleRes = Electro.BubbleRes(BubbleCoverage,dB,kbath)
    Vbub = InputVector[0]*BubbleRes
    Vbub_label.config(text=f"Vbub: {Vbub:.4f}")

    V_an = InputVector[0]*Glob.Ran
    V_an_label.config(text=f"V_an: {V_an:.4f}")

    V_ca = InputVector[0]*Glob.Rca
    V_ca_label.config(text=f"V_ca: {V_ca:.4f}")

    V_ext = InputVector[0]*Glob.Rext
    V_ext_label.config(text=f"V_ext: {V_ext:.4f}")


def plot_cell_volt():

    # Define a range of C_Al2O3 values from 2 to 8
    C_Al2O3_values = np.linspace(1.9, 8, 100)
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



    for C_Al2O3 in C_Al2O3_values:
        BR = Bath.Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2)[0]
        Erev=Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]
        Erev_values.append(Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0])

        Esa=Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0]
        Esa_values.append(Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0])
        BCover = Electro.BubbleCoverage(C_Al2O3, BR, I_line)[0]

        Vbub=I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line),
                                                      Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[
                                                          0])
        Vbub_values.append(I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line),
                                                      Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[
                                                          0]))

        Eca=Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0]

        Eca_values.append(Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0])

        Ecc=Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0]
        Ecc_values.append(Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0])

        dB = Electro.BubbleThickness(I_line)
        kbath = Bath.Conductivity(C_Al2O3, ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2],
                                  ConstVector[4])[0]
        Rbath = Electro.BathRes(ACD, dB, kbath)

        Vbath=I_line * Rbath
        Vbath_values.append(I_line * Rbath)

        V_ca=I_line * Glob.Ran
        V_an = I_line * Glob.Rca
        V_ext = I_line * Glob.Rext

        V_ca_values.append(I_line * Glob.Ran)
        V_an_values.append(I_line * Glob.Rca)
        V_ext_values.append(I_line * Glob.Rext)

        V_cell_values.append(Erev+Esa+Vbub+Eca+Ecc+Vbath+V_ca+V_an+V_ext)

    # Create a figure and axis for the plot
    # fig = Figure(figsize=(8, 6))
    # ax = fig.add_subplot(111)
    fig, (ax2, ax1) = plt.subplots(2, 1, figsize=(8, 8), sharex=False)

    #plt.subplots_adjust(right=0.6, bottom=0.2)
    # Plot the calculated values against C_Al2O3
    # ax.plot(C_Al2O3_values, Erev_values, label="E_rev")
    # ax.plot(C_Al2O3_values, Esa_values, label="E_sa")
    # ax.plot(C_Al2O3_values, Vbub_values, label="V_bub")
    # ax.plot(C_Al2O3_values, Eca_values, label="E_ca")
    # ax.plot(C_Al2O3_values, Ecc_values, label="E_cc")
    # ax.plot(C_Al2O3_values, Vbath_values, label="V_bath")
    # ax.plot(C_Al2O3_values, Vbub_values, label="V_bub")
    # ax.plot(C_Al2O3_values, V_ca_values, label="V_ca")
    # ax.plot(C_Al2O3_values, V_an_values, label="V_an")
    # ax.plot(C_Al2O3_values, V_ext_values, label="V_ext")

    # Plot the first set of data on the first subplot (ax1)
    ax1.plot(C_Al2O3_values, Erev_values, label="E_rev")
    ax1.plot(C_Al2O3_values, Esa_values, label="E_sa")
    ax1.plot(C_Al2O3_values, Eca_values, label="E_ca")
    ax1.plot(C_Al2O3_values, Ecc_values, label="E_cc")
    ax1.plot(C_Al2O3_values, Vbub_values, label="V_bub")
    ax1.plot(C_Al2O3_values, Vbath_values, label="V_bath")
    #ax1.plot(C_Al2O3_values, Vbub_values, label="V_bub")
    ax1.plot(C_Al2O3_values, V_ca_values, label="V_ca")
    ax1.plot(C_Al2O3_values, V_an_values, label="V_an")
    ax1.plot(C_Al2O3_values, V_ext_values, label="V_ext")
    ax1.set_xlabel("Al2O3 %wt.")
    ax1.set_ylabel("Components of cell voltage [V]")
    ax1.legend()

    # Plot the second set of data on the second subplot (ax2)
    ax2.plot(C_Al2O3_values, V_cell_values, label="Total cell voltage")
    ax2.set_xlabel("Al2O3 %wt.")
    ax2.set_ylabel("Total cell voltage [V]")
    ax2.legend()
    #
    # ax.plot(C_Al2O3_values, V_cell_values, label="V_cell")
    #
    #
    #
    # # Set labels and legend
    # ax.set_xlabel("C_Al2O3")
    # ax.set_ylabel("Voltage")
    # ax.legend()
    # Adjust spacing between the subplots to avoid overlapping labels
    plt.tight_layout()

    return fig

def plot_cell_volt_dynamic():



    #C_Al2O3_values = np.linspace(1.9, 8, 100)
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
    time = np.linspace(1000, 1500, 501)
    C_Al2O3_values = []
    for t in time:
        C_Al2O3 = 7.5977e-18 * pow(t, 6) - 3.8081e-14 * pow(t, 5) + 7.1245e-11 * pow(t, 4)\
                  - 6.0775e-08 * pow(t, 3) + 2.2483e-05 * pow(t, 2) - 1.2591e-03 * pow(t, 1) + 7.2423e-01
        C_Al2O3_values.append(C_Al2O3)
        BR = Bath.Ratio(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2)[0]
        Erev=Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0]
        Erev_values.append(Electro.ReversiblePotential(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[0])

        Esa=Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0]
        Esa_values.append(Electro.AnodeSurfOverVolt(C_Al2O3, T_Bath, I_line)[0])
        BCover = Electro.BubbleCoverage(C_Al2O3, BR, I_line)[0]

        Vbub=I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line),
                                                      Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[
                                                          0])
        Vbub_values.append(I_line * Electro.BubbleRes(BCover, Electro.BubbleThickness(I_line),
                                                      Bath.Conductivity(C_Al2O3, C_AlF3, C_CaF2, C_LiF, C_MgF2, T_Bath)[
                                                          0]))

        Eca=Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0]

        Eca_values.append(Electro.AnodeConcOverVolt(C_Al2O3, T_Bath, I_line)[0])

        Ecc=Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0]
        Ecc_values.append(Electro.CathodeConcOverVolt(BR, T_Bath, I_line)[0])

        dB = Electro.BubbleThickness(I_line)
        kbath = Bath.Conductivity(C_Al2O3, ConstVector[0], ConstVector[1], ConstVector[3], ConstVector[2],
                                  ConstVector[4])[0]
        Rbath = Electro.BathRes(ACD, dB, kbath)

        Vbath=I_line * Rbath
        Vbath_values.append(I_line * Rbath)

        V_ca=I_line * Glob.Ran
        V_an = I_line * Glob.Rca
        V_ext = I_line * Glob.Rext

        V_ca_values.append(I_line * Glob.Ran)
        V_an_values.append(I_line * Glob.Rca)
        V_ext_values.append(I_line * Glob.Rext)

        V_cell_values.append(Erev+Esa+Vbub+Eca+Ecc+Vbath+V_ca+V_an+V_ext)

    # Create a figure and axis for the plot
    # fig = Figure(figsize=(8, 6))
    # ax = fig.add_subplot(111)
    fig, (ax2, ax1) = plt.subplots(2, 1, figsize=(8, 8), sharex=False)

    #plt.subplots_adjust(right=0.6, bottom=0.2)
    # Plot the calculated values against C_Al2O3
    # ax.plot(C_Al2O3_values, Erev_values, label="E_rev")
    # ax.plot(C_Al2O3_values, Esa_values, label="E_sa")
    # ax.plot(C_Al2O3_values, Vbub_values, label="V_bub")
    # ax.plot(C_Al2O3_values, Eca_values, label="E_ca")
    # ax.plot(C_Al2O3_values, Ecc_values, label="E_cc")
    # ax.plot(C_Al2O3_values, Vbath_values, label="V_bath")
    # ax.plot(C_Al2O3_values, Vbub_values, label="V_bub")
    # ax.plot(C_Al2O3_values, V_ca_values, label="V_ca")
    # ax.plot(C_Al2O3_values, V_an_values, label="V_an")
    # ax.plot(C_Al2O3_values, V_ext_values, label="V_ext")

    # Plot the first set of data on the first subplot (ax1)
    ax1.plot(time, Erev_values, label="E_rev")
    ax1.plot(time, Esa_values, label="E_sa")
    ax1.plot(time, Eca_values, label="E_ca")
    ax1.plot(time, Ecc_values, label="E_cc")
    ax1.plot(time, Vbub_values, label="V_bub")
    ax1.plot(time, Vbath_values, label="V_bath")
    #ax1.plot(time, Vbub_values, label="V_bub")
    ax1.plot(time, V_ca_values, label="V_ca")
    ax1.plot(time, V_an_values, label="V_an")
    ax1.plot(time, V_ext_values, label="V_ext")
    ax1.set_xlabel("time [s]")
    ax1.set_ylabel("Components of cell voltage [V]")
    ax1.legend()

    # Plot the second set of data on the second subplot (ax2)
    ax2.plot(time, V_cell_values, label="Total cell voltage")
    ax2.set_xlabel("time [s]")
    ax2.set_ylabel("Total cell voltage [V]")
    ax2.legend()
    #
    # ax.plot(C_Al2O3_values, V_cell_values, label="V_cell")
    #
    #
    #
    # # Set labels and legend
    # ax.set_xlabel("C_Al2O3")
    # ax.set_ylabel("Voltage")
    # ax.legend()
    # Adjust spacing between the subplots to avoid overlapping labels
    plt.tight_layout()

    return fig

def plot_button_cmd():
    #    # Call the plot_cell_volt function to get the figure
    #fig = plot_cell_volt()
    fig = plot_cell_volt_dynamic()
    # Set the size of the Figure
    #fig.set_size_inches(10, 4)
    plt.subplots_adjust(right=0.6, bottom=0.2)

    # Create a new frame for the plot
    plot_frame = ttk.LabelFrame(root, text="Plot", width=300, height=300)
    plot_frame.grid(row=0, column=8, columnspan=2, pady=10)
    #plot_frame.pack()

    # Create a canvas for the plot
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    # Display the plot in the new frame
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    toolbar = NavigationToolbar2Tk(canvas, plot_frame)
    toolbar.update()
    #canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)




update_button = ttk.Button(root, text="Update Values", command=update_values)
update_button.grid(row=1, column=0, columnspan=4, pady=10)

# Modify the command for the "Plot" button to call the plot_button_cmd function
plot_button = ttk.Button(root, text="Plot", command=plot_button_cmd)
plot_button.grid(row=1, column=1, columnspan=5, pady=10)

root.mainloop()
