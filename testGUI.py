import tkinter as tk
from tkinter import ttk

##Initialize vars

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

# def update_variable(value, idx):
#     slider_vars[idx].set(float(value))
def update_values():
    # Update StateVector and ConstVector with the current slider values
    StateVector[0] = C_Al2O3.get()
    ConstVector[0] = C_AlF3.get()
    ConstVector[1] = C_CaF2.get()
    ConstVector[2] = C_MgF2.get()
    ConstVector[3] = C_LiF.get()
    ConstVector[4] = T_Bath.get()

update_button = ttk.Button(root, text="Update Values", command=update_values)
update_button.grid(row=len(slider_labels), column=0, columnspan=2, pady=10)

root.mainloop()
