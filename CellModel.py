import Global as Glob
import BathProperties as Bath
import ElectroChemistry as Electro
import pandas as pd
"""
Vcell = Erev + Eca + Esa + Ecc + I_line*Rpath + I_line*Glob.Rext
"""

C_Al2O3 = 4.2
C_AlF3 =10.3
C_CaF2 =7
C_LiF =0
C_MgF2 =0.1
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
# anodebottomarea = Glob.AAnode/Glob.nAnode
# ibottomcd = ibottom/anodebottomarea
# Vbath = ibottomcd*(ACD-db)*(1/kbath)
Vdroppath = Van+Vbath+Vbub+Vca

Vcell = Erev + Eca + Esa + Ecc + Vdroppath + Vdropext
Volt_table = pd.DataFrame([{'Erev': Erev, 'Eca': Eca, 'Esa': Esa, 'Ecc': Ecc, 'Vdroppath': Vdroppath, 'Vdropext': Vdropext}])
print(f"Cell voltage: {Volt_table.sum(axis=1)[0]}")
