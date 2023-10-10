import ElectroChemistry
import BathProperties
import Global as Glob
import numpy as np
import pandas as pd
from sensitivity import SensitivityAnalyzer
#ElectroChemistry.BubbleThickness(100000)
I_cell = 126000
T_bath_step = 20
Report_structure = {'Cell component': ['Equilibrium potential Ee', 'Anode conc. overvoltage',
                                       'Anode surf. overvoltage', 'Cathode overvoltage', 'Bubbles voltage drop',
                                       'Bath voltage drop'], 'Min': [0,0,0,0,0,0], 'Max': [0,0,0,0,0,0], 'Range': [0,0,0,0,0,0]}
Report = pd.DataFrame(Report_structure)
#reversible potential sensitivity analysis
Erev_sensitivity_dict = {
    'C_Al2O3' : (np.arange(1, 10, 1)).tolist(),
    'C_AlF3' : (np.arange(0, 12, 1)).tolist(),
    'C_CaF2' : (np.arange(0, 5, 1)).tolist(),
    'C_LiF' : (np.arange(0, 5, 1)).tolist(),
    'C_MgF2' : (np.arange(0, 5, 1)).tolist(),
    'T_Bath' : (np.arange(900, 1000, T_bath_step)).tolist()
}
Erev_labels = {
    'C_Al2O3' : 'C_Al2O3',
    'C_AlF3' : 'C_AlF3',
    'C_CaF2' : 'C_CaF2',
    'C_LiF' : 'C_LiF',
    'C_MgF2' : 'C_MgF2',
    'T_Bath' : 'T_Bath'
}
Erev_sa = SensitivityAnalyzer(sensitivity_values=Erev_sensitivity_dict, func=ElectroChemistry.ReversiblePotential, grid_size=3, labels=Erev_labels)
Erev_plot = Erev_sa.plot()

Erev_results=pd.DataFrame()
Erev_max=Erev_sa.df[Erev_sa.df.Result == Erev_sa.df.Result.max()]
Erev_min=Erev_sa.df[Erev_sa.df.Result == Erev_sa.df.Result.min()]
Erev_max = Erev_max.reset_index(drop=True)
Erev_min = Erev_min.reset_index(drop=True)
Erev_max_value = Erev_max['Result'][0][0]
Erev_min_value = Erev_min['Result'][0][0]
Erev_range_diff = Erev_max_value - Erev_min_value

Report.loc[0, 'Min'] = Erev_min_value
Report.loc[0, 'Max'] = Erev_max_value
Report.loc[0, 'Range'] = Erev_range_diff


Erev_results=pd.concat([Erev_results, Erev_max, Erev_min])

###Anode concentration overvoltage sensitivity
Eca_sensitivity_dict = {
    'C_Al2O3' : (np.arange(1, 10, 1)).tolist(),
    'T_Bath' : (np.arange(900, 1000, T_bath_step)).tolist(),
    'I_line' : (np.arange(115000, 135000, 10000)).tolist()
}

Eca_labels = {
    'C_Al2O3' : 'C_Al2O3',
    'T_Bath' : 'T_Bath',
    'I_line' : 'I_line'
}
Eca_sa = SensitivityAnalyzer(sensitivity_values=Eca_sensitivity_dict, func=ElectroChemistry.AnodeConcOverVolt, grid_size=3, labels=Eca_labels)
Eca_plot = Eca_sa.plot()
####SUMMARY
Eca_results=pd.DataFrame()
Eca_max=Eca_sa.df[Eca_sa.df.Result == Eca_sa.df.Result.max()]
Eca_min=Eca_sa.df[Eca_sa.df.Result == Eca_sa.df.Result.min()]
Eca_max = Eca_max.reset_index(drop=True)
Eca_min = Eca_min.reset_index(drop=True)
Eca_max_value = Eca_max['Result'][0][0]
Eca_min_value = Eca_min['Result'][0][0]
Eca_range_diff = Eca_max_value - Eca_min_value
Report.loc[1, 'Min'] = Eca_min_value
Report.loc[1, 'Max'] = Eca_max_value
Report.loc[1, 'Range'] = Eca_range_diff
Eca_results=pd.concat([Eca_results, Eca_max, Eca_min])

####
###Anode surface overvoltage sensitivity
Esa_sensitivity_dict = {
    'C_Al2O3' : (np.arange(1, 10, 1)).tolist(),
    'T_Bath' : (np.arange(900, 1000, T_bath_step)).tolist(),
    'I_line' : (np.arange(115000, 135000, 10000)).tolist()
}

Esa_labels = {
    'C_Al2O3' : 'C_Al2O3',
    'T_Bath' : 'T_Bath',
    'I_line' : 'I_line'
}
Esa_sa = SensitivityAnalyzer(sensitivity_values=Esa_sensitivity_dict, func=ElectroChemistry.AnodeSurfOverVolt, grid_size=3, labels=Esa_labels)
Esa_plot = Esa_sa.plot()
####SUMMARY
Esa_results=pd.DataFrame()
Esa_max=Esa_sa.df[Esa_sa.df.Result == Esa_sa.df.Result.max()]
Esa_min=Esa_sa.df[Esa_sa.df.Result == Esa_sa.df.Result.min()]
Esa_max = Esa_max.reset_index(drop=True)
Esa_min = Esa_min.reset_index(drop=True)
Esa_max_value = Esa_max['Result'][0][0]
Esa_min_value = Esa_min['Result'][0][0]
Esa_range_diff = Esa_max_value - Esa_min_value
Report.loc[2, 'Min'] = Esa_min_value
Report.loc[2, 'Max'] = Esa_max_value
Report.loc[2, 'Range'] = Esa_range_diff
Esa_results=pd.concat([Esa_results, Esa_max, Esa_min])




###Cathode surface overvoltage sensitivity
Ecc_sensitivity_dict = {
    'Ratio' : (np.arange(0.9, 1.6, 0.1)).tolist(),
    'T_Bath' : (np.arange(900, 1000, T_bath_step)).tolist(),
    'I_line' : (np.arange(115000, 135000, 10000)).tolist()
}

Ecc_labels = {
    'Ratio' : 'Ratio',
    'T_Bath' : 'T_Bath',
    'I_line' : 'I_line'
}
Ecc_sa = SensitivityAnalyzer(sensitivity_values=Ecc_sensitivity_dict, func=ElectroChemistry.CathodeConcOverVolt, grid_size=3, labels=Ecc_labels)
Ecc_plot = Ecc_sa.plot()
####SUMMARY
Ecc_results=pd.DataFrame()
Ecc_max=Ecc_sa.df[Ecc_sa.df.Result == Ecc_sa.df.Result.max()]
Ecc_min=Ecc_sa.df[Ecc_sa.df.Result == Ecc_sa.df.Result.min()]
Ecc_max = Ecc_max.reset_index(drop=True)
Ecc_min = Ecc_min.reset_index(drop=True)
Ecc_max_value = Ecc_max['Result'][0][0]
Ecc_min_value = Ecc_min['Result'][0][0]
Ecc_range_diff = Ecc_max_value - Ecc_min_value
Report.loc[3, 'Min'] = Ecc_min_value
Report.loc[3, 'Max'] = Ecc_max_value
Report.loc[3, 'Range'] = Ecc_range_diff
Ecc_results=pd.concat([Ecc_results, Ecc_max, Ecc_min])
#bath conductivity sensitivity analysis

kbath_sensitivity_dict = {
    'C_Al2O3' : (np.arange(1, 10, 1)).tolist(),
    'C_AlF3' : (np.arange(0, 12, 1)).tolist(),
    'C_CaF2' : (np.arange(0, 5, 1)).tolist(),
    'C_LiF' : (np.arange(0, 5, 1)).tolist(),
    'C_MgF2' : (np.arange(0, 5, 1)).tolist(),
    'T_Bath' : (np.arange(900, 1000, T_bath_step)).tolist()
}
kbath_labels = {
    'C_Al2O3' : 'C_Al2O3',
    'C_AlF3' : 'C_AlF3',
    'C_CaF2' : 'C_CaF2',
    'C_LiF' : 'C_LiF',
    'C_MgF2' : 'C_MgF2',
    'T_Bath' : 'T_Bath'
}
kbath_sa = SensitivityAnalyzer(sensitivity_values=kbath_sensitivity_dict, func=BathProperties.Conductivity, grid_size=6, labels=kbath_labels)
kbath_plot = kbath_sa.plot()
####SUMMARY
kbath_results=pd.DataFrame()
kbath_max=kbath_sa.df[kbath_sa.df.Result == kbath_sa.df.Result.max()]
kbath_min=kbath_sa.df[kbath_sa.df.Result == kbath_sa.df.Result.min()]
kbath_results=pd.concat([kbath_results, kbath_max, kbath_min])
kbath_plot.show()
#####Coverage

Coverage_sensitivity_dict = {
    'C_Al2O3': (np.arange(1, 10, 1)).tolist(),
    'Ratio': (np.arange(0.9, 1.6, 0.1)).tolist(),
    'I_line': (np.arange(115000, 135000, 10000)).tolist()

}
Coverage_labels = {
    'C_Al2O3': 'C_Al2O3',
    'Ratio': 'Ratio',
    'I_line': 'I_line'
}
Coverage_sa = SensitivityAnalyzer(sensitivity_values=Coverage_sensitivity_dict, func=ElectroChemistry.BubbleCoverage,
                                  grid_size=3, labels=Coverage_labels)
Coverage_plot = Coverage_sa.plot()
####SUMMARY
Coverage_results = pd.DataFrame()
Coverage_max = Coverage_sa.df[Coverage_sa.df.Result == Coverage_sa.df.Result.max()]
Coverage_min = Coverage_sa.df[Coverage_sa.df.Result == Coverage_sa.df.Result.min()]
Coverage_results = pd.concat([Coverage_results, Coverage_max, Coverage_min])




##bubble thickness sensitivity analysis
db=[]
I_line = (np.arange(115000, 135000, 5000)).tolist()
for I_line in I_line:
    db.append(ElectroChemistry.BubbleThickness(I_line))
db_min=min(db)
db_max=max(db)
# db_sensitivity_dict = {
#     'I_line' : (np.arange(115000, 135000, 5000)).tolist(),
#     'f' : (np.arange(1.05, 1.15, 0.01)).tolist()
#
# }
# db_labels = {
#     'I_line' : 'I_line',
#     'f' : 'f'
# }
# db_sa = SensitivityAnalyzer(sensitivity_values=db_sensitivity_dict, func=ElectroChemistry.BubbleThickness, grid_size=3, labels=db_labels)
# db_plot = db_sa.plot()
# ####SUMMARY
# db_results=pd.DataFrame()
# db_max=db_sa.df[db_sa.df.Result == db_sa.df.Result.max()]
# db_min=db_sa.df[db_sa.df.Result == db_sa.df.Result.min()]
# db_results=pd.concat([db_results, db_max, db_min])


##Bubble resistance sensitivity analysis
R_bub_sensitivity_dict = {
    'Coverage' : (np.arange(0.46e0, 0.9e0, 0.1e0)).tolist(),
    'db' : (np.arange(0.502e0, 0.504e0, 0.00025e0)).tolist(),
    'k_bath' : (np.arange(1.76e0, 2.92e0, 0.1e0)).tolist()
}
R_bub_labels = {
    'Coverage' : 'Coverage',
    'db' : 'db',
    'k_bath' : 'k_bath',

}
R_bub_sa = SensitivityAnalyzer(sensitivity_values=R_bub_sensitivity_dict, func=ElectroChemistry.BubbleRes, grid_size=3, labels=R_bub_labels)
R_bub_plot = R_bub_sa.plot()

####SUMMARY
R_bub_results=pd.DataFrame()
R_bub_max=R_bub_sa.df[R_bub_sa.df.Result == R_bub_sa.df.Result.max()]
R_bub_min=R_bub_sa.df[R_bub_sa.df.Result == R_bub_sa.df.Result.min()]
R_bub_results=pd.concat([R_bub_results, R_bub_max, R_bub_min])
Vbub_min=min(R_bub_results['Result'].to_numpy())*115000
Vbub_max=max(R_bub_results['Result'].to_numpy())*135000

Vbub_range_diff = Vbub_max - Vbub_min
Report.loc[4, 'Min'] = Vbub_min
Report.loc[4, 'Max'] = Vbub_max
Report.loc[4, 'Range'] = Vbub_range_diff

##Bath resistance sensitivity analysis
Rbath_sensitivity_dict = {
    'ACD' : (np.arange(2e0, 4e0, 0.2e0)).tolist(),
    'db' : (np.arange(0.2e0, 2e0, 0.2e0)).tolist(),
    'k_bath' : (np.arange(1.76e0, 2.92e0, 0.1e01)).tolist()
}
Rbath_labels = {
    'ACD' : 'ACD',
    'db' : 'db',
    'k_bath' : 'k_bath',

}
Rbath_sa = SensitivityAnalyzer(sensitivity_values=Rbath_sensitivity_dict, func=ElectroChemistry.BathRes, grid_size=3, labels=Rbath_labels)
Rbath_plot = Rbath_sa.plot()

####SUMMARY
Rbath_results=pd.DataFrame()
Rbath_max=Rbath_sa.df[Rbath_sa.df.Result == Rbath_sa.df.Result.max()]
Rbath_min=Rbath_sa.df[Rbath_sa.df.Result == Rbath_sa.df.Result.min()]
Rbath_max_value = Rbath_max['Result'].item()
Rbath_min_value = Rbath_min['Result'].item()
Rbath_results=pd.concat([Rbath_results, Rbath_max, Rbath_min])
# Rbath_results["Result"] = Rbath_results["Result"].map("{:,.4e}".format)
# Rbath_results=Rbath_results.reset_index(drop=True)
#Vbath=Rbath_results['Result'].to_numpy()*I_cell
Vbath_min = Rbath_min_value*115000
Vbath_max = Rbath_max_value*135000

Vbath_range_diff = Vbath_max - Vbath_min
Report.loc[5, 'Min'] = Vbath_min
Report.loc[5, 'Max'] = Vbath_max
Report.loc[5, 'Range'] = Vbath_range_diff