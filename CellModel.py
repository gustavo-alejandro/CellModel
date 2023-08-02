import CellVoltage as Volt
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

Vcell = Volt.Vcell(StateVector, InputVector, ConstVector)[0]