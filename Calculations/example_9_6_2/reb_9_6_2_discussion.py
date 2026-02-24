"""Calculations for Reaction Engineering Basics Example 9.6.2 Discussion"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import pandas as pd
import reb_9_6_2

# define local constants
CA_0 = reb_9_6_2.CA_0
CB_0 = reb_9_6_2.CB_0
V = reb_9_6_2.V
Tex0 = reb_9_6_2.Tex0
fA_f = reb_9_6_2.fA_f

# set a lower initial temperature
T_0 = 55 + 273.15 # K

ind_0 = 0.0
dep_0 = np.array([CA_0*V, CB_0*V, 0.0, 0.0, 0.0, T_0, Tex0])

# define the stopping criterion
f_var = 1
f_val = CA_0*V*(1-fA_f)
    
# solve the IVODEs
t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , reb_9_6_2.derivatives)

# check that a solution was found
if not(success):
    print(f"An IVODE solution was NOT obtained: {message}")

nX = dep[2,:]
nZ = dep[4,:]
sel_X_Z = nX[-1]/nZ[-1]

# read in the results from the assignment
results_df = pd.read_csv('results.csv')

# add the new results
n_rows = len(results_df.index)
results_df.loc[n_rows] = ['lower_T',f'{T_0 - 273.15}','°C']
results_df.loc[n_rows+1] = ['sel_lower_T',f'{sel_X_Z}','mol X per mol Z']
results_df.loc[n_rows+2] = ['t_lower_T',f'{t[-1]}', 'min']

# display the results
print(' ')
print(f'Results with an initial temperature of {T_0 - 273.15:.3g} °C')
print(f'Selectivity: {sel_X_Z:.3g} mol X per mol Z')
print(f'Reaction Time: {t[-1]:.3g} min')

# save the results
results_df.to_csv('results.csv', index=False)
