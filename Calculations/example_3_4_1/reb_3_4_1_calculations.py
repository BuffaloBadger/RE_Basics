""" Calculations for REB Example 3.4.1 """

import pandas as pd
from scipy import optimize

# Global constants
P = 1 # atm
T = 150 + 273.15 # K
y_C2H6_0 = 0.9
y_air_0 = 0.1
f_O2 = 0.5
S_C2H4_CO2 = 3.0

# Basis
n_total_0 = 1.0 # mol

# Calculation of apparent extents of reaction
n_O2_0 = 0.21*y_air_0*n_total_0
# equations to solve, as residuals
def residuals(x):
    return [f_O2 * n_O2_0 - x[0] - 7 * x[1],
        S_C2H4_CO2 * (4 * x[1]) - 2 * x[0]]
# guess for solution
guess = [0.5, 0.5]
# try to solve
soln = optimize.root(residuals, guess)
# check for success
if not(soln.success):
    print("A solution was NOT obtained:")
    print(soln.message)
# extract the solution
extent1 = soln.x[0]
extent2 = soln.x[1]

# Calculation of final CO2 mole fraction
y_CO2_f = 4*extent2/(n_total_0 + extent1 + extent2)

# Display the result
print(' ')
print('CO2 Mole Fraction:',y_CO2_f)
print(' ')

# Save the results to a .csv file
data = [['extent 1', extent1, 'mol'],
        ['extent 2', extent2, 'mol'],
        ['y_CO2_final', y_CO2_f, '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv("reb_3_4_1_results.csv", index=False)
