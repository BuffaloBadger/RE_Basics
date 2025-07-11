import pandas as pd
from scipy import optimize

# Given or known
P = 1 # atm
T = 150 + 273.15 # K
y_C2H6_0 = 0.9
y_O2_0 = 0.21 * 0.1
f_O2 = 0.5
s_C2H4_CO2 = 3.0

# Basis
n_total_0 = 1.0 # mol

# Initial mole fractions
y_C2H4_0 = 0 # eqn 7
y_CO2_0 = 0 # eqn 8

# Initial moles
n_O2_0 = y_O2_0 * n_total_0 # eqn 4
n_C2H4_0 = y_C2H4_0 * n_total_0 # eqn 5
n_CO2_0 = y_CO2_0 * n_total_0 # eqn 6

# Solve eqns 3 and 4
# equations to solve, as residuals
def residuals(x):
    return [f_O2 * n_O2_0 - x[0] - 7 * x[1],
        s_C2H4_CO2 * (n_CO2_0 + 4 * x[1]) - n_C2H4_0 - 2 * x[0]]
# guess for solution
guess = [0.5, 0.5]
# try to solve
soln = optimize.root(residuals, guess)
# check for success
if not(soln.success):
    print("A solution was NOT obtained:")
    print(soln.message)
# extract the solution
extent = soln.x

# Calculate the mole fraction of CO2
y_CO2 = (n_CO2_0 + 4 * extent[1]) / (n_total_0 + extent[0] + extent[1])
y_CO2 = float('%.3g' % y_CO2)

# Display the result
print(' ')
print('CO2 Mole Fraction:',y_CO2)
print(' ')

# Save the result to a .csv file
data = [['y_CO2', y_CO2]]
result = pd.DataFrame(data, columns=['item','value'])
result.to_csv("reb_3_4_1_results.csv", index=False)
