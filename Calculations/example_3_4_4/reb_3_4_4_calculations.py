import pandas as pd
from scipy import optimize

# Given or known
n_CO_0 = 34.0 # mol
n_H2_0 = 66.0 # mol
n_CH3OH_0 = 0.0
n_CH4_0 = 0.0
n_H2O_0 = 0.0
n_CO2_0 = 0.0
f_CO = 0.4
Y_CH3OH = 0.107
S_CH3OH_CH4 = 0.38

# Solve eqns 7 through 9 for the apparent extents of reaction
# equations to solve, as residuals
def residuals(x):
    return [x[0] + x[1] + x[2] - n_CO_0*f_CO,
        n_CH3OH_0 + x[0] - Y_CH3OH*n_CO_0,
        n_CH3OH_0 + x[0] - S_CH3OH_CH4*(n_CH4_0 + x[1])]
# guess for solution
guess = [10, 10, 10]
# try to solve
soln = optimize.root(residuals, guess)
# check for success
if not(soln.success):
    print("A solution was NOT obtained:")
    print(soln.message)
# extract the solution
extent = soln.x

# Calculate the final moles using eqns 10 through 15
n_CO = n_CO_0 - extent[0] - extent[1] - extent[2]
n_H2 = n_H2_0 - 2*extent[0] - 3*extent[1] + extent[2]
n_CH3OH = n_CH3OH_0 + extent[0]
n_CH4 = n_CH4_0 + extent[1]
n_H2O = n_H2O_0 + extent[1] - extent[2]
n_CO2 = n_CO2_0 + extent[2]

# Check the results by calculating the conversion, yield, and selectivity
X_CO_check = (n_CO_0 - n_CO)/n_CO_0
Y_CH3OH_check = n_CH3OH/n_CO_0
S_CH3OH_CH4_check = n_CH3OH/(n_CH4)

# Display the result
print(' ')
print('Final Moles')
print(' CO:','%.3g' % n_CO)
print(' H2:','%.3g' % n_H2)
print(' CH3OH:','%.3g' % n_CH3OH)
print(' CH4:','%.3g' % n_CH4)
print(' H2O:','%.3g' % n_H2O)
print(' CO2:','%.3g' % n_CO2)
print(' ')
print('Check')
print(' X_CO:','%.3g' % X_CO_check)
print(' Y_CH3OH:','%.3g' % Y_CH3OH_check)
print(' S_CH3OH_CH4:','%.3g' % S_CH3OH_CH4_check)

# Save the result to a .csv file
data = [['n_CO', '%.3g' % n_CO, 'mol'], \
    ['n_H2', '%.3g' % n_H2, 'mol'], \
    ['n_CH3OH', '%.3g' % n_CH3OH, 'mol'], \
    ['n_CH4', '%.3g' % n_CH4, 'mol'], \
    ['n_H2O', '%.3g' % n_H2O, 'mol'], \
    ['n_CO2', '%.3g' % n_CO2, 'mol'], \
    ['X_CO', '%.3g' % X_CO_check, ''], \
    ['Y_CH3OH', '%.3g' % Y_CH3OH_check, ''], \
    ['S_CH3OH_CH4', '%.3g' % S_CH3OH_CH4_check, '']
]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv("reb_3_4_4_results.csv", index=False)
