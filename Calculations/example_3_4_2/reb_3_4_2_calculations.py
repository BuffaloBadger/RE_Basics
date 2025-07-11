import pandas as pd
from scipy import optimize

# Given or known
nDot_N2O5_in = 0.5 # mol/min
nDot_N2_in = 0.5 # mol/min
T = 600.0 # K
P = 5.0 # MPa
y_N2 = 0.3
R = 8314E-6 # L MPa /mol /K

# Solve eqn 7
# equations to solve, as residuals
def residuals(x):
    return [y_N2*(nDot_N2O5_in + nDot_N2_in + 3*x[0]) - nDot_N2_in]
# guess for solution
guess = [0.5]
# try to solve
soln = optimize.root(residuals, guess)
# check for success
if not(soln.success):
    print("A solution was NOT obtained:")
    print(soln.message)
# extract the solution
extent = soln.x

# Calculate the molar flow rate of N2O5, eqn 8
nDot_N2O5 = nDot_N2O5_in - 2*extent

# Calculate the concentration of NO2, eqn 6
C_NO2 = 4*P*(nDot_N2O5_in - nDot_N2O5) \
    /(R*T*(5*nDot_N2O5_in + 2*nDot_N2_in - 3*nDot_N2O5))
C_NO2 = float('%.3g' % C_NO2)

# Display the result
print(' ')
print('Concentration of NO2:',C_NO2)
print(' ')

# Save the result to a .csv file
data = [['C_NO2', C_NO2, 'mol/L']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv("reb_3_4_2_results.csv", index=False)
