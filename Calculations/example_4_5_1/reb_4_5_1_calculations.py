import pandas as pd

# Given
rho_bed = 155.0 * 0.0283 # lbm/m^3
k0_NH3 = 1.54e15 # kmol NH3 / m^3 / h

# Calculate k0 for generalized rate normalized per pound
k0 = k0_NH3/2/rho_bed

# Display the result
print(' ')
print('k0 (kmol/lbm/h):','%.3g' % k0)
print(' ')

# Save the result to a .csv file
data = [['k0', '%.3g' % k0]]
result = pd.DataFrame(data, columns=['item','value'])
result.to_csv("reb_4_5_1_results.csv", index=False)
