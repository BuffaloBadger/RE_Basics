from score_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

# an Arrhenius parameters function named Arrhenius_parameters was
# written and added to the rebutils library

# function that performs the calculations
def perform_the_calculations():

    # define the ideal gas constant
    R = 8.314E-3 # kJ/mol/K

    # read the data file and convert the temperature to K
    df = pd.read_csv('reb_4_5_4_data.csv')
    df.columns = ['i','T','k']
    T = df['T'].to_numpy() + 273.15
    k = df['k'].to_numpy()

    # calculate the Arrhenius parameters and statistics
    k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k,T,R)

    # generate and save a model plot
    y_pred = k0*np.exp(-E/R/T)
    plt.figure()
    plt.semilogy(1/T,k,color='k',marker='o', ls='none')
    plt.semilogy(1/T,y_pred,color='r')
    plt.xlabel('T$^{-1}$ (K$^{-1}$)')
    plt.ylabel('k (L mol$^{-1}$ min$^{-1}$)')
    # save and show the figure
    plt.savefig('Arrhenius_plot.png')
    plt.show()

    # show the results
    print(' ')
    print(f'k0: {k0:.3g} L/mol/min, 95% CI [{k0_ci[0]:.3g}, {k0_ci[1]:.3g}]')
    print(f'E: {E:.3g} kJ/mol, 95% CI [{E_ci[0]:.3g}, {E_ci[1]:.3g}]')
    print(f'R-squared: {r_squared:.3g}')
    print(' ')

    # tabulate and save the results
    data = [['k0', f'{k0:.3g}', 'L mol^-1^ min^-1^'],
        ['k0_lower_limit', f'{k0_ci[0]:.3g}', 'L mol^-1^ min^-1^'],
        ['k0_upper_limit', f'{k0_ci[1]:.3g}', 'L mol^-1^ min^-1^'],
        ['E', f'{E:.3g}', 'kJ mol^-1^'],
        ['E_lower_limit', f'{E_ci[0]:.3g}', 'kJ mol^-1^'],
        ['E_upper_limit', f'{E_ci[1]:.3g}', 'kJ mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv("reb_4_5_4_results.csv", index=False)

    return

if __name__=="__main__":
    perform_the_calculations()
