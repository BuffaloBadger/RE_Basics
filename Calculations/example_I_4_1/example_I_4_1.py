"""Calculations for the SCoRE Video, Solving ATEs using Python"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd

# global constants
CA_0 = 1.0 # mol /L
V = 50 # L
CB_0 = 1.2 # mol /L
rho = 1.0E3 # g /L
Cp = 1.0 # cal /g /K
T_0 = 303 # K
dH = -10700 # cal /mol
k0 = 8.72E5 # L /mol /min
E = 7200 # cal /mol
Vdot_case = np.array([75,100]) # L /min
R = 1.987 # cal /mol /K

# global variables
g_Vdot = float('nan')

# residuals function
def residuals(guess):
    # extract the individual guesses
    nDotA_1 = guess[0]
    nDotB_1 = guess[1]
    nDotY_1 = guess[2]
    nDotZ_1 = guess[3]
    T_1 = guess[4]

    # calculate r
    k = k0*np.exp(-E/R/T_1)
    CA = nDotA_1/g_Vdot
    CB = nDotB_1/g_Vdot
    r = k*CA*CB

    # evaluate the residuals
    epsilon_1 = g_Vdot*CA_0 - nDotA_1 - r*V
    epsilon_2 = g_Vdot*CB_0 - nDotB_1 - r*V
    epsilon_3 = - nDotY_1 + r*V
    epsilon_4 = - nDotZ_1 + r*V
    epsilon_5 = rho*g_Vdot*Cp*(T_1 - T_0) + r*V*dH

    # return the residuals
    epsilon = np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4, epsilon_5])
    return epsilon

# CSTR function
def unknowns(Vdot):
    # guess the solution
    nA1_guess = 0.9*Vdot*CA_0
    nB1_guess = 0.9*Vdot*CB_0
    nY1_guess = 0.0
    nZ1_guess = 0.0
    T1_guess = T_0 + 5.0
    guess = np.array([nA1_guess, nB1_guess, nY1_guess, nZ1_guess, T1_guess])

    # make Vdot globally available
    global g_Vdot
    g_Vdot = Vdot

    # solve the ATEs
    soln = sp.optimize.root(residuals,guess)

    # check that the solution is converged
    if not(soln.success):
        print("")
        print(f"The solver did NOT converge: {soln.message}")

    # return the solution
    nA = soln.x[0]
    nB = soln.x[1]
    nY = soln.x[2]
    nZ = soln.x[3]
    T = soln.x[4]
    return nA, nB, nY, nZ, T

def quantities_of_interest():
    # solve the ATEs for both cases
    nA1, nB1, nY1, nZ1, T1 = unknowns(Vdot_case[0])
    nA2, nB2, nY2, nZ2, T2 = unknowns(Vdot_case[1])

    # tabulate the results
    data = [['nA',f'{nA1:.1f}',f'{nA2:.1f}','mol/min'],
            ['nB',f'{nB1:.1f}',f'{nB2:.1f}','mol/min'],
            ['nY',f'{nY1:.1f}',f'{nY2:.1f}','mol/min'],
            ['nZ',f'{nZ1:.1f}',f'{nZ2:.1f}','mol/min']
            ,['T',f'{T1:.1f}',f'{T2:.1f}','K']]
    results_df = pd.DataFrame(data, columns=['item',
                                             f'Vdot = {Vdot_case[0]:.0f}',
                                             f'Vdot = {Vdot_case[1]:.0f}',
                                             'units'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('Calculations/example_I_4_1/results.csv',index=False)
    return

# perform the analysis
if __name__=="__main__":
    quantities_of_interest()