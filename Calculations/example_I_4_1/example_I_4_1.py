"""Calculations for Example I.4.1 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd

# global constants
Vdot = 75
CA_0 = 1.0
V = 50
CB_0 = 1.2
rho = 1.0E3
Cp = 1.0
T_0 = 303
dH = -10700
k0 = 8.72E5
E = 7200
R = 1.987

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
    CA = nDotA_1/Vdot
    CB = nDotB_1/Vdot
    r = k*CA*CB

    # evaluate the residuals
    epsilon_1 = Vdot*CA_0 - nDotA_1 - r*V
    epsilon_2 = Vdot*CB_0 - nDotB_1 - r*V
    epsilon_3 = - nDotY_1 + r*V
    epsilon_4 = - nDotZ_1 + r*V
    epsilon_5 = rho*Vdot*Cp*(T_1 - T_0) + r*V*dH

    # return the residuals
    epsilon = np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4, epsilon_5])
    return epsilon

# CSTR function
def unknowns():
    # guess the solution
    nA1_guess = 0.9*Vdot*CA_0
    nB1_guess = 0.9*Vdot*CB_0
    nY1_guess = 0.0
    nZ1_guess = 0.0
    T1_guess = T_0 + 5.0
    guess = np.array([nA1_guess, nB1_guess, nY1_guess, nZ1_guess, T1_guess])

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
    # solve the design equations
    nA, nB, nY, nZ, T = unknowns()

    # tabulate the results
    data = [['nA',f'{nA}'],['nB',f'{nB}'],['nY',f'{nY}'],['nZ',f'{nZ}']
            ,['T',f'{T}']]
    results_df = pd.DataFrame(data, columns=['item','value'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('Calculations/example_I_4_1/results_python.csv'
                      ,index=False)
    return

# perform the analysis
if __name__=="__main__":
    quantities_of_interest()