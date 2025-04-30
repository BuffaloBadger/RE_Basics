""" Calculations for Solving IVODEs using Python """

# import libraries
import numpy as np
from score_utils import solve_ivodes
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

# global constants
V = 1.0 # L
tf = 30.0 # min
CA0 = 0.5 # mol/L
T = 338 # K
R = 8.314E-3 # kJ/mol/K

# global variables
g_k0 = float('nan')
g_E = float('nan')

# derivative function
def derivatives(t,dep):
    # get the dependent variables
    nA = dep[0]

    # calculate the rate coefficient
    k = g_k0*np.exp(-g_E/R/T)

    # calculate the concentration
    CA = nA/V

    # calculate the rate
    r = k*CA

    # calculate the time derivatives of the dependent variables
    dnAdt = -r*V
    dnZdt = r*V

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnZdt])
    return ddt

# BSTR function
def profiles(k0,E):
    # make the rate expression parameters available
    global g_k0, g_E
    g_k0 = k0
    g_E = E

    # set initial values and stopping criterion
    t0 = 0
    nA0 = CA0*V
    n0 = np.array([nA0, 0.0])
    stop_var = 0

    # solve the BSTR reactor design equations
    t, dep, success, message = solve_ivodes(t0, n0, stop_var, tf, derivatives
                                            ,True)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print("WARNING: The ODE solution may not be accurate!")
        print(f"         {message}")
    
    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]

    # return the profiles
    return t, nA, nZ

# quantities of interest function
def quantities_of_interest():
    # set the rate expression parameters
    k0 = 3.6E8 # L/mol/min
    E = 67.5 # kJ/mol

    # get the profiles
    t, nA, nZ = profiles(k0,E)

    # plot the results
    plt.figure()
    plt.plot(t,nA,label='A')
    plt.plot(t,nZ,label='Z')
    plt.xlabel('Time (min)')
    plt.ylabel('Amount (moles)')
    plt.legend()
    plt.savefig('profiles.png')
    plt.show()

if __name__=="__main__":
    quantities_of_interest()