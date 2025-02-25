"""Calculations for Example 12.7.1 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi = 300)

# constants available to all functions
# given
D = 8.0E-6
us = 0.01
k = 0.012
K = 1.0
CA0 = 1.0
L = 1.25

# derivatives function
def derivatives(z, y):
    # extract the mesh size and dependent variables
    nMesh = z.size
    y1 = y[0,:]
    y2 = y[1,:]

    # allocate space for dydz
    dydz = np.zeros((y.shape))

    # evaluate the derivatives at each mesh point
    for i in range(0,nMesh):
        dydz[0,i] = y2[i]
        dydz[1,i] = (1/D)*((k+k/K)*y1[i]+us*y2[i]-k*CA0/K)

    # return the derivatives
    return dydz

# boundary conditions residuals function
def residuals(ya, yb):
    # extract the boundary values
    y1a = ya[0]
    y2a = ya[1]
    y1b = yb[0]
    y2b = yb[1]

    # evaluate the residuals
    epsilon_1 = us*y1a - D*y2a - us*CA0
    epsilon_2 = y2b

    # return the residuals
    return np.array([epsilon_1, epsilon_2])

# reactor function
def profiles():
	# set the initial mesh
    z = np.linspace(0, L, 20)

    # set the guess
    y = np.zeros((2,20))

    # solve the bvodes
    soln = solve_bvp(derivatives, residuals, z, y)

    # extract and return the profiles
    return soln.x, soln.y[0,:], soln.y[1,:]

# perform the analysis
def quantities_of_interest():
	# solve the design equations
    z, y1, y2 = profiles()

    # report the initial and final mesh sizes
    print(f"Size of initial mesh: 20, size of final mesh: {len(z)}")
    
    # plot the results
    plt.figure(1) 
    plt.plot(z, y1, color = 'k', label = 'y1') 
    plt.plot(z, y2, color = 'b', label = 'y2')
    plt.xlabel("z")
    plt.ylabel("y")
    plt.legend()

    # save and show the figure
    plt.savefig('Calculations/example_M_4_1/results.png')
    plt.show()
    return

if __name__=="__main__":
    quantities_of_interest()
