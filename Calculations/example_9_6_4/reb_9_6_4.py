"""Calculations for Reaction Engineering Basics Example 9.6.4"""

# import libraries
import math
import numpy as np
from score_utils import solve_ivodes
import pandas as pd
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

# global constants
k_0 = 2.59E9 # /min
E = 16500. # cal /mol
dH = -22200. # cal /mol
CA_0 = 2. # mol /L
T_0 = 23  + 273.15 # K
Cp = 440. # cal /L /K
V = 4.0 # L
V_shell = 0.5 # L
A_shell = 0.6 # ft^2
U_shell = 1.13E4/60 # cal /ft^2 /min /K
Tex_in = 20 + 273.15 # K
rho_ex = 1.0 # g /cm^3
Cp_ex = 1.0 # cal /g /K
U_coil = 3.8E4/60 # cal /ft^2 /min /K
A_coil = 0.23 # ft^2
T_coil = 120 + 273.15 # K
Tex_0 = 23 + 273.15 # K
T_1 = 50 + 273.15 # K
T_f = 25 + 273.15 # K
t_turn = 25 # min
Re = 1.987 # cal /mol /K

# global variable
g_mDot_ex = float('nan')

# derivatives function for the first stage of operation
def first_stage_derivatives(ind, dep):
	# extract necessary dependent variables for this integration step
    nA = dep[0]
    T = dep[2]
    Tex = dep[3]

	# calculate the rate
    CA = nA/V
    k = k_0*math.exp(-E/Re/T)
    r = k*CA
    
    # calculate the rates of heat exchange
    Qdot_shell = U_shell*A_shell*(Tex - T)
    Qdot_coil = U_coil*A_coil*(T_coil - T)

	# evaluate the derivatives
    dnAdt = -V*r
    dnZdt = V*r
    dTdt = (Qdot_shell + Qdot_coil - V*r*dH)/V/Cp
    dTexdt = -Qdot_shell/rho_ex/V_shell/Cp_ex
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTexdt]

# derivatives function for the second stage of operation
def second_stage_derivatives(ind, dep):
	# extract necessary dependent variables for this integration step
    nA = dep[0]
    T = dep[2]
    Tex_out = dep[3]

	# calculate the rate
    CA = nA/V
    k = k_0*math.exp(-E/Re/T)
    r = k*CA

    # calculate the rates of heat exchange
    Qdot_shell = U_shell*A_shell*(Tex_out - T)

	# evaluate the derivatives
    dnAdt = -V*r
    dnZdt = V*r
    dTdt = (Qdot_shell - V*r*dH)/V/Cp
    dTexdt = -(Qdot_shell + g_mDot_ex*Cp_ex*(Tex_out-Tex_in))/rho_ex/V_shell/Cp_ex
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTexdt]

# BSTR function for the first stage of operation
def first_stage_BSTR_variables():
	# set the initial values
    ind_0 = 0.0
    nA_0 = CA_0*V
    dep_0 = np.array([nA_0, 0.0, T_0, Tex_0])

	# define the stopping criterion
    stop_var = 3
    stop_val = T_1
    
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, stop_var, stop_val
                                        , first_stage_derivatives)
    
    # check that a solution was found
    if not(success):
        print(f"First stage IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    Tex = dep[3,:]

    # return all profiles
    return t, nA, nZ, T, Tex

# BSTR function for the second stage of operation
def second_stage_BSTR_variables(mDot_ex, t_0, nA_0, nZ_0, T_0, Tex_0):
    # make the coolant flow rate globally available
    global g_mDot_ex
    g_mDot_ex = mDot_ex

	# set the initial values
    ind_0 = t_0
    dep_0 = np.array([nA_0, nZ_0, T_0, Tex_0])

	# define the stopping criterion
    stop_var = 3
    stop_val = T_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, stop_var, stop_val
                                        , second_stage_derivatives)

    # check that a solution was found
    if not(success):
        print(f"Second stage IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    Tex_out = dep[3,:]

    # return all profiles
    return t, nA, nZ, T, Tex_out

# deliverables function
def deliverables():
    # solve the reactor design equations for the first stage
    t_1, nA_1, nZ_1, T_1, Tex_1 = first_stage_BSTR_variables()

    # choose a range of coolant flow rates
    coolant_flows = np.linspace(100.0, 250.0, 100)

    # calculate the net rate for each coolant flow rate
    net_rates = np.zeros(100)
    for i in range(0,100):
        # make the coolant flow rate available
        mDot_ex = coolant_flows[i]

        # solve the reactor design equations
        t, _, nZ, _, _ = second_stage_BSTR_variables(mDot_ex, t_1[-1], nA_1[-1]
                                    , nZ_1[-1], T_1[-1], Tex_1[-1])
        
        # calculate the net rate
        net_rates[i] = nZ[-1]/(t[-1] + t_turn)

    # find the coolant flow where the net rate is maximized
    i_max = np.argmax(net_rates)
    mDot_max = coolant_flows[i_max]

    # solve the reactor design equations using the optimum coolant flow
    t_2, nA_2, nZ_2, T_2, _ = second_stage_BSTR_variables(mDot_max, t_1[-1]
                                    , nA_1[-1], nZ_1[-1], T_1[-1], Tex_1[-1])
    
    # combine the profiles
    t = np.concatenate((t_1, t_2))
    nA = np.concatenate((nA_1, nA_2))
    nZ = np.concatenate((nZ_1, nZ_2))
    T = np.concatenate((T_1, T_2))
    
    # calculate the conversion vs time at the optimum coolant flow rate
    nA_0 = CA_0*V
    pct_conversion = 100*(nA_0 - nA)/nA_0
    
    # tabulate the results
    max_net_rate = nZ[-1]/(t[-1] + t_turn)
    data = [["Optimum Coolant Flow", f"{mDot_max}", "g min^-1^"],
            ["Maximum Net Rate", f"{max_net_rate}", "mol min^-1^"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(" ")
    print(f"Optimum Coolant Flow: {mDot_max} g /min")
    print(f"Maximum Net Rate: {max_net_rate} mol /min")

    # save the results
    results_df.to_csv('results.csv', index=False)

    # display and save the graphs
    plt.figure() # net rate vs coolant flow
    plt.plot(coolant_flows, net_rates)
    plt.xlabel("Coolant Flow (g/min)")
    plt.ylabel("Net Rate (mol/min)")
    plt.savefig(
        'net_rate_vs_coolant_flow.png')
    plt.show()

    plt.figure() # conversion profile
    plt.plot(t,pct_conversion)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Conversion (%)")
    plt.savefig('conversion_profile.png')
    plt.show()

    plt.figure() # temperature profile
    plt.plot(t,T - 273.15)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Temperature (°C)")
    plt.savefig('temperature_profile.png')
    plt.show()

    return

if __name__=="__main__":
    deliverables()
