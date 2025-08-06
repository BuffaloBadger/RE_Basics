"""Calculations for Reaction Engineering Basics Example 9.6.4"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import pandas as pd
import matplotlib.pyplot as plt
import reb_9_6_4

# read and set the optimum coolant flow rate
results_df = pd.read_csv('results.csv')
mDot_ex_opt = results_df.iat[0,1]
reb_9_6_4.g_mDot_ex=mDot_ex_opt

# solve the reactor design equations using the optimum coolant flow
t_1, nA_1, nZ_1, T_1, Tex_1 = reb_9_6_4.first_stage_BSTR_variables()
t_2, nA_2, _, T_2, _ = reb_9_6_4.second_stage_BSTR_variables(mDot_ex_opt,
                            t_1[-1], nA_1[-1], nZ_1[-1], T_1[-1], Tex_1[-1])

# combine the profiles
t = np.concatenate((t_1, t_2))
nA = np.concatenate((nA_1, nA_2))
T = np.concatenate((T_1, T_2))

# calculate the conversion and instantaneous rate
r_inst = reb_9_6_4.k_0*np.exp(-reb_9_6_4.E/reb_9_6_4.Re/T)\
    *nA/reb_9_6_4.V*1000

# calculate and plot instantaneous rate profile for discussion
plt.figure() # instantaneous rate profile
plt.plot(t,r_inst)
plt.xlabel("Reaction Time (min)")
plt.ylabel("Instantaneous Rate (mmol/L/min)")
plt.savefig('inst_rate_profile.png')
plt.show()