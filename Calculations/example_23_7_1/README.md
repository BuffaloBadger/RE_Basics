This example involves an adiabatic, liquid-phase reaction. The reactor operates at constant pressure and steady-state. The reactor is modeled using the early-mixing segregated flow reactor model. 

The cumulative age distribution function from Example 22.5.1 is used. The space time was 3.3 min in that example.

The liquid-phase reaction in example 20.5.1 is A + B --> Y + Z. kinetics data were generated using r = k*CA*CB/CZ^2 with k0 = 5.29E9 mol/L/min and E = 12100 cal/mol
I changed the reaction to A --> Z and the rate expression to k*CA^2 with k0 = 5.29E9 L/mol/min and E = 12100 cal/mol.

If the feed was at 300K and 0.5 M in A, the isothermal conversion was ca. 76% at a space time of 3.3 min. At 310 K, the isothermal conversion was ca. 82%, so a 10 K temperature rise should work. If the heat capacity is 1.0 cal/ml/K, a heat of reaction of ca. 24.4 kcal/mol should give at 10 K temperature rise, so I used those values.

