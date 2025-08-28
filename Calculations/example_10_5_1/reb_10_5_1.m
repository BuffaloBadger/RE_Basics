function reb_10_5_1
% Calculations for Reaction Engineering Basics Example 10.5.1
    % global constants
    CA_0 = 15.0; % mol/L
    T_0 = 20 + 273.15; % K
    CB_in = 15.0; % mol/L
    T_in = T_0;
    V_0 = 5.0; % L
    Vdot_in = 0.25/60; % L/s
    V_B = 5.0; % L
    V_ex = 2.5; % L
    Tex_in = T_0;
    Vdot_ex = 2.5/60; % L/s
    A = 2150; % cm^2
    U = 73.0/60/60; % cal/cm^2/s/K
    Tex_0 = T_0;
    rho = 1.0E3; % g/L
    Cp = 1.0; % cal/g/K
    dH_1 = -13700; % cal/mol
    k0_1 = 8.11E12; % L/mol/s
    E_1 = 17700; % cal/mol
    P = 1.0; % atm
    Re = 1.987; % cal/mol/K
    Rw = 0.08206; % L*atm/mol/K
    M_B = 40.0; % g/mol

    % derivatives function
    function ddt = derivatives(t, dep)
        % extract dependent variables that are needed
        nA = dep(1);
        nB = dep(2);
        T = dep(5);
        T_ex = dep(6);
        
        % calculate the reacting fluid volume and rate coefficient
        V = Vdot_in*t + V_0;
        k = k0_1*exp(-E_1/Re/T);

        % calculate the rate
        CA = nA/V;
        CB = nB/V;
        r = k*CA*CB;

        % calculate the rate of heat exchange
        Qdot = U*A*(T_ex-T);

        % calculate the flow rates
        nDotB_in = Vdot_in*CB_in;
        mDotB_in = nDotB_in*M_B;
        mDot_ex = Vdot_ex*rho;

        % evaluate the derivatives
        dAdt = -V*r;
        dBdt = nDotB_in - V*r;
        dSdt = V*r;
        dWdt = V*r;
        dTdt = (Qdot - mDotB_in*Cp*(T-T_in) - r*V*dH_1 ...
            + P*Vdot_in*Re/Rw)/rho/V/Cp;
        dTexdt = (-Qdot - mDot_ex*Cp*(T_ex - Tex_in))...
            /(rho*V_ex*Cp);

        % return the derivatives
        ddt = [dAdt; dBdt; dSdt; dWdt; dTdt; dTexdt];
    end

    % SBSTR function
    function [t, nA, nB, nS, nW, T, Tex] = SBSTR_variables()
        % set the initial values
        ind_0 = 0.0;
        nA_0 = CA_0*V_0;
        dep_0 = [nA_0; 0.0; 0.0; 0.0; T_0; Tex_0];

        % define the stopping criterion
        t_f = V_B/Vdot_in;
        stop_var = 0;
        stop_val = t_f;
        
        % solve the IVODEs
        odes_are_stiff = true;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, stop_var...
            , stop_val, @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The IVODE solution may not be accurate!')
        end

        % extract corresponding sets of dependent variables
        nA = dep(:,1);
        nB = dep(:,2);
        nS = dep(:,3);
        nW = dep(:,4);
        T = dep(:,5);
        Tex = dep(:,6);
    end

    % deliverables function
    function deliverables()
        % solve the reactor design equations
        [t, nA, ~, ~, ~, T, ~] = SBSTR_variables();
    
        % calculate the other quantities of interest
        V = Vdot_in*t + V_0;
        CA = nA./V;
        T_C = T - 273.15;
        t_min = t/60.0;
    
        % display and save the graphs
        figure; % acid concentration profile
        plot(t_min,CA,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Time (min)','FontSize', 14)
        ylabel('Concentration of A (mol/L)','FontSize', 14)
        saveas(gcf, 'matlab_concentration_profile.png')
    
        figure; % temperature profile
        plot(t_min,T_C,'k','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Time (min)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        saveas(gcf, 'matlab_temperature_profile.png')
    end

    % perform the analysis
    deliverables();
end