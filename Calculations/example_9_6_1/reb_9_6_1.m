function reb_9_6_1
% Calculations for Example 9.6.1 of Reaction Engineering Basics
    % global constants
    dH_1 = -101.2E3; % J/mol
    k_0_1 = 5.11e4 * 3600; % L/mol/h
    E_1 = 74.8e3; % J/mol
    T_0 = 180 + 273.15; % K
    V = 1900.0; % L
    CA_0 = 2.9; % mol/L
    CB_0 = 3.2; % mol/L
    Cp = 1.23 * 4.184; % J/g/K
    rho = 1.02 * 1000.0; % g/L
    t_f = 2.0; % h
    Re = 8.314; % J/mol/K

    % derivatives function
    function derivs = derivatives(~, dep) % ind not needed
        % extract necessary dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        T = dep(5);
        
        % calculate the rate
        CA = nA/V;
        CB = nB/V;
        k = k_0_1*exp(-E_1/Re/T);
        r = k*CA*CB;

        % evaluate the derivatives
        dAdt = -V*r;
        dBdt = -V*r;
        dYdt = V*r;
        dZdt = V*r;
        dTdt = -r*dH_1/rho/Cp;

        % return the derivatives
        derivs = [dAdt; dBdt; dYdt; dZdt; dTdt];
    end

    % reactor model function
    function [t, nA, nB, nY, nZ, T] = BSTR_variables()
        % set the initial values
        ind_0 = 0.0;
        nA_0 = CA_0*V;
        nB_0 = CB_0*V;
        dep_0 = [nA_0; nB_0; 0.0; 0.0; T_0];

        % define the stopping criterion
        stop_var = 0;
        stop_val = t_f;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, stop_var, stop_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The IVODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % deliverables function
    function deliverables()
        % solve the reactor design equations
        [t, nA, nB, nY, nZ, T] = BSTR_variables();
    
        % calculate the other quantities of interest
        CA = nA/V;
        CB = nB/V;
        CY = nY/V;
        CZ = nZ/V;
        k = k_0_1*exp(-E_1/Re./T);
        r = k.*CA.*CB;
        T_C = T - 273.15;
    
        % tabulate the results
        results_table = table(t,nA,nB,nY,nZ,T_C,r);

        % save the results
        writetable(results_table, 'results_matlab.csv');
    
        % display and save the graphs
        figure; % concentration profiles
        hold("on")
        plot(t,CA,t,CB,t,CY,'LineWidth',2)
        plot(t,CZ,':','LineWidth',4)
        set(gca, 'FontSize', 14);
        xlabel('Time (h)','FontSize', 14)
        ylabel('Concentration (mol/L)','FontSize', 14)
        legend({'A','B','Y','Z'}, 'Location', 'east', 'FontSize', 14)
        saveas(gcf, 'concentrations_matlab.png')
    
        figure; % temperature profile
        plot(t,T_C,'k','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Time (h)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        saveas(gcf, 'temperature_matlab.png')
    
        figure; % rate profile
        plot(t,r,'k','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Time (h)','FontSize', 14)
        ylabel('Rate (mol/L/h)','FontSize', 14)
        saveas(gcf, 'rate_matlab.png')
    end

    % execution statement
    deliverables();
end