function example_23_7_1
% Calculations for Reaction Engineering Basics Example 23.7.1
    % make given and known constants available to all functions
    VFR_feed = 3.0; % L/min
    T_feed = 300.0; % K
    CA_feed = 0.5; % mol/L
    k0 = 5.29E9; % L/mol/min
    E = 12100; % cal/mol
    dH = -24400; % cal/mol
    Cp = 1000.0; % cal/L/K
    % basis
    delta_tB = 1.0; % min
    % universal constants
    Ren = 1.987; % cal/mol/K
    
    % read F vs. lambda and make it available to all functions
    data_table = readtable(...
        '../example_22_5_1/cum_age_dist_fcn.csv'...
        ,'VariableNamingRule','preserve');
    lambda = table2array(data_table(:,1)); % min
    F = table2array(data_table(:,2));

    % global variables to be made available to all functions
    N_lambda = length(lambda);
    gV0 = nan;
    gT = nan(N_lambda,1);
    gdF = nan(N_lambda,1);
    

    % early mixing segregated flow reactor function
    function [nA_prod, nZ_prod, T_prod] = products()
        % define the hypothetical sample
        V0 = VFR_feed*delta_tB;

        % make V0 available to the derivatives and residuals functions
        gV0 = V0;

        % calculate CA, CZ, and T for each age
        CA = nan(N_lambda,1);
        CA(1) = CA_feed;
        CZ = nan(N_lambda,1);
        CZ(1) = 0;
        T = nan(N_lambda,1);
        T(1) = T_feed;
        for i = 2:N_lambda
            % solve the BSTR design equations
            [~, nA, nZ, Temp] = profiles(lambda(i));
            CA(i) = nA(end)/V0;
            CZ(i) = nZ(end)/V0;
            T(i) = Temp(end);
        end
        
        % calculate dF/d_lambda for each age
        dF = nan(N_lambda,1);
        for i=1:N_lambda - 1
            dF(i) = (F(i+1) - F(i))/(lambda(i+1) - lambda(i));
        end
        dF(N_lambda) = 0.0;

        % calculate the product molar flow rates
        integrandA = nan(N_lambda,1);
        integrandZ = nan(N_lambda,1);
        for i=1:N_lambda
            integrandA(i) = CA(i)*dF(i);
            integrandZ(i) = CZ(i)*dF(i);
        end
        nA_prod = 0.5*V0/delta_tB*trapz(lambda, integrandA);
        nZ_prod = 0.5*V0/delta_tB*trapz(lambda, integrandZ);

        % make T, and dF available to the residuals function
        gT = T;
        gdF = dF;

        % solve the implicit equation for T_prod
        init_guess = T_feed + 10.0;
        [T_prod, flag, message] = solve_ates(@residuals, init_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end

    end

    % residuals function
    function epsilon = residuals(unknown)
        integrand = (unknown - gT).*gdF;
        epsilon = 0.5*gV0*Cp*trapz(lambda,integrand);
    end

    % BSTR function
    function [t, nA, nZ, T] = profiles(tr)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [VFR_feed*CA_feed*delta_tB; 0; T_feed];

        % define the stopping criterion
        f_var = 0;
        f_val = tr;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract the dependent variable profiles
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
    end

    % derivatives function
    function ddt = derivatives(~, dep)
        % extract the necessary dependent variables
        nA = dep(1);
        T = dep(3);

        % calculate the rate
        CA = nA/gV0;
        k = k0*exp(-E/Ren/T);
        r = k*CA^2;

        % evaluate the derivatives
        dnAdt = -gV0*r;
        dnZdt = gV0*r;
        dTdt = -r*dH/Cp;

        % combine the derivatives in a vector and return
        ddt = [dnAdt; dnZdt; dTdt];
    end

    % quantities of interest function
    function quantities_of_interest()
        [nA_prod, ~, T_prod] = products();
        fA = 100*(VFR_feed*CA_feed - nA_prod)/(VFR_feed*CA_feed);
        
        item = ["Conversion"; "Temperature"];
        value = [fA; T_prod];
        units = ["%"; " K"];

        results_table = table(item, value, units);
        disp(' ')
        disp(results_table)
        writetable(results_table, 'results.csv');
    end

    % calculate the quantities of interest
    quantities_of_interest();
end