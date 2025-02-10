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
    Vfe = 1.0; % L
    % universal constants
    Ren = 1.987; % cal/mol/K
    
    % read F vs. lambda and make it available to all functions
    data_table = readtable(...
        '../example_22_5_1/cum_age_dist_fcn.csv'...
        ,'VariableNamingRule','preserve');
    lambda = table2array(data_table(:,1)); % min
    F = table2array(data_table(:,2));

    % global variables
    gT = nan;
    gdF = nan;

    % derivatives function
    function ddt = derivatives(~, dep)
        % extract the necessary dependent variables
        nA = dep(1);
        T = dep(3);

        % calculate the rate
        CA = nA/Vfe;
        k = k0*exp(-E/Ren/T);
        r = k*CA^2;

        % evaluate the derivatives
        dnAdt = -Vfe*r;
        dnZdt = Vfe*r;
        dTdt = -r*dH/Cp;

        % combine the derivatives in a vector and return
        ddt = [dnAdt; dnZdt; dTdt];
    end
    
    % fluid element function
    function [CA, CZ, T] = profiles()
        % get the number of integrand evaluation points
        N_lambda = length(F);

        % set the initial values
        ind_0 = 0.0;
        dep_0 = [Vfe*CA_feed; 0.0; T_feed];

        % define the stopping variable
        f_var = 0;
        
        % set the ODE type
        odes_are_stiff = false;

        % allocate storage for the return values
        CA = nan(N_lambda,1);
        CZ = nan(N_lambda,1);
        T = nan(N_lambda,1);

        % set the initial values
        CA(1) = CA_feed;
        CZ(1) = 0.0;
        T(1) = T_feed;

        % loop through the integrand evaluation points
        for i = 2:N_lambda
            % set the stopping criterion
            f_val = lambda(i);
            
            % solve the design equations
            [~, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
                , @derivatives, odes_are_stiff);
    
            % check that the solution was found
            if flag <= 0
                disp(' ')
                disp('WARNING: The ODE solution may not be accurate!')
            end

            % Add to the return value vectors
            CA(i) = dep(end,1)/Vfe;
            CZ(i) = dep(end,2)/Vfe;
            T(i) = dep(end,3);
        end
    end

    % residual function
    function epsilon = residual(unknown)
        integrand = (unknown - gT).*gdF;
        epsilon = 0.5*Vfe*Cp*trapz(lambda,integrand);
    end

    % segregated flow reactor function
    function [nA_prod, nZ_prod, T_prod] = products()
        % get the number of integrand evaluation points
        N_lambda = length(F);

        % calculate CA, CZ, and T for each age
        [CA, CZ, T] = profiles();
        
        % calculate dF/d_lambda for each age
        dF = nan(N_lambda,1);
        for i=1:N_lambda - 1
            dF(i) = (F(i+1) - F(i))/(lambda(i+1) - lambda(i));
        end
        dF(N_lambda) = 0.0;

        % calculate the product molar flow rates
        integrandA = CA.*dF;
        integrandZ = CZ.*dF;
        nA_prod = 0.5*VFR_feed*trapz(lambda, integrandA);
        nZ_prod = 0.5*VFR_feed*trapz(lambda, integrandZ);

        % make T, and dF available to the residuals function
        gT = T;
        gdF = dF;

        % solve the implicit equation for T_prod
        init_guess = T_feed + 10.0;
        [T_prod, flag, message] = solve_ates(@residual, init_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    end

    % quantities of interest function
    function quantities_of_interest()
        [nA_prod, ~, T_prod] = products();
        fA = 100*(VFR_feed*CA_feed - nA_prod)/(VFR_feed*CA_feed);
        
        item = ["Conversion"; "Temperature"];
        value = [fA; T_prod];
        units = ["%"; "K"];

        results_table = table(item, value, units);
        disp(' ')
        disp(results_table)
        writetable(results_table, 'results.csv');
    end

    % calculate the quantities of interest
    quantities_of_interest();
end