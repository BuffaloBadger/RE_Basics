function example_23_7_2
% Calculations for Reaction Engineering Basics Example 23.7.2
    % constants available to all functions
    % given
    dH = -7200; % cal /mol
    k0 = 1.37E5; % m^3 /mol /min
    E = 11100; % cal /mol
    yA_feed = 0.1;
    yB_feed = 0.65;
    yI_feed = 0.25;
    VFR_feed = 1.0; % m^3 /min
    T_feed = 165 + 273.15; % K
    P = 5; % atm
    Cp_A = 7.6; % cal /mol K
    Cp_B = 8.2; % cal /mol K
    Cp_Z = 13.8; % cal /mol K
    Cp_I = 4.3; % cal /mol K
    % known
    Re = 1.987; % cal /mol K 
    Rw = 8.206E-5; % m^3 atm /mol /K
    % basis
    Vfe = 1.0; % m^3^

    % read the cumulative age distribution function table
    data_table = readtable('cum_age_dist_fcn.csv'...
        ,'VariableNamingRule','preserve');
    lambda = table2array(data_table(:,1))./60.0; % min
    F = table2array(data_table(:,2));

    % global variables to make quantities available
    gT = nan;
    gdF = nan;
    gCA = nan;
    gCB = nan;
    gCZ = nan;
    gCI = nan;
    ggamma = nan;

    % derivatives function
    function ddt = derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        nZ = dep(3);
        nI = dep(4);
        T = dep(5);
        V = dep(6);

        % calculate the rate
        k = k0*exp(-E/Re/T);
        CA = nA/V;
        CB = nB/V;
        r = k*CA*CB;

        % calculate the mass matrix
        M = zeros(6);
        for i=1:4
            M(i,i) = 1.0;
            M(6,i) = Rw*T;
        end
        M(5,5) = nA*Cp_A + nB*Cp_B + nZ*Cp_Z + nI*Cp_I;
        M(5,6) = -P*Re/Rw;
        M(6,5) = Rw*(nA + nB + nZ + nI);
        M(6,6) = -P;

        % calculate the right-hand side vector
        g = zeros(6,1);
        g(1) = -r*V;
        g(2) = -r*V;
        g(3) = r*V;
        g(4) = 0.0;
        g(5) = -r*V*dH;
        g(6) = 0.0;

        % evaluate and return the derivatives
        ddt = M\g;
    end

    % fluid element function
    function [CA, CB, CZ, CI, T, gamma] = profiles()
        % get the number of integrand evaluation points
        N_lambda = length(F);

        % set the initial values
        ind_0 = 0.0;
        n0 = P*Vfe/Rw/T_feed;
        dep_0 = [yA_feed*n0; yB_feed*n0; 0.0; yI_feed*n0; T_feed; Vfe];

        % define the stopping variable
        f_var = 0;
        
        % set the ODE type
        odes_are_stiff = false;

        % allocate storage for the return values
        CA = nan(N_lambda,1);
        CB = nan(N_lambda,1);
        CZ = nan(N_lambda,1);
        CI = nan(N_lambda,1);
        T = nan(N_lambda,1);
        gamma = nan(N_lambda,1);

        % set the initial values
        CA(1) = yA_feed*n0;
        CB(1) = yB_feed*n0;
        CZ(1) = 0.0;
        CI(1) = yI_feed*n0;
        T(1) = T_feed;
        gamma(1) = 1.0;

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
            CA(i) = dep(end,1)/dep(end,6);
            CB(i) = dep(end,2)/dep(end,6);
            CZ(i) = dep(end,3)/dep(end,6);
            CI(i) = dep(end,4)/dep(end,6);
            T(i) = dep(end,5);
            gamma(i) = dep(end,6)/Vfe;
        end
    end

    % residual function
    function epsilon = residual(unknown)
        integrand = (unknown - gT).*gdF.*ggamma ...
            .*(gCA*Cp_A + gCB*Cp_B + gCZ*Cp_Z + gCI*Cp_I);
        epsilon = 0.5*VFR_feed*trapz(lambda,integrand);
    end

    function [nA, nB, nZ, nI, Tprod] = products()
        % get the fluid element profiles
        [CA, CB, CZ, CI, T, gamma] = profiles();
        
        % calculate dF/d_lambda for each age
        N_lambda = length(F);
        dF = nan(N_lambda,1);
        for i=1:N_lambda - 1
            dF(i) = (F(i+1) - F(i))/(lambda(i+1) - lambda(i));
        end
        dF(N_lambda) = 0.0;

        % calculate the product molar flow rates
        integrandA = gamma.*CA.*dF;
        integrandB = gamma.*CB.*dF;
        integrandZ = gamma.*CZ.*dF;
        integrandI = gamma.*CI.*dF;
        nA = 0.5*VFR_feed*trapz(lambda, integrandA);
        nB = 0.5*VFR_feed*trapz(lambda, integrandB);
        nZ = 0.5*VFR_feed*trapz(lambda, integrandZ);
        nI = 0.5*VFR_feed*trapz(lambda, integrandI);

        % make T, and dF available to the residuals function
        gT = T;
        gdF = dF;
        gCA = CA;
        gCB = CB;
        gCZ = CZ;
        gCI = CI;
        ggamma = gamma;

        % solve the implicit equation for T_prod
        init_guess = T_feed + 10.0;
        [Tprod, flag, message] = solve_ates(@residual, init_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    end

    function quantities_of_interest()
        [nA_prod, ~, ~, ~, T_prod] = products();
        n0 = P*VFR_feed/Rw/T_feed;
        fA = 100*(yA_feed*n0 - nA_prod)/(yA_feed*n0);
        
        item = ["Conversion"; "Temperature"];
        value = [fA; T_prod - 273.15];
        units = ["%"; "°C"];

        results_table = table(item, value, units);
        disp(' ')
        disp(results_table)
        writetable(results_table, 'results.csv');
    end

    quantities_of_interest();
end