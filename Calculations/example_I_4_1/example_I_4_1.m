function example_I_4_1
% Calculations for Example I.4.1 of Reaction Engineering Basics
    % global constants
    Vdot = 75;
    CA_0 = 1.0;
    V = 50.0;
    CB_0 = 1.2;
    rho = 1.0E3;
    Cp = 1.0;
    T_0 = 303;
    dH = -10700;
    k0 = 8.72E5;
    E = 7200;
    R = 1.987;

    % residuals function
    function epsilon = residuals(guess)
        % extract the individual guesses
        nDotA_1 = guess(1);
        nDotB_1 = guess(2);
        nDotY_1 = guess(3);
        nDotZ_1 = guess(4);
        T_1 = guess(5);

        % calculate r
        k = k0*exp(-E/R/T_1);
        CA = nDotA_1/Vdot;
        CB = nDotB_1/Vdot;
        r = k*CA*CB;

        % evaluate the residuals
        epsilon_1 = Vdot*CA_0 - nDotA_1 - r*V;
        epsilon_2 = Vdot*CB_0 - nDotB_1 - r*V;
        epsilon_3 = - nDotY_1 + r*V;
        epsilon_4 = - nDotZ_1 + r*V;
        epsilon_5 = rho*Vdot*Cp*(T_1 - T_0) + r*V*dH;

        % return the residuals
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4; epsilon_5];
    end

    % CSTR function
    function [nA1, nB1, nY1, nZ1, T1] = unknowns()
        % guess the solution
        nA1_guess = 0.9*Vdot*CA_0;
        nB1_guess = 0.9*Vdot*CB_0;
        nY1_guess = 0.0;
        nZ1_guess = 0.0;
        T1_guess = T_0 + 5.0;
        guess = [nA1_guess; nB1_guess; nY1_guess; nZ1_guess; T1_guess];
        
        % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end

        % extract and return the outlet molar flow rates and temperature
        nA1 = soln(1);
        nB1 = soln(2);
        nY1 = soln(3);
        nZ1 = soln(4);
        T1 = soln(5);
    end

    % quantities of interest function
    function quantities_of_interest()
        % solve the design equations
        [nA, nB, nY, nZ, T] = unknowns();

        % tabulate the results
        item = ["nA";"nB";"nY";"nZ";"T"];
        value = [nA; nB; nY; nZ; T];
        results_table = table(item,value);

        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,'results_matlab.csv');
    end

    % perform the analysis
    quantities_of_interest();
end