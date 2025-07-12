function reb_3_4_1_calculations
    % Global constants
    P = 1.0; % atm
    T = 150 + 273.15; % K
    y_C2H6_0 = 0.9;
    y_air_0 = 0.1;
    conv_O2 = 0.5;
    sel_Ethene_CO2 = 3.0;

    % Basis
    n_Tot_0 = 1.0; % mol

    % calculation of apparent extents of reaction
    n_O2_0 = 0.21*y_air_0*n_Tot_0;

    function epsilon = residuals(guess)
        epsilon = [conv_O2*n_O2_0 - guess(1) - 7*guess(2);
            sel_Ethene_CO2*4*guess(2) - 2*guess(1)];
    end
    guess = [0.5, 0.5];
    options = optimoptions('fsolve','Display','off');
    [soln, ~, flag, details] = fsolve(@residuals,guess,options);

    % check that a solution was found
    if flag <= 0
        disp(' ')
        disp('WARNING: The ATE solver did not converge')
        disp(['         ' , details.message])
    end

    % extract the solution
    extent1 = soln(1);
    extent2 = soln(2);

    % Calculate the mole fraction of CO2
    y_CO2 = (4 * extent2) / (n_Tot_0 + extent1 + extent2);

    % Display the result
    disp(' ')
    disp(['Final CO2 mole fraction: ', num2str(y_CO2,3)])

    % Save the result to a .csv file
    item = ["extent 1"; "extent 2"; "y_CO2_final"];
    value = [extent1; extent2; y_CO2];
    units = ["mol"; "mol"; ""];
    results_table = table(item,value,units);
    writetable(results_table,"reb_3_4_1_Matlab_results.csv");
end