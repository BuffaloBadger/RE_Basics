function reb_3_4_4_calculations
    % Given or known
    n_CO_0 = 34.0; % mol
    n_H2_0 = 66.0; % mol
    n_CH3OH_0 = 0;
    n_CH4_0 = 0;
    n_H2O_0 = 0;
    n_CO2_0 = 0;
    f_CO = 0.4;
    Y_CH3OH = 0.107;
    S_CH3OH_CH4 = 0.38;

    % Solve eqns 7 through 9 for the apparent extents of reaction
    % equations to solve, as residuals
    function residual = calcResid(x_test)
        residual = [
            f_CO*n_CO_0 - x_test(1) - x_test(2) - x_test(3)
            Y_CH3OH*n_CO_0 - n_CH3OH_0 - x_test(1)
            S_CH3OH_CH4 *(n_CH4_0 + x_test(2)) - n_CH3OH_0 - x_test(1)
        ];
    end
    % guess for the solution
    x_guess = [10, 10, 10];
    % try to solve
    options = optimoptions('fsolve','Display','off');
    [x, ~, flag, details] = fsolve(@calcResid,x_guess,options);
    % check for success
    if flag <= 0
        disp('The equations were NOT solved.')
        disp(['    Matlab provided the following message: ',details.message])
    end

    % Calculate the molar amounts, eqns 10 through 15
    n_CO = n_CO_0 - x(1) - x(2) - x(3);
    n_H2 = n_H2_0 - 2*x(1) - 3*x(2) + x(3);
    n_CH3OH = n_CH3OH_0 + x(1);
    n_CH4 = n_CH4_0 + x(2);
    n_H2O = n_H2O_0 + x(2) - x(3);
    n_CO2 = n_CO2_0 + x(3);

    % Check the results by calculating the conversion, yield, and
    % selectivity
    X_CO_check = (n_CO_0 - n_CO)/n_CO_0;
    Y_CH3OH_check = n_CH3OH/n_CO_0;
    S_CH3OH_CH4_check = n_CH3OH/n_CH4;

    % Display the result
    disp(' ')
    disp('Final Moles')
    disp([' CO: ', num2str(n_CO,3)])
    disp([' H2: ', num2str(n_H2,3)])
    disp([' CH3OH: ', num2str(n_CH3OH,3)])
    disp([' CH4: ', num2str(n_CH4,3)])
    disp([' H2O: ', num2str(n_H2O,3)])
    disp([' CO2: ', num2str(n_CO2,3)])
    disp(' ')
    disp('Check')
    disp([' X_CO: ', num2str(X_CO_check,3)])
    disp([' Y_CH3OH: ', num2str(Y_CH3OH_check,3)])
    disp([' S_CH3OH_CH4_check: ', num2str(S_CH3OH_CH4_check,3)])

    % Save the result to a .csv file
    results_file ="reb_3_4_4_Matlab_results.csv";
    item = ["n_CO";"n_H2";"n_CH3OH";"n_CH4";"n_H2O";"n_CO2";"X_CO";...
        "Y_CH3OH";"S_CH3OH_CH4"];
    value = [round(n_CO,3,'significant'); round(n_H2,3,'significant'); ...
        round(n_CH3OH,3,'significant'); round(n_CH4,3,'significant'); ...
        round(n_H2O,3,'significant'); round(n_CO2,3,'significant'); ...
        round(X_CO_check,3,'significant'); ...
        round(Y_CH3OH_check,3,'significant'); ...
        round(S_CH3OH_CH4_check,3,'significant'); ...
        ];
    units = ["mol";"mol";"mol";"mol";"mol";"mol";"";"";""];
    results_table = table(item,value,units);
    writetable(results_table,results_file);
end