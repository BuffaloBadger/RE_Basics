function reb_3_4_1_calculations
    % Given or known
    y_O2_0 = 0.21 * 0.1;
    conv_O2 = 0.5;
    sel_Ethene_CO2 = 3.0;

    % Basis
    n_Tot_0 = 1.0; % mol

    % Initial mole fractions
    y_Ethene_0 = 0.0; % eqn 7
    y_CO2_0 = 0.0; % eqn 8

    % Initial Moles
    n_O2_0 = y_O2_0 * n_Tot_0; % eqn 4
    n_Ethene_0 = y_Ethene_0 * n_Tot_0; % eqn 5
    n_CO2_0 = y_CO2_0 * n_Tot_0; % eqn 6

    % Solve eqns 3 and 4
    % equations to solve, as residuals
    function residual = calcResid(x_test)
        residual = [
            (x_test(1) + 7 * x_test(2)) - conv_O2 * n_O2_0
            (n_Ethene_0 + 2 * x_test(1)) - sel_Ethene_CO2 * ...
                (n_CO2_0 + 4 * x_test(2))
        ];
    end
    % guess for the solution
    x_guess = [0.5; 0.5];
    % try to solve
    [x, flag, message] = solveATEs(@calcResid, x_guess);
    % check for success
    if flag <= 0
        disp('The equations were NOT solved.')
        disp(['    Matlab provided the following message: ',message])
    end

    % Calculate the mole fraction of CO2
    y_CO2 = (n_CO2_0 + 4 * x(2)) / (n_Tot_0 + x(1) + x(2));

    % Display the result
    disp(' ')
    disp(['Final CO2 mole fraction: ', num2str(y_CO2,3)])

    % Save the result to a .csv file
    results_file ="reb_3_4_1_Matlab_results.csv";
    item = "y_CO2";
    value = round(y_CO2,3,'significant');
    results_table = table(item,value);
    writetable(results_table,results_file);


end