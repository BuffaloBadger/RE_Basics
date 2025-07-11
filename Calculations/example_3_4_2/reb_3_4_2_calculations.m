function reb_3_4_2_calculations
    % Given or known
    nDot_N2O5_in = 0.5; % mol/min
    nDot_N2_in = 0.5; % mol/min
    T = 600.0; % K
    P = 5.0; % MPa
    y_N2 = 0.3;
    R = 8314E-6; % L MPa /mol /K

    % Solve eqn7 for the extent
    % equations to solve, as residuals
    function residual = calcResid(x_test)
        residual = [
            y_N2 * (nDot_N2O5_in + nDot_N2_in + 3*x_test) - nDot_N2_in
        ];
    end
    % guess for the solution
    x_guess = [0.5];
    % try to solve
    [x, flag, message] = solveATEs(@calcResid, x_guess);
    % check for success
    if flag <= 0
        disp('The equations were NOT solved.')
        disp(['    Matlab provided the following message: ',message])
    end

    % Calculate the molar flow rate of N2O5, eqn 8
    nDot_N2O5 = nDot_N2O5_in - 2*x;

    % Calculate the concentration of N2O5, eqn 6
    C_NO2 = 4*P*(nDot_N2O5_in - nDot_N2O5)/R/T/...
        (5*nDot_N2O5_in + 2*nDot_N2_in - 3*nDot_N2O5);

    % Display the result
    disp(' ')
    disp(['Outlet concentration of NO2: ', num2str(C_NO2,3),' mol/L'])

    % Save the result to a .csv file
    results_file ="reb_3_4_2_Matlab_results.csv";
    item = "C_NO2";
    value = round(C_NO2,3,'significant');
    units = "mol/L";
    results_table = table(item,value,units);
    writetable(results_table,results_file);
end