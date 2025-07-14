function reb_4_5_1_calculations
    % Given
    rho_bed = 155 * 0.0283; % lbm/m^3
    k0_NH3 = 1.54E15; % kmol NH3/ m^3 /h

    % Calculate k0 for reaction rate expression normalized per pound
    k0 = k0_NH3/2/rho_bed;

    % Display the results
    disp(' ')
    disp(['k0: ', num2str(k0,3),' kmol/lbm/h'])

    % Save the result to a .csv file
    results_file ="reb_4_5_1_Matlab_results.csv";
    item = "k0";
    value = round(k0,3,"significant");
    results_table = table(item,value);
    writetable(results_table,results_file);
end