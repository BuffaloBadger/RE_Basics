function reb_4_5_4

    % an Arrhenius parameters function named Arrhenius_parameters was
    % written and saved in a folder in the Matlab path

    % function that performs the calculations
	function perform_the_calculations()
    
        % define the ideal gas constant
        R = 8.314E-3; % kJ/mol/K
    
        % read the data file
        data_file = "Calculations/example_4_5_4/reb_4_5_4_data.csv";
        data_table = readtable(data_file,'VariableNamingRule','preserve');
        T = table2array(data_table(:,2)); % °C
        k = table2array(data_table(:,3)); % L/mol/min
    
        % create a temperature data set in absolute units
        T = T + 273.15;
    
        % calculate the Arrhenius parameters and statistics
        [k0, k0_ci, E, E_ci, r_squared] = Arrhenius_parameters(k,T,R);

        % generate, show, and save a model plot 
        k_pred = k0*exp(-E/R./T);
        figure
        hold on
        plot(log(k),1./T,'or','MarkerFaceColor','r','MarkerSize',10)
        plot(log(k_pred),1./T,'k',LineWidth=2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('1/T (K^-^1)','FontSize', 14)
        ylabel('k (L mol^-^1 min^-^1)','FontSize', 14)
        saveas(gcf,"Matlab_Arrhenius_plot.png")
    
        % show the results
        disp(' ')
        disp(['k0 = ',num2str(k0,3),' [',num2str(k0_ci(1),3),...
            ', ',num2str(k0_ci(2),3),'] (95%CI) L/mol/min'])
        disp(['E = ',num2str(E,3),' [',num2str(E_ci(1),3),...
            ', ',num2str(E_ci(2),3),'] (95%CI) kJ/mol'])
        disp(['R_squared = ',num2str(r_squared,3)])
    
        % tabulate and save the results
        item = ["k0";"k0_lower_limit";"k0_upper_limit"
            "E";"E_lower_limit";"E_upper_limit";"R_squared"];
        value = round([k0;k0_ci(1);k0_ci(2);E;E_ci(1);
            E_ci(2);r_squared],3,'significant');
        units = ["L/mol/min";"L/mol/min";"L/mol/min";"kJ/mol";"kJ/mol"
            "kJ/mol";""];
        results_table = table(item,value,units);
        writetable(results_table,'Matlab_results.csv');
    end

    % perform the calculations
    perform_the_calculations()
end