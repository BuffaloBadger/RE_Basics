function reb_4_5_2_calculations
    % Given
    k1 = 3.37; % lbmol/h/ft^3/atm^-0.55
    K1 = 12.0;
    y_H2O_0 = 0.75;
    y_CO_0 = 0.24;
    y_CO2_0 = 0.01;
    P = 10; % atm
    T = 675; % K
    R = 0.08206; % L atm/mol/K

    % Basis
    V = 1.0; % L

    % Calculate the initial moles of each reagent, eqns 4 - 7
    n_CO_0 = y_CO_0*P*V/R/T;
    n_H2O_0 = y_H2O_0*P*V/R/T;
    n_CO2_0 = y_CO2_0*P*V/R/T;
    n_H2_0 = 0;

    % Create vectors for plot data
    f_CO = linspace(0.0,1.0,100).';
    r_m0 = nan(100,1);
    r_m1 = nan(100,1);

    % calculate the predicted rates at each conversion
    for i = 1:100
        % calculate the moles of CO, eqn 8
        n_CO = n_CO_0*(1-f_CO(i));
        % calculate the apparent extent, eqn 9
        extent = n_CO_0 - n_CO;
        % calculate the moles of H2O, CO2 and H2, eqns 10 - 12
        n_H2O = n_H2O_0 - extent;
        n_CO2 = n_CO2_0 + extent;
        n_H2 = n_H2_0 + extent;
        % calculate the partial pressures, eqns 13 - 16
        P_CO = n_CO*R*T/V;
        P_H2O = n_H2O*R*T/V;
        P_CO2 = n_CO2*R*T/V;
        P_H2 = n_H2*R*T/V;
        % calculate the rates
        r_m0(i) = k1*P_CO^0.9*P_H2O^0.25*P_CO2^-0.6;
        factor = 1 - P_CO2*P_H2/K1/P_CO/P_H2O;
        r_m1(i) = r_m0(i)*factor;
    end

    % plot the results
    figure; % predicted rates vs conversion
    plot(100*f_CO,r_m0,'r',100*f_CO,r_m1,'b','LineWidth',2)
    yline(0,'LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('CO Conversion (%)','FontSize', 14)
    ylabel('Predicted Rate (lbmol h^-^1 ft^-^3 atm^-^0^.^5^5)',...
        'FontSize', 14)
    legend({'without equilibrium factor','with equilibrium factor'},'Location','northeast','FontSize',14)
    % save the figure
    saveas(gcf,"reb_4_5_2_Matlab_fig_1.png")
    
    figure; % predicted rates vs conversion
    plot(100*f_CO(85:100),r_m0(85:100),'r',100*f_CO(85:100),...
        r_m1(85:100),'b','LineWidth',2)
    yline(0,'LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('CO Conversion (%)','FontSize', 14)
    ylabel('Predicted Rate (lbmol h^-^1 ft^-^3 atm^-^0^.^5^5)',...
        'FontSize', 14)
    legend({'without equilibrium factor','with equilibrium factor'},'Location','northeast','FontSize',14)
    % save the figure
    saveas(gcf,"reb_4_5_2_Matlab_fig_2.png")

    % save the plot data to a .csv file
    data_file = 'reb_4_5_2_Matlab_results.csv';
    data_table = table(100*f_CO,r_m0,r_m1);
    data_table.Properties.VariableNames = ["Conversion", "r2", "r3"];
    writetable(data_table, data_file)
end