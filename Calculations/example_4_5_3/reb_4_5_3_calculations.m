function reb_4_5_3_calculations
    % Given
    mu_max = 1.0; % /h
    Ks = 0.2; % g/L

    % Create vectors for plot data
    C_S = linspace(0.0,5.0,100);
    mu = nan(100,1);

    % calculate the predicted rates at each conversion
    for i = 1:100
        % calculate the rates
        mu(i) = mu_max*C_S(i)/(Ks + C_S(i));
    end

    % plot the results
    figure; % specific rate vs. substrate concentration
    plot(C_S,mu,'b','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Substrate Mass Concentration (g L^-^1)','FontSize', 14)
    ylabel('Specific Growth Rate (h^-^1)', 'FontSize', 14)
    % save the figure
    saveas(gcf,"reb_4_5_3_Matlab_fig_1.png")
    
end