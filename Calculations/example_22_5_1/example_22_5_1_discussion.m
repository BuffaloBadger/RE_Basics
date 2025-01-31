function example_22_5_1_discussion()
% Calculations for the Reaction Engineering Basics Example 22.5.1
% discussion
    % given and known constants
    VFR = 3.0; % L /min
    m_impulse = 5.0; % mg

    % read the experimental data
    data_table = readtable('data.csv','VariableNamingRule','preserve');
    t = table2array(data_table(:,1)); % min
    Cout = table2array(data_table(:,2)); % mg /L

    % only use the first 10 data
    t = t(1:10);
    Cout = Cout(1:10);

    % integrate using the trapezoid rule
    F = nan(length(t),1);
    F(1) = 0.0;
    for i=2:length(t)
        F(i) = VFR/m_impulse*trapz(t(1:i),Cout(1:i));
    end

    % repeat using only the first 10 data
    m10 = VFR*trapz(t(1:10),Cout(1:10));
    F10 = nan(10,1);
    F10(1) = 0.0;
    for i=2:10
        F10(i) = VFR/m10*trapz(t(1:i),Cout(1:i));
    end
    
    % report and save the results
    item = ["True Impulse Mass"; "Eqn 22.4 Impluse Mass"];
    value = [m_impulse; m10];
    units = ["mg"; "mg"];
    results_table = table(item,value,units);
    disp(results_table)
    writetable(results_table,'results.csv');

    % plot F vs. lambda (t)
    figure;
    hold on
    plot(t,F,'ok','LineWidth',2)
    plot(t,F10,'xr','MarkerSize',8,'LineWidth',2)
    hold off
    set(gca, 'FontSize', 14);
    xlabel('Age (min)','FontSize', 14)
    ylabel('F','FontSize', 14)
    legend({'True Impulse Mass','Eqn 22.4 Impluse Mass'},'Location',...
        'southeast','FontSize',14)
    saveas(gcf,"cum_age_dist_comp.png")
end