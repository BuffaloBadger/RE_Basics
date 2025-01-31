function example_22_5_1()
% Calculations for Reaction Engineering Basics Example 22.5.1
    % given and known constants
    V = 10.0; % L
    VFR = 3.0; % L /min
    m_impulse = 5.0; % mg

    % read the experimental data
    data_table = readtable('data.csv','VariableNamingRule','preserve');
    t = table2array(data_table(:,1)); % min
    Cout = table2array(data_table(:,2)); % mg /L

    % integrate using the trapezoid rule
    F = nan(length(t),1);
    F(1) = 0.0;
    for i=2:length(t)
        F(i) = VFR/m_impulse*trapz(t(1:i),Cout(1:i));
    end

    % save the cumulative age distribution function
    results_table = table(t,F);
    results_table.Properties.VariableNames = ["lambda (min)", "F"];
    writetable(results_table,'cum_age_dist_fcn.csv');

    % calculate F vs. lambda for an equivalent ideal CSTR
    t_bar = V/VFR;
    F_cstr = nan(length(t),1);
    F_cstr(1) = 0.0;
    for i=2:length(t)
        F_cstr(i) = 1 - exp(-t(i)/t_bar);
    end

    % plot F vs. lambda (t)
    figure;
    plot(t,F,'ok',t,F_cstr,'b','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Age (min)','FontSize', 14)
    ylabel('F','FontSize', 14)
    legend({'From Response','Ideal CSTR'},'Location',...
        'southeast','FontSize',14)
    saveas(gcf,"cum_age_dist.png")

    % plot Cout vs. t
    figure;
    plot(t,Cout,'ok','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Elapsed Time (min)','FontSize', 14)
    ylabel('Tracer Concentration (mg L^-^1','FontSize', 14)
    saveas(gcf,"response.png")
end