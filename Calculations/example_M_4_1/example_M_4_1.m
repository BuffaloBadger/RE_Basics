function example_M_4_1
% Calculations for Example M.4.1 of Reaction Engineering Basics
    % constants available to all functions
    % given
    D = 8.0E-6;
    us = 0.01;
    k = 0.012;
    K = 1.0;
    CA0 = 1.0;
    L = 1.25;

    % boundary conditions residuals function
    function epsilon = BC_Residuals(ya, yb)
        % evaluate the residuals
        epsilon_1 = us*ya(1) - D*ya(2) - us*CA0;
        epsilon_2 = yb(2);

        % combine the residuals as a vector and return
        epsilon = [epsilon_1; epsilon_2];
    end

    % derivatives function
    function ddz = derivatives(~, y)
        % evaluate the derivatives
        dy1dz = y(2);
        dy2dz = (1/D)*((k+k/K)*y(1)+us*y(2)-k*CA0/K);
   
        % combine the derivatives in a vector and return
        ddz = [dy1dz; dy2dz];
    end

    % reactor function
    function [z, y] = profiles()
        % set the range
        za = 0.0;
        zb = L;

        % guess the average values of the dependent variables
        yGuess = [CA0; -CA0/L];

        [z, y] = solve_bvodes(za, zb, yGuess, @derivatives...
            , @BC_Residuals);
    end

    % quantities of interest
    function quantities_of_interest()
        % get the reactor profiles
        [z, y] = profiles();

        % plot the results
        figure;
        plot(z,y(1,:),'b',z,y(2,:),'k','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('z','FontSize', 14)
        ylabel('y','FontSize', 14)
        legend({'y_1','y_2'},'Location','southeast','FontSize',14)
        saveas(gcf,"results.png")
    end

    % perform the analysis
    quantities_of_interest();
end