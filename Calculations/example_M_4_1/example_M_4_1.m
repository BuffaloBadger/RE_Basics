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

    % derivatives function
    function ddz = derivatives(~, y)
        % evaluate the derivatives
        dy1dz = y(2);
        dy2dz = (1/D)*((k+k/K)*y(1)+us*y(2)-k*CA0/K);
   
        % combine the derivatives in a vector and return
        ddz = [dy1dz; dy2dz];
    end

    % boundary conditions residuals function
    function epsilon = BC_Residuals(ya, yb)
        % extract the boundary values
        y1a = ya(1);
        y2a = ya(2);
        y2b = yb(2);

        % evaluate the residuals
        epsilon_1 = us*y1a - D*y2a - us*CA0;
        epsilon_2 = y2b;

        % combine the residuals as a vector and return
        epsilon = [epsilon_1; epsilon_2];
    end

    % reactor function
    function [z, y1, y2] = profiles()
        % set the initial mesh
        z = linspace(0, L, 20);

        % set the guess
        yGuess = [0.0; 0.0];
        solinit=bvpinit(z,yGuess);

        % solve the BVODEs
        soln = bvp4c(@derivatives, @BC_Residuals, solinit);

        % extract and return the profiles
        z = soln.x;
        y1 = soln.y(1,:);
        y2 = soln.y(2,:);
    end

    % quantities of interest
    function quantities_of_interest()
        % get the reactor profiles
        [z, y1, y2] = profiles();

        % plot the results
        figure;
        plot(z,y1,'k',z,y2,'b','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('z','FontSize', 14)
        ylabel('y','FontSize', 14)
        legend({'y_1','y_2'},'Location','northeast','FontSize',14)
        saveas(gcf,"results_matlab.png")
    end

    % perform the analysis
    quantities_of_interest();
end