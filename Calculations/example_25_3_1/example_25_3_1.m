function example_25_3_1
    % Calculations for Example 25.3.1 of Reaction Engineering Basics
    % constants available to all functions
    % given
    LNI = 6; % m
    DNI = 7 / 100; % m
    T = 450 + 273.15; % K
    P = 5; % atm
    VFR_0 = 200 * 0.02832; % m^3/h
    yA0 = 0.15;
    yB0 = 0.15;
    yI0 = 0.7;
    yZ0 = 0.0;
    fVR2 = 0.05;
    k_R1 = 2160.0; % mol/(h m^3 atm^0.5)
    k_R2 = 1785; % mol/(h m^3 atm^0.5)
    % known
    Rpv = 0.08206E-3; % m3 atm/(mol K)

    % global variable to make k available to derivatives function
    gk = nan;
    
    % derivatives function
    function ddV = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        nZ = dep(3);
        nI = dep(4);

        % calculate the rate
        PA = nA/(nA + nB + nZ + nI)*P;
        PB = nB/(nA + nB + nZ + nI)*P;
        if PA > 0
            r = gk*PB*sqrt(PA);
        else
            r = 0.0;
        end

        % evaluate the derivatives
        dnAdV = -2*r;
        dnBdV = -r;
        dnZdV = 2*r;
        dnIdV = 0.0;

        % return the derivatives
        ddV = [dnAdV; dnBdV; dnZdV; dnIdV];
    end

    % PFR zone model
    function [V, nA, nB, nZ, nI] = profiles(nAin, nBin, nZin, nIin...
            , Vf, k)
        % make k available to the derivatives function
        gk = k;

        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nAin; nBin; nZin; nIin];

        % define the stopping criterion
        f_var = 0;
        f_val = Vf;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [V, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nA = dep(:,1);
        nB = dep(:,2);
        nZ = dep(:,3);
        nI = dep(:,4);
    end

    % zoned reactor function
    function [nA5, nB5, nZ5, nI5] = zoned_reactor(f_feedR2)
        nA0 = yA0*P*VFR_0/Rpv/T;
        nB0 = yB0*P*VFR_0/Rpv/T;
        nZ0 = yZ0*P*VFR_0/Rpv/T;
        nI0 = yI0*P*VFR_0/Rpv/T;
        nA1 = (1-f_feedR2)*nA0;
        nB1 = (1-f_feedR2)*nB0;
        nZ1 = (1-f_feedR2)*nZ0;
        nI1 = (1-f_feedR2)*nI0;
        VNI=pi*DNI^2*LNI/4.0;
        Vf = (1-fVR2)*VNI;
        k = k_R1;
        [~, nA, nB, nZ, nI] = profiles(nA1, nB1, nZ1, nI1, Vf, k);
        nA3 = nA(end);
        nB3 = nB(end);
        nZ3 = nZ(end);
        nI3 = nI(end);
        nA2 = f_feedR2*nA0;
        nB2 = f_feedR2*nB0;
        nZ2 = f_feedR2*nZ0;
        nI2 = f_feedR2*nI0;
        Vf = fVR2*VNI;
        k = k_R2;
        [~, nA, nB, nZ, nI] = profiles(nA2, nB2, nZ2, nI2, Vf, k);
        nA5 = nA3 + nA(end);
        nB5 = nB3 + nB(end);
        nZ5 = nZ3 + nZ(end);
        nI5 = nI3 + nI(end);
    end

    % function that performs the analysis
    function quantities_of_interest()
        % set a range of f_feedR2 values
        f_feedR2 = linspace(0.0, 0.25);

        % allocate storage for the conversion
        fA = nan(100,1);

        % calculate the inlet molar flow of A
        nA0 = yA0*P*VFR_0/Rpv/T;

        % loop through those values
        for i = 1:100
            % get the outlet molar flow of A
            [nA5, ~, ~, ~] = zoned_reactor(f_feedR2(i));

            % calculate the conversion
            fA(i) = 100*(nA0 - nA5)/nA0;
        end

        % plot the results
        figure; 
        plot(100*f_feedR2,fA,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Percent of Feed to Less Dense Bed','FontSize', 14)
        ylabel('Percent Conversion','FontSize', 14)
        saveas(gcf,"example_25_3_1_results.png")
    end

    % perform the analysis
    quantities_of_interest()
end