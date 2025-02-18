function example_25_3_2
% Calculations for Example 25.3.2 of Reaction Engineering Basics
    % constants available to all functions
    % given
    k01 = 10.2; %  gal /mol /min
    k02 = 17.0; % gal /mol /min
    E1 = 15300; % J /mol
    E2 = 23700; % kJ /mol
    CA0 = 10; % mol /gal
    CB0 = 12; % mol /gal
    T0 = 350; % K
    V_ZR = 25; % gal
    Vdot0 = 12.5; % gal /min
    dH1 = -12000;% J /mol
    dH2 = -21300; % J /mol
    CpA = 85; % J /mol /K
    CpB = 125; % J /mol /K
    CpD = 200; % J /mol /K
    CpU = 170; % J /mol /K
    fR2 = 0.05;
    Vdot3 = 0.5; % gal /min
    Vdot4 = Vdot3;
    % known
    R = 8.314; % J /mol /K

    % global variables for making quantities available
    gnAin = nan;
    gnBin = nan;
    gnDin = nan;
    gnUin = nan;
    gTin = nan;
    gV = nan;
    gVdot = nan;

    % CSTR residuals function
    function epsilon = CSTR_Residuals(guess)
        % extract the guesses
        nA = guess(1);
        nB = guess(2);
        nD = guess(3);
        nU = guess(4);
        T = guess(5);

        % calculate the rates
        CA = nA/gVdot;
        CB = nB/gVdot;
        k1 = k01*exp(-E1/R/T);
        k2 = k02*exp(-E2/R/T);
        r1 = k1*CA*CB;
        r2 = k2*CA*CB;

        % evaluate the residuals
        epsA = gnAin - nA - gV*(r1 + r2);
        epsB = gnBin - nB - gV*(r1 + r2);
        epsD = gnDin - nD + gV*r1;
        epsU = gnUin - nU + gV*r2;
        epsT = (gnAin*CpA + gnBin*CpB + gnDin*CpD + gnUin*CpU)...
            * (T-gTin) + gV*(r1*dH1 + r2*dH2);

        % combine the residuals as a vector and return
        epsilon = [epsA; epsB; epsD; epsU; epsT];
    end

    % CSTR zone function
    function unknowns = CSTR_zone(nAin, nBin, nDin, nUin, Tin, V, Vdot)
        % make inputs available to the CSTR residuals function
        gnAin = nAin;
        gnBin = nBin;
        gnDin = nDin;
        gnUin = nUin;
        gTin = Tin;
        gV = V;
        gVdot = Vdot;

        % guess the solution
        guess = [0.9*nAin; 0.9*nBin; 0.05*nAin; 0.05*nAin; Tin + 10];

        % solve the ATEs
        [unknowns, flag, message] = solve_ates(@CSTR_Residuals, guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    end

    % stream merge residuals function
    function epsilon = stream_merge_residuals(guess)
        % extract the guess
        nA1 = guess(1);
        nB1 = guess(2);
        nD1 = guess(3);
        nU1 = guess(4);
        T1 = guess(5);

        % solve the CSTR zone design equations
        VR1 = (1-fR2)*V_ZR;
        Vdot1 = Vdot0 + Vdot4;
        Vdot2 = Vdot1;
        unknowns = CSTR_zone(nA1, nB1, nD1, nU1, T1, VR1, Vdot1);

        % calculate the stream 3 flows and T
        fExch = Vdot3/Vdot2;
        nA3 = fExch*unknowns(1);
        nB3 = fExch*unknowns(2);
        nD3 = fExch*unknowns(3);
        nU3 = fExch*unknowns(4);
        T3 = unknowns(5);

        % solve the well-mixed stagnant zone design equations
        VR2 = fR2*V_ZR;
        unknowns = CSTR_zone(nA3, nB3, nD3, nU3, T3, VR2, Vdot3);

        % evaluate the residuals
        nA0 = CA0*Vdot0;
        nB0 = CB0*Vdot0;
        eps21 = nA0 + unknowns(1) - nA1;
        eps22 = nB0 + unknowns(2) - nB1;
        eps23 = unknowns(3) - nD1;
        eps24 = unknowns(4) - nU1;

        % combine the residuals as a vector and return
        epsilon = [eps21; eps22; eps23; eps24];
    end

    % stream merge function
    function [stream1] = stream_merge()
        % generate a guess
        nA = Vdot0*CA0 + 0.05*Vdot0*CA0;
        nB = Vdot0*CB0 + 0.05*Vdot0*CA0;
        nD = 0.05*Vdot0*CA0;
        nU = 0.05*Vdot0*CA0;
        T = T0 + 5;
        guess = [nA; nB; nD; nU; T];

        % solve the stream merge design equations
        [stream1, flag, message] = solve_ates(@stream_merge_residuals...
            , guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    end

    % zoned reactor function
    function [nA5, nB5, nD5, nU5, T5] = zoned_reactor()
        % solve the stream merge balances
        stream1 = stream_merge();
        nA1 = stream1(1);
        nB1 = stream1(2);
        nD1 = stream1(3);
        nU1 = stream1(4);
        T1 = stream1(5);
        VR1 = (1-fR2)*V_ZR;
        Vdot1 = Vdot0 + Vdot4;
        Vdot2 = Vdot1;

        % solve the CSTR zone design equations
        unknowns = CSTR_zone(nA1, nB1, nD1, nU1, T1, VR1, Vdot1);
        fExch = Vdot3/Vdot2;
        nA5 = (1-fExch)*unknowns(1);
        nB5 = (1-fExch)*unknowns(2);
        nD5 = (1-fExch)*unknowns(3);
        nU5 = (1-fExch)*unknowns(4);
        T5 = unknowns(5);
    end

    % quantities of interest
    function quantities_of_interest()
        [nA5, ~, nD5, nU5, T5] = zoned_reactor();
        fA = 100*(Vdot0*CA0 - nA5)/(Vdot0*CA0);
        SelDU = nD5/nU5;
        item = ["Conversion";"Selectivity";"Temperature"];
        value = [fA; SelDU; T5];
        units = ["%";"mol D per mol U";"K"];
        resultsTable = table(item,value,units);
        writetable(resultsTable,'results.csv');

        disp(resultsTable)
    end

    % perform the analysis
    quantities_of_interest();
end