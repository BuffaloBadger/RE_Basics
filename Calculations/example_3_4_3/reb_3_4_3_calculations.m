function reb_3_4_3_calculations
    % Reaction matrix with column order: CO, H2, CH3OH, CH4, H2O, CO2
    rxn_matrix = [
        -1 -2 1 0 0 0
        -1 -3 0 1 1 0
        0 1 1 -1 -1 0
        -1 1 0 0 -1 1
        0 -4 0 1 2 -1
        0 -3 1 0 1 -1];

    % Calculate the number of independent reactions
    n_ind = rank(rxn_matrix);

    % Create a working matrix from the first reaction
    working_matrix = rxn_matrix(1,:);
    % set its rank
    working_matrix_rank = 1;
    % add it to the complete independent reaction set
    ind_rxn_set = 1;
    % set the next reaction to test
    next_rxn = 2;

    % Find remaining independent reactions
    while working_matrix_rank < n_ind
        % create a test matrix
        test_matrix = [working_matrix; rxn_matrix(next_rxn,:)];
        % calculate the rank
        test_matrix_rank = rank(test_matrix);
        if (test_matrix_rank > working_matrix_rank)
            % the added reaction is independent; the test matrix becomes the new
            % working matrix
            working_matrix = test_matrix;
            working_matrix_rank = test_matrix_rank;
            ind_rxn_set = [ind_rxn_set, next_rxn];
        end
        next_rxn = next_rxn + 1;
    end

    % Display the results
    disp(' ')
    disp(['There are ', num2str(n_ind),' mathematically independent reactions'])
    disp(['One complete mathematically independent subset: ',num2str(ind_rxn_set)])

    % Save the result to a .csv file
    results_file ="reb_3_4_3_Matlab_results.csv";
    item = "n_ind";
    value = round(n_ind,0);
    results_table = table(item,value);
    writetable(results_table,results_file);
end