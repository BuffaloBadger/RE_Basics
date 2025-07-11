import pandas as pd
import numpy as np

# Reaction Reaction matrix with column order: CO, H2, CH3OH, CH4, H2O, CO2
rxn_matrix = np.array([ \
    [-1, -2, 1, 0, 0, 0], \
    [-1, -3, 0, 1, 1, 0], \
    [0, 1, 1, -1, -1, 0], \
    [-1, 1, 0, 0, -1, 1], \
    [0, -4, 0, 1, 2, -1], \
    [0, -3, 1, 0, 1, -1] \
])

# % Calculate the number of independent reactions
n_ind = np.linalg.matrix_rank(rxn_matrix)

# Create a working matrix from the first reaction
working_matrix = rxn_matrix[0,]
# set its rank
working_matrix_rank = 1
# add it to the complete independent reaction set
ind_rxn_set = np.array([1])
# set the next reaction to test
next_rxn = 1

# Find remaining independent reactions
while working_matrix_rank < n_ind:
    # create a test matrix
    test_matrix = np.vstack([working_matrix, rxn_matrix[next_rxn,]])
    # calculate the rank
    test_matrix_rank = np.linalg.matrix_rank(test_matrix)
    if test_matrix_rank > working_matrix_rank:
        # the added reaction is independent
        # the test matrix becomes the new working matrix
        working_matrix = test_matrix
        working_matrix_rank = test_matrix_rank
        ind_rxn_set = np.append(ind_rxn_set, next_rxn + 1)
    next_rxn = next_rxn + 1

# Display the result
print(' ')
print('There are', n_ind,'mathematically independent reactions')
print('One complete mathematically independent subset:',ind_rxn_set)
print(' ')

# Save the result to a .csv file
data = [['Number of Independent Reactions', n_ind]]
result = pd.DataFrame(data, columns=['item','value'])
for i in range(0,len(ind_rxn_set)):
    item = 'Independent Reaction {:d}'.format(i + 1)
    value = ind_rxn_set[i]
    result.loc[i+1]=[item,value]
result.to_csv("reb_3_4_3_results.csv", index=False)
