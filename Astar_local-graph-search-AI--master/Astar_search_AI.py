# Author: Gurpreet Gosal
# University of Waterloo
# ID 20620838
# CS686

# Local Search for Graph Optimization problem

# Stopping criteria for neighbour r: c(r) - c(s) < 0
# Starting Point random
# Neighbourhood function: select to elements and swap

import numpy as np
import re
import string
import operator
import itertools

fh = open("a1_data5.txt", "r")

pattern1 = r'\d'
pattern2 = r'\d+.\d*'
pattern3 = '\{(\d+.\d)(,\s*\d+.\d)*\}'

data_raw = []
for line in fh:
    data_raw.append(re.findall(pattern2, line))

num_vertices = int(data_raw[0][0])

# Determine how many parent sets each vertex has
num_parent_sets = np.zeros((num_vertices,), dtype = np.int)
store_pos1 = np.zeros((num_vertices,), dtype = np.int)
for i in range(0,num_vertices):
    if i == 0:
        num_parent_sets[i] = int(data_raw[1][0].split(' ')[1])
        store_pos1[i] =2
    else:
        num_parent_sets[i] = int(data_raw[store_pos1[i-1]+num_parent_sets[i-1]][0].split(' ')[1])
        store_pos1[i] = store_pos1[i-1]+num_parent_sets[i-1] + 1


vert_parent_sets = []
vert_parent_costs = []
store_pos2 = np.copy(store_pos1)
store_pos2 = store_pos2 - 1
store_pos3 = np.copy(store_pos1)
store_pos3[0] = store_pos3[len(store_pos1)-1] + num_parent_sets[len(store_pos1)-1] + 1

for i in range(2,len(data_raw)):
    if (0 == sum(i==store_pos2)) & (sum(i == store_pos3-2) != 1):
        if len(data_raw[i]) == 2:
            
            if  '}' in data_raw[i][0]:
                if len(data_raw[i][0]) == 2:
                    vert_parent_sets.append([data_raw[i][0][0]])
                else:
                    vert_parent_sets.append(''.join([data_raw[i][0][0],data_raw[i][0][1]]))
                
                vert_parent_costs.append(data_raw[i][1])
            
            else:
                vert_parent_sets.append( data_raw[i][0].split(',') )
                vert_parent_costs.append(data_raw[i][1])

        elif len(data_raw[i]) == 3:
            if '}' in data_raw[i][1]:
                l_ob1 =  data_raw[i][0].split(',')
                
                if len(data_raw[i][1]) == 2:
                    l_ob2 =  list(data_raw[i][1])[0]
                elif len(data_raw[i][1]) == 3:
                    l_ob2 = ''.join([data_raw[i][1][0],data_raw[i][1][1]])
                
                vert_parent_sets.append([l_ob1[0],l_ob1[1], l_ob2])
                vert_parent_costs.append(data_raw[i][2])
            
            else:
                l_ob1 = data_raw[i][0].split(',')
                l_ob2 = data_raw[i][1].split(',')
                vert_parent_sets.append([l_ob1[0], l_ob1[1], l_ob2[0], l_ob2[1]])
                vert_parent_costs.append(data_raw[i][2])

        elif len(data_raw[i]) == 4:
            l_ob1 = data_raw[i][0].split(',')
            l_ob2 = data_raw[i][1].split(',')
            
            if len(data_raw[i][2]) == 2:
                l_ob3 = list(data_raw[i][2])[0]
            elif len(data_raw[i][2]) == 3:
                l_ob3 = ''.join([data_raw[i][2][0], data_raw[i][2][1]])
    
            vert_parent_sets.append([l_ob1[0], l_ob1[1], l_ob2[0], l_ob2[1], l_ob3])
            vert_parent_costs.append(data_raw[i][3])

    elif (0 == sum(i==store_pos2)) & (sum(i == store_pos3-2) == 1):
        vert_parent_sets.append([])
        vert_parent_costs.append(data_raw[i][0])



# Now we have all the data regarding the graph, lets implement local search

# begin with an initial ordering as a random permutation of  {1,2,3,......n}

iterations = 30
j= 0
initial_order = np.zeros((num_vertices,), dtype = np.int)
initial_order = np.random.permutation(range(1,num_vertices+1))
curr_order = np.copy(initial_order) # current order is equal to initial ordering at the beginning

neighbourhood = 1

cost_record_iter = []

# compute total cost of initial ordering

mc_current = np.zeros((num_vertices,), dtype=np.float)
r = 0
for k in curr_order:
    
    parent_index = np.zeros((num_vertices,), dtype=np.int)
    for p in range(1, num_vertices):
        parent_index[p] = parent_index[p - 1] + num_parent_sets[p - 1]
    
    if k == curr_order[0]:  # for the first vertex in the current ordering the mc is the cost of empty set
        mc_current[k - 1] = vert_parent_costs[parent_index[k - 1] - 1 + num_parent_sets[k - 1]]
    else:
        # go through each parent set of k_th vertex in the ordering and check if it is consistent
        # if consistent raise the flag
        flag_set_k = np.zeros((num_parent_sets[k - 1],), dtype=np.int)
        set_no = 0
        set_k = []
        # for set_k in vert_parent_sets[parent_index[k-1]:parent_index[k-1]+num_parent_sets[k-1]]:
        for ind_set_k in range(0, len(vert_parent_sets[parent_index[k - 1]:parent_index[k - 1]
                                                       + num_parent_sets[k - 1]])):
            if isinstance(vert_parent_sets[ind_set_k + parent_index[k - 1]], str):
                set_k = ''.join(vert_parent_sets[ind_set_k + parent_index[k - 1]])
                set_k = [set_k]
            else:
                set_k = vert_parent_sets[ind_set_k + parent_index[k - 1]]
            
                if not set_k:
                    flag_set_k[set_no] = 1
                elif (len(set_k) == 1):
                    if int(set_k[0]) in curr_order[0:r]:
                        flag_set_k[set_no] = 1
                elif (len(set_k) == 2):
                    if (int(set_k[0]) in curr_order[0:r]) & (int(set_k[1]) in curr_order[0:r]):
                        flag_set_k[set_no] = 1
                elif (len(set_k) == 3):
                    if (int(set_k[0]) in curr_order[0:r]) & (int(set_k[1]) in curr_order[0:r]) \
                            & (int(set_k[2]) in curr_order[0:r]):
                        flag_set_k[set_no] = 1
                elif (len(set_k) == 4):
                    if (int(set_k[0]) in curr_order[0:r]) & (int(set_k[1]) in curr_order[0:r]) \
                            & (int(set_k[2]) in curr_order[0:r]) & (int(set_k[3]) in curr_order[0:r]):
                        flag_set_k[set_no] = 1
                elif (len(set_k) == 5):
                    if (int(set_k[0]) in curr_order[0:r]) & (int(set_k[1]) in curr_order[0:r]) \
                            & (int(set_k[2]) in curr_order[0:r]) & (int(set_k[3]) in curr_order[0:r]) \
                            & (int(set_k[4]) in curr_order[0:r]):
                        flag_set_k[set_no] = 1

            set_no += 1

        store_cmp = []
        for t in range(0, len(flag_set_k)):
            if flag_set_k[t] == 1:
                store_cmp.append(vert_parent_costs[t + parent_index[k - 1]])
        
        mc_current[k - 1] = min(store_cmp)
    
    r += 1

# Now compute the total cost for the active node

total_cost_current = sum(mc_current)

# begin iterations
while (j <= iterations):
    j += 1

    
    previous_order = np.copy(curr_order)
    
    # generate neighbours of current node using neighbourhood function
    if neighbourhood == 1: # generate neighbours by swapping two elements in the ordering
        nbhd_nodes_comb = list(itertools.combinations(curr_order,2))
        
        # for each neighbourhood node first determine its ordering and then compute minimum cost
        total_cost = []
        ii=0
        all_neighbours = []
        r = 0
        for node_comb in nbhd_nodes_comb:
            # first obtain the ordering for current node
            branch_node_curr = np.copy(curr_order)
            comb_index1 = np.where(curr_order == node_comb[0])
            comb_index2 = np.where(curr_order == node_comb[1])
            branch_node_curr[comb_index1] = node_comb[1]
            branch_node_curr[comb_index2] = node_comb[0]
            all_neighbours.append([branch_node_curr])
            # compute cost
            mc = np.zeros((num_vertices,), dtype = np.float)
            
            for k in branch_node_curr:
                
                parent_index = np.zeros((num_vertices,), dtype = np.int)
                for p in range(1,num_vertices):
                    parent_index[p] = parent_index[p-1] + num_parent_sets[p-1]
                
                if k == branch_node_curr[0]: # for the first vertex in the current ordering the mc is the cost of empty set
                    mc[k-1] = vert_parent_costs[ parent_index[k-1] - 1 + num_parent_sets[k-1]]
                else:
                    # go through each parent set of k_th vertex in the ordering and check if it is consistent
                    # if consistent raise the flag
                    flag_set_k = np.zeros((num_parent_sets[k-1],),dtype=np.int)
                    set_no = 0
                    set_k = []
                    #for set_k in vert_parent_sets[parent_index[k-1]:parent_index[k-1]+num_parent_sets[k-1]]:
                    for ind_set_k in range(0,len(vert_parent_sets[parent_index[k - 1]:parent_index[k - 1]
                                                                  + num_parent_sets[k - 1]] ) ):
                        if isinstance(vert_parent_sets[ind_set_k + parent_index[k - 1]],str):
                            set_k = ''.join(vert_parent_sets[ind_set_k + parent_index[k - 1]])
                            set_k = [set_k]
                        else:
                            set_k = vert_parent_sets[ind_set_k + parent_index[k - 1]]
                        
                        if not set_k:
                                flag_set_k[set_no] = 1
                        elif (len(set_k) == 1):
                            if int(set_k[0]) in branch_node_curr[0:r]:
                                    flag_set_k[set_no] = 1
                        elif (len(set_k) == 2):
                            if (int(set_k[0]) in branch_node_curr[0:r]) & (int(set_k[1]) in branch_node_curr[0:r]):
                                flag_set_k[set_no] = 1
                        elif (len(set_k) == 3):
                            if (int(set_k[0]) in branch_node_curr[0:r]) & (int(set_k[1]) in branch_node_curr[0:r])\
                                         & (int(set_k[2]) in branch_node_curr[0:r]):
                                flag_set_k[set_no] = 1
                        elif (len(set_k) == 4):
                            if (int(set_k[0]) in branch_node_curr[0:r]) & (int(set_k[1]) in branch_node_curr[0:r]) \
                                        & (int(set_k[2]) in branch_node_curr[0:r])\
                                        & (int(set_k[3]) in branch_node_curr[0:r]) :
                                flag_set_k[set_no] = 1
                        elif (len(set_k) == 5):
                            if (int(set_k[0]) in branch_node_curr[0:r]) & (int(set_k[1]) in branch_node_curr[0:r]) \
                                        & (int(set_k[2]) in branch_node_curr[0:r])& (int(set_k[3]) in branch_node_curr[0:r]) \
                                        & (int(set_k[4]) in branch_node_curr[0:r]):
                                flag_set_k[set_no] = 1

                        set_no += 1
                    store_cmp = []
                    for t in range(0,len(flag_set_k)):
                        if flag_set_k[t] == 1:
                            store_cmp.append(vert_parent_costs[t + parent_index[k-1]])
                                                        
                    mc[k-1] = min(store_cmp)
                                                    
                    # Now compute the total cost for the active node
            total_cost.append(sum(mc))
                # selection of neighbour check
            if total_cost[r] - total_cost_current < 0:
                curr_order = branch_node_curr
                break
            r += 1
                                                                
        if r< len(nbhd_nodes_comb):
            cost_record_iter.append(total_cost[r])
            total_cost_current = total_cost[r]
                                                                    
    if r>= len(nbhd_nodes_comb):
       break


print(cost_record_iter)
print(curr_order)














