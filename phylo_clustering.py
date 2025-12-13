#!/usr/bin/env python
import copy

class UPGMA():
    """
    Uses UPGMA clustering to build a phylogenetic tree given a distance matrix
    and a list of taxon labels.
    """
    
    def __init__(self, matrix, labels):
        self.dm = copy.deepcopy(matrix)
        self.labels = labels[:]
    
    def run(self, output_filename):
        internal_nodes = []
        node_index = 1
        
        while (len(self.labels) > 1):
            # Find the minimum in the matrix
            min_distance = self.dm[0][1]
            child_one_index = 0
            child_two_index = 1
            for i in range(len(self.labels)):
                for j in range(i+1, len(self.labels)):
                    if(self.dm[j][i] < min_distance):
                        min_distance = self.dm[j][i]
                        child_one_index = j
                        child_two_index = i
                        
            # Create a parent node
            internal_nodes.append([node_index, self.labels[child_two_index], self.labels[child_one_index], (min_distance / 2)])
            
            # Update Tree
            # To prevent indices being messed up, the child with a greater index must be deleted first.
            child_to_delete_first = max(child_one_index, child_two_index)
            child_to_delete_second = min(child_one_index, child_two_index)
            del self.dm[child_to_delete_first]
            del self.dm[child_to_delete_second]
            del self.labels[child_to_delete_first]
            del self.labels[child_to_delete_second]
            new_row = []
            for i in range(len(self.labels)):
                dist_one = self.dm[i].pop(child_to_delete_first)
                dist_two = self.dm[i].pop(child_to_delete_second)
                new_dist = (dist_one + dist_two)/2
                self.dm[i].append(new_dist)
                new_row.append(new_dist)
            new_row.append(0)
            self.dm.append(new_row)
            self.labels.append(node_index)
            node_index += 1
        
        # Write to output file
        with open(output_filename, "w") as f:
            for node in internal_nodes:
                output_string = str(node[1]) + "\t" + str(node[0]) + "\t" + str(node[3]) + "\n"
                f.write(output_string)
                output_string = str(node[2]) + "\t" + str(node[0]) + "\t" + str(node[3]) + "\n"
                f.write(output_string)
        
        return

class NeighborJoining():
    """
    Uses Neighbor Joining clustering to build a phylogenetic tree given a 
    distance matrix and a list of taxon labels.
    """
    def __init__(self, matrix, labels):
        self.dm = copy.deepcopy(matrix)
        self.labels = labels[:]
    
    def run(self, output_filename):
        internal_nodes = []
        node_index = 1
        
        while (len(self.dm) > 2):
            #calculate row sums
            r = []
            n = len(self.dm)
            for i in range(n):
                row_sum = 0
                for j in range(n):
                    row_sum += self.dm[i][j]
                r.append(row_sum)
            
            #calculate neighbor-joining matrix
            Q = [[0 for i in range(n)] for j in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    nj_value = (n - 2) * self.dm[i][j] - r[i] - r[j]
                    Q[i][j] = nj_value
                    Q[j][i] = nj_value
            
            # Find minimum in Q
            min_distance = Q[0][1]
            i_min = 0
            j_min = 1
            for i in range(n):
                for j in range(i+1, n):
                    if(Q[i][j] < min_distance):
                        min_distance = Q[i][j]
                        j_min = j
                        i_min = i
            
            # Create a parent node
            delta_i_j = (r[i_min] - r[j_min]) / (n - 2)
            limb_length_i = (self.dm[i_min][j_min] + delta_i_j) / 2
            limb_length_j = (self.dm[i_min][j_min] - delta_i_j) / 2
            internal_nodes.append([node_index, self.labels[i_min], limb_length_i])
            internal_nodes.append([node_index, self.labels[j_min], limb_length_j])
            
            # Update Tree
            # To prevent indices being messed up, the child with a greater index must be deleted first.
            children_dist = self.dm[i][j]
            child_to_delete_first = max(i_min, j_min)
            child_to_delete_second = min(i_min, j_min)
            del self.dm[child_to_delete_first]
            del self.dm[child_to_delete_second]
            del self.labels[child_to_delete_first]
            del self.labels[child_to_delete_second]
            new_row = []
            for i in range(len(self.labels)):
                dist_one = self.dm[i].pop(child_to_delete_first)
                dist_two = self.dm[i].pop(child_to_delete_second)
                new_dist = (dist_one + dist_two - children_dist)/2
                self.dm[i].append(new_dist)
                new_row.append(new_dist)
            new_row.append(0)
            self.dm.append(new_row)
            self.labels.append(node_index)
            node_index += 1
            for i in range(len(self.dm)):
                print(self.dm[i])
        
        # Join last 2 nodes
        limb_length = self.dm[0][1] / 2
        internal_nodes.append([node_index, self.labels[0], limb_length])
        internal_nodes.append([node_index, self.labels[1], limb_length])
        
        # Write output to file
        with open(output_filename, "w") as f:
            for node in internal_nodes:
                output_string = str(node[1]) + "\t" + str(node[0]) + "\t" + str(node[2]) + "\n"
                f.write(output_string)
        
        return