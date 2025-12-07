#!/usr/bin/env python

class UPGMA():
    """
    Uses UPGMA clustering to build a phylogenetic tree given a distance matrix
    and a list of taxon labels.
    """
    
    def __init__(self, matrix, labels):
        self.dm = matrix
        self.labels = labels
    
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
    