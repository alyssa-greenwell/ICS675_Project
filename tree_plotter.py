import matplotlib.pyplot as plt
import argparse
import sys
import os

# --- Configuration ---
TAXON_A = "marbled_cat"
TAXON_B = "mainland_clouded_leopard"

class Node:
    """A simple Node class to hold tree structure and plot coordinates."""
    def __init__(self, name=""):
        self.name = name
        self.children = []  # List of (Node, distance)
        self.x = 0.0        # Horizontal position (distance from root)
        self.y = 0.0        # Vertical position

def parse_edge_file(filepath):
    """Parses a text file containing Child, Parent, Length columns."""
    edges = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                
                parts = line.split()
                
                # Skip header lines or malformed lines
                if parts[0].lower() == "child" or parts[1].lower() == "parent":
                    continue
                
                if len(parts) < 3:
                    continue
                    
                child = parts[0]
                parent = parts[1]
                try:
                    length = float(parts[2])
                    edges.append((child, parent, length))
                except ValueError:
                    continue
    except FileNotFoundError:
        print(f"Error: Input file '{filepath}' not found.")
        sys.exit(1)
        
    return edges

def visualize_tree(edge_list, output_file):
    """
    Builds tree from edges, roots at midpoint between TAXON_A and TAXON_B, and saves plot.
    """
    # 1. Build Undirected Adjacency Graph
    adj = {}
    nodes = {}

    def get_node(name):
        if name not in nodes:
            nodes[name] = Node(name)
        return nodes[name]

    for child_name, parent_name, dist in edge_list:
        u = get_node(child_name)
        v = get_node(parent_name)
        
        if u not in adj: adj[u] = []
        if v not in adj: adj[v] = []
        
        adj[u].append((v, dist))
        adj[v].append((u, dist))


    if TAXON_A not in nodes or TAXON_B not in nodes:
        print("Error: midpoint taxa not found in dataset.")
        sys.exit(1)

    start = nodes[TAXON_A]
    end = nodes[TAXON_B]

    # Find the unique path between start and end (BFS parent tracing)
    from collections import deque

    parent = {start: None}
    parent_edge_length = {start: 0.0}
    q = deque([start])

    while q:
        curr = q.popleft()
        if curr == end:
            break
        for nxt, dist in adj[curr]:
            if nxt not in parent:
                parent[nxt] = curr
                parent_edge_length[nxt] = dist
                q.append(nxt)

    # Reconstruct the path (list of nodes)
    path_nodes = []
    path_edges = []
    x = end

    while x is not None:
        path_nodes.append(x)
        if parent[x] is not None:
            path_edges.append(parent_edge_length[x])
        x = parent[x]

    # Reverse to go start → end
    path_nodes = path_nodes[::-1]
    path_edges = path_edges[::-1]

    # Total distance
    total_dist = sum(path_edges)
    midpoint = total_dist / 2.0

    # Locate the exact insertion point on the path
    cumulative = 0.0
    root = Node("ROOT")
    root_edge_children = []  # will store two children + branch lengths

    for i in range(len(path_edges)):
        e = path_edges[i]
        n1 = path_nodes[i]
        n2 = path_nodes[i + 1]

        if cumulative + e >= midpoint:
            # midpoint is on this edge
            d1 = midpoint - cumulative
            d2 = e - d1

            # ROOT children on each side
            root_edge_children.append((n1, d1))
            root_edge_children.append((n2, d2))
            break

        cumulative += e

    # Build the directed, rooted tree by orienting edges away from ROOT
    new_root = root
    new_root.children = root_edge_children

    visited = {new_root}
    for child, _ in root_edge_children:
        visited.add(child)

    queue = [c for c, _ in root_edge_children]

    while queue:
        curr = queue.pop(0)
        for nxt, dist in adj[curr]:
            if nxt not in visited:
                curr.children.append((nxt, dist))
                visited.add(nxt)
                queue.append(nxt)

    # 4. Calculate Coordinates (Layout)
    def get_leaves(node):
        if not node.children:
            return [node]
        leaves = []
        for child, _ in node.children:
            leaves.extend(get_leaves(child))
        return leaves

    def assign_coordinates(node, counter, curr_x=0):
        # Set X (Depth)
        node.x = curr_x
        
        # Set Y (Vertical Layout)
        if not node.children:
            node.y = counter[0]
            counter[0] += 1
        else:
            # Ladderize: Sort children by size of subtree
            node.children.sort(key=lambda x: len(get_leaves(x[0])), reverse=True)
            
            child_ys = []
            for child, dist in node.children:
                assign_coordinates(child, counter, curr_x + dist)
                child_ys.append(child.y)
            
            # Internal node Y is average of children Ys
            node.y = sum(child_ys) / len(child_ys)

    counter = [0]
    assign_coordinates(new_root, counter)

    # 5. Plotting
    fig, ax = plt.subplots(figsize=(12, 8))
    
    def draw_lines(node):
        x_start, y_start = node.x, node.y
        
        # Label leaves
        if not node.children and node.name != "ROOT":
            ax.text(x_start + 0.01, y_start, f" {node.name}", va='center', fontsize=10)
            
        for child, _ in node.children:
            x_end, y_end = child.x, child.y
            
            # Draw Square Branches
            ax.plot([x_start, x_start], [y_start, y_end], color='black', lw=1.5)
            ax.plot([x_start, x_end], [y_end, y_end], color='black', lw=1.5)
            draw_lines(child)

    draw_lines(new_root)

    # Styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_ticks([]) 
    ax.set_xlabel("Genetic Distance")
    ax.set_title("Phylogenetic Tree (Midpoint Rooted)", fontsize=14, pad=20)
    
    leaves = get_leaves(new_root)
    max_x = max(n.x for n in leaves)
    ax.set_xlim(0, max_x * 1.25)
    ax.set_ylim(-0.5, len(leaves) - 0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Success: Tree plot saved to '{output_file}'")

def main():
    parser = argparse.ArgumentParser(description="Plot a phylogenetic tree from an edge list file.")
    
    # Input file is a positional argument
    parser.add_argument("input_file", help="Path to the text file containing edges (Child Parent Length)")
    
    # Output file is a REQUIRED named argument
    parser.add_argument("-o", "--output", required=True, help="Filename for the output image (e.g., tree.png)")
    
    args = parser.parse_args()
    
    edges = parse_edge_file(args.input_file)
    
    if not edges:
        print("Error: No valid edges found in input file.")
        sys.exit(1)
        
    visualize_tree(edges, args.output)

if __name__ == "__main__":
    main()
