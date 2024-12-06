import os
import sys
import json
from Bio import Phylo
from io import StringIO
import argparse
import re
from ete3 import Tree

# Script to label also the internal nodes as test or reference for relax. It uses the Clade info just to use internal nodes that have all test, or reference species from the same clades (convergent cases of trait evolution, this way we do not merge species from independent convergent events)

parser = argparse.ArgumentParser(description='Script that adds label to newick branches given a trait table. It is intended for use before hyphy RELAX.')
parser.add_argument('tree', type = str,
                     help = 'Tree file')
parser.add_argument('traits_table', type = str,
                     help = 'Name of the table cointaining the GAGA_ID column (default 2) and the desired trait in another column.')
parser.add_argument('output', type = str,
                     help = 'Output file name')

args = parser.parse_args()

file = args.tree
traits_table = args.traits_table
output_name = args.output

tree = open(file,'r').read()
traits_table = [i.replace('\n','').split('\t') for i in open(traits_table,'r').readlines()]

##### Read trait table

species_trait = {}
# traits que seran labeled como 'test'
# (en el ejemplo tengo 'N' y 'R' que son 'no' y 'reduced'),
# los demas estaran labeled como 'reference'. Cambia el 7 por el numero
# de columna (empezando por zero) de tu trait de interes.
# Los nombres de test y reference tambien los puedes cambiar a mas especifico
# ej: Polymorphism y Monomorphism.
for i in traits_table[0:]:
    if i[3] in ['test']:
        species_trait[i[2]]='{test}'
    if i[3] in ['reference']:
        species_trait[i[2]]='{reference}'

##### Read tree, and add internal node ids

# Read the Newick tree string
t = Tree(file, format=0)

# Function to add internal node IDs
def add_internal_node_ids(node, internal_node_id):
    if not node.is_leaf():
        node.name = f"Node{internal_node_id[0]}"
        internal_node_id[0] += 1
    for child in node.children:
        add_internal_node_ids(child, internal_node_id)


# Initialize internal node ID
internal_node_id = [1]

# Call the function to add IDs
add_internal_node_ids(t, internal_node_id)

# Write the modified tree to an output file
output_filetree = output_name+"tmpinternalnodes.tree"
with open(output_filetree, "w") as f:
    f.write(t.write(format=1))

#print(f"Modified tree written to {output_filetree}")


##### Get now a table with the Nodes and the species it contains

# Create a dictionary to store leaf nodes for each internal node
internal_nodes_to_leaves = {}

# Create a list to store external leaf nodes
external_leaf_nodes = []

# Function to recursively collect leaf nodes for internal nodes
def collect_leaf_nodes(node):
    if node.is_leaf():
        return [node.name]
    else:
        leaf_names = []
        for child in node.children:
            leaf_names.extend(collect_leaf_nodes(child))
        return leaf_names

# Iterate through internal nodes and collect leaf nodes for each
for node in t.traverse("postorder"):
    if node.is_leaf():
        external_leaf_nodes.append(node.name)
    elif not node.is_root():
        leaf_names = collect_leaf_nodes(node)
        internal_nodes_to_leaves[node.name] = leaf_names

# Write the collected information to the output file
output_file = output_name + "internal_node_leaf_table.tsv" 
with open(output_file, "w") as f:
    f.write("Internal_Node\tLeaf_Nodes\n")
    for internal_node, leaf_nodes in internal_nodes_to_leaves.items():
        leaf_nodes_str = ", ".join(leaf_nodes)
        f.write(f"{internal_node}\t{leaf_nodes_str}\n")

    # Add external leaf nodes to the output file
    #f.write("External_Leaf_Nodes\n")
    for leaf_node in external_leaf_nodes:
        f.write(f"{leaf_node}\t{leaf_node}\n")

#print(f"Internal node leaf table written to {output_file}")


#### Now get the trait info for all leafs and nodes

# Read the first file
with open(output_name+"internal_node_leaf_table.tsv", "r") as f:
    first_file_lines = f.readlines()

# Read the second file and store traits in a dictionary, also the clade
traits_dict = {}
clade_dict = {}
with open(args.traits_table, "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        leaf_id = parts[2]
        trait = parts[3]
        clade = parts[0]
        traits_dict[leaf_id] = trait
        clade_dict[leaf_id] = clade


# Print the traits dictionary for debugging
#print("Traits Dictionary:", traits_dict)

# Prepare the output lines with the third column
output_lines = []
for line in first_file_lines:
    parts = line.strip().split("\t")
    internal_node = parts[0]
    leaf_nodes = parts[1].split(", ")

    # Check traits for all leaf nodes in the row
    #leaf_traits = [traits_dict.get(leaf_node, "Unknown") for leaf_node in leaf_nodes]
    leaf_traits = []
    leaf_clades = []
    for leaf_node in leaf_nodes:
        common_part = leaf_node.split("_")[0]  # Extract the common part of the leaf node
        if common_part in traits_dict:
            leaf_traits.append(traits_dict[common_part])
            leaf_clades.append(clade_dict[common_part])
        else:
            #print(f"Leaf node {leaf_node} not found in traits_dict when using {common_part}")
            leaf_traits.append("Unknown")
            leaf_clades.append("Unknown")

    #print("Leaf Nodes:", leaf_nodes) # Debugging
    #print("Leaf Traits:", leaf_traits) # Debugging
 
    # Determine the value for the third column
    third_column_value = "Unknown"
    if len(set(leaf_clades)) == 1:  # Check if all leaf nodes are from the same clade
        if all(trait == "test" for trait in leaf_traits):
            third_column_value = "test"
        elif all(trait == "reference" for trait in leaf_traits):
            third_column_value = "reference"
        else:
            third_column_value = "mixed"
    else:
            third_column_value = "mixed"        

    #print("Third Column Value:", third_column_value) # Debugging

    output_line = f"{line.strip()}\t{third_column_value}"
    output_lines.append(output_line)

# Write the output to a new file
with open(output_name+"internal_node_leaf_tablec3.txt", "w") as f:
    f.write("\t".join(["Internal_Node", "Leaf_Nodes", "Trait"]) + "\n")
    for line in output_lines:
        f.write(line + "\n")

#print(f"Internal node leaf table with traits written to {output_file}internal_node_leaf_tablec3.txt")


###### Now tag the tree for relax

# Read the tree and new trait table with internal node info
tree = open(output_filetree,'r').read()
output_filec3 = output_name + "internal_node_leaf_tablec3.txt"
traits_table_wnodes = [i.replace('\n','').split('\t') for i in open(output_filec3,'r').readlines()]


species_trait = {}
for i in traits_table_wnodes[0:]:
    #print("Current 'i' value:", i)
    #print("Length of 'i':", len(i))
    if i[2] in ['test']:
        species_trait[i[0]] = '{test}'
    if i[2] in ['reference']:
        species_trait[i[0]] = '{reference}'


tree= tree.split(':')
for i in species_trait.keys():
    for j in range(len(tree)):
        #print(tree[j]) # Debug
        #if i in tree[j]: # like contains, but Node9 could involve Node99
        #if tree[j] == i: # Exact match, as the table with trait and node/leafs here is already containing full names, but splitting with : leaves more data so this does not apply
        if tree[j].endswith(i): # Needs to match at the end
            try:
                tree[j]+=species_trait[i]
            except:
                print('error in: ', i, j)
output = open(output_name, 'w')
output.write(':'.join(tree))
output.close()


# Delete the temporary internal node tree file
#if os.path.exists(output_filetree):
#    os.remove(output_filetree)

# Delete the temporary internal node leaf table file
if os.path.exists(output_name + "internal_node_leaf_table.tsv"):
    os.remove(output_name + "internal_node_leaf_table.tsv")
    #print(f"Temporary internal node leaf table file {output_file}internal_node_leaf_table.tsv deleted.")
