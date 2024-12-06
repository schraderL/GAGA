import os
import sys
import argparse

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


species_trait = {}
# Joel aqui tienes que poner los traits que seran labeled como 'test'
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

tree= tree.split(':')
for i in species_trait.keys():
    for j in range(len(tree)):
        if i in tree[j]:
            try:
                tree[j]+=species_trait[i]
            except:
                print('error in: ', i, j)
output = open(output_name, 'w')
output.write(':'.join(tree))
output.close()
