import random
import sys
import os
import argparse
from ete3 import Tree


def addAncestor(t):

    for node in t.traverse("levelorder"):
        node_id = node.name
        node_branch = node.dist
        if node_id in dict_stat.keys():
            stat = dict_stat[node_id]
            stat = int(stat)
            stat2 = stat/node_branch
            print(node_id + "\t" + str(stat2) + "\t" + str(node_branch))
            node.add_features(gene=stat2)
        else:
            NA="NA"
            node.add_features(gene=NA)
    print (t.write(features=["gene"]))

            
parser = argparse.ArgumentParser()
parser.add_argument("treeFile",help="input a tree file", type=str)
parser.add_argument("statFile",help="input a stat file", type=str)
args=parser.parse_args()
tree_file = args.treeFile
stat_file = args.statFile

dict_stat = {}
with open(stat_file,'r') as fin:
    for line in fin:
        l = line.strip().split()
        dict_stat[l[0]] = l[1]

with open(tree_file,'r') as f_in:
    for line in f_in:
        nhx = line
        t = Tree(nhx,format=1)
        addAncestor(t)




