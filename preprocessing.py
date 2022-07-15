import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
from scipy.special import comb as choose
from scipy.special import logsumexp as fix
import math

def preprocessTree(tree,k):
    post_order_nodes = list(tree.find_clades(order='postorder'))
    node_lookup = dict((j, i) for (i, j) in enumerate(post_order_nodes))  # indexes all nodes
    v = len(post_order_nodes)
    countterm= [None]*(v+1)

    for node in post_order_nodes:
        countterm[node_lookup[node]] = node.count_terminals()

    choosearr = np.full((v+1, k + 1), -1.0)
    for i in range(v+1):
        for j in range(k+1):
            choosearr[i][j]= choose(i,j)

    return(post_order_nodes,node_lookup,choosearr,countterm)