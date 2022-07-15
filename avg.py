import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
#from scipy.special import comb as choose
#from scipy.special import logsumexp as fix
import math

def GetnewAvgPD(tree, k,post_order_nodes,node_lookup,choosearr,countterm):


    # create array of (|v|k) size with all values = inf
    arr = np.full((len(post_order_nodes), k + 1), 0.0)
    arr1 = np.full((len(post_order_nodes), k + 1), 0.0)

    arr[:, 0] = 0  # base case 1- for k=0, d(v,k)=0 regardless of vertex
    arr1[:, 0] = 0

    # calculate D(v,k) for all nodes
    for node in post_order_nodes:
        # base case 2- D(v,k) = 0 for all leaf nodes
        if node.is_terminal():
            arr[node_lookup[node]] = 0
            arr1[node_lookup[node]] = 0
            # D(v.k) for all k = 0
        else:
            # recursive case
            x = node.clades[0]  # right tree
            y = node.clades[1]  # left tree
            v_x = tree.distance(node, x)  # edge (v,x)
            v_y = tree.distance(node, y)  # edge (v,y)
            for p in range(min(k,
                               countterm[node_lookup[node]]) + 1):  # cant select more than |v| nodes at that vertex , thus min(k,|v|)
                # for average the extreme cases always exist that is right and left subtree are always chosen
                # print ("p =", p)

                sum_avg = 0
                for r in range(max(0, p - countterm[node_lookup[y]]), min(countterm[node_lookup[x]],
                                                                    p) + 1):  # goes over range for possible r values such that r+l =k
                    l = p - r  # l = k-r
                    # sum_avg += choose_dyn[x.count_terminals()][r] * choose_dyn[y.count_terminals()][l] * (
                    # arr[node_lookup[x]][r] + arr[node_lookup[y]][l] + min(1,r) *v_x +min(1,l)* v_y)
                    sum_avg += choosearr[countterm[node_lookup[x]]][r] * choosearr[countterm[node_lookup[y]]][l] * (
                                arr[node_lookup[x]][r] + arr[node_lookup[y]][l] + min(1, r) * v_x + min(1, l) * v_y)

                # print(sum_avg/(choose(node.count_terminals(), p)))
                arr[node_lookup[node]][p] = sum_avg / choosearr[countterm[node_lookup[node]]][p]
                arr1[node_lookup[node]][p] = sum_avg  # normalize over set siz

    return (arr[node_lookup[tree.root]], arr1[node_lookup[tree.root]], arr1)