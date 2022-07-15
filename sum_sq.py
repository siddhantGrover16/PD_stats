import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
from scipy.special import comb as choose
from scipy.special import logsumexp as fix
import math
from math import log
from scipy.special import gammaln

def Getsumsq(tree, k,beta,post_order_nodes,node_lookup,comb,countterm):
    v = len(post_order_nodes)
    k = k

    #create array of (|v|k) size with all values = inf
    arr = np.full((v, k + 1), 0.0)
    arr1 = np.full((v, k + 1), 0.0)

    arr[:, 0] = 0 #base case 1- for k=0, d(v,k)=0 regardless of vertex
    arr1[:, 0] = 0


    #calculate D(v,k) for all nodes
    for node in post_order_nodes:
        #base case 2- D(v,k) = 0 for all leaf nodes
        if node.is_terminal():
            arr[node_lookup[node]] = 0
            arr1[node_lookup[node]] = 0
            # D(v.k) for all k = 0
        else:
            #recursive case
            x = node.clades[0] #right tree
            y = node.clades[1] #left tree
            v_x = tree.distance(node, x) #edge (v,x)
            v_y = tree.distance(node, y) #edge (v,y)
            lamb = v_x+v_y
            lamb2= lamb*lamb
            for p in range(min(k,countterm[node_lookup[node]])+1): # cant select more than |v| nodes at that vertex , thus min(k,|v|)
                #for average the extreme cases always exist that is right and left subtree are always chosen
               # print ("p =", p)
                sums_x = 0
                sums_y = 0
                sums = 0
                for r in range(max(0, p - countterm[node_lookup[y]]), min(countterm[node_lookup[x]],p) + 1):  # goes over range for possible r values such that r+l =k
                    l = p - r  # l = k-r
                    sums += comb[countterm[node_lookup[y]]][l]*(arr[node_lookup[x]][r]) + comb[countterm[node_lookup[x]]][r]*(arr[node_lookup[y]][l]) + comb[countterm[node_lookup[y]]][l] *comb[countterm[node_lookup[x]]][r]*(lamb2) + 2*(beta[node_lookup[x]][r]*beta[node_lookup[y]][l] +lamb*(comb[countterm[node_lookup[y]]][l]*beta[node_lookup[x]][r]+comb[countterm[node_lookup[x]]][r]*beta[node_lookup[y]][l]))
                    #print("1",arr[node_lookup[x]][r])
                    #print("2", arr[node_lookup[y]][l])
                    #print("3", beta[node_lookup[y]][l])
                    #print("4", arr[node_lookup[x]][r])

                sums += arr[node_lookup[x]][p] + comb[countterm[node_lookup[x]]][p]*(v_x)*(v_x) + 2*(v_x)*(beta[node_lookup[x]][p])
                sums += arr[node_lookup[y]][p] + comb[countterm[node_lookup[y]]][p] * (v_y) * (v_y) + 2 * (v_y) * (beta[node_lookup[y]][p])


               # print(sum_avg/(choose(node.count_terminals(), p)))
                arr[node_lookup[node]][p] = sums


    return (arr[node_lookup[tree.root]])