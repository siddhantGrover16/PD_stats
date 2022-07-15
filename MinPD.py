import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
import pdb
import time

def GetMinPD(tree, k,post_order_nodes,node_lookup,countterm):

    v=len(post_order_nodes)
    arr = np.full((v, k + 1),np.inf)
    #taxa_arr = np.full((v,k+1),set(),set)
    arr[:, 0] = 0#base case 1- for k=0, d(v,k)=0 regardless of vertex

    #calculate D(v,k) for all nodes
    for node in post_order_nodes:
        #base case 2- D(v,k) = 0 for all leaf nodes
        if node.is_terminal():
            arr[node_lookup[node]][1:] = 0
            #for i in range(k+1):
               # taxa_arr[node_lookup[node]][1:] = set([node])
            # D(v.k) for all k = 0
        else:
            #recursive case
            x = node.clades[0] #right tree
            y = node.clades[1] #left tree
            v_x = tree.distance(node, x) #edge (v,x)
            v_y = tree.distance(node, y) #edge (v,y)
            for p in range(1,min(k,countterm[node_lookup[node]]) + 1): # cant select more than |v| nodes at that vertex , thus min(k,|v|)
                mindist = np.inf
                mintaxa = set()
                for r in range(max(0, p - countterm[node_lookup[y]]), min(countterm[node_lookup[x]], p) + 1): #goes over range for possible r values such that r+l =k
                    l = p - r #l = k-r
                    dist = arr[node_lookup[x]][r] + arr[node_lookup[y]][l] + v_x * (min(r, 1)) + v_y * (min(l, 1))#if,r is 0, wont chosose (v,x), same for l=0
                    if dist < mindist:
                        mindist = dist
                        #mintaxa = taxa_arr[node_lookup[x]][r].union(taxa_arr[node_lookup[y]][l])


                arr[node_lookup[node]][p] = mindist # fill array with value
                #taxa_arr[node_lookup[node]][p] = mintaxa

    end = time.time()

  #  print(arr[node_lookup[tree.root]][p])
    #print(arr)
    #print(arr[node_lookup[tree.root]])
    return arr[node_lookup[tree.root]]













