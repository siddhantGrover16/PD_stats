import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
from scipy.special import comb as choose
from decimal import *

def GetSums_Avg(tree, k,post_order_nodes,node_lookup,choosearr,countterm):
    """
    Given a Bio.Phylo tree and the preprocessed data, computes the average PD statistics on a tree
    :param tree: input tree
    :param k: user selction for choice of k leaves
    [the following is the preprocessed data on the tree]
    :param post_order_nodes:
    :param node_lookup:
    :param choosearr:
    :param countterm:
    :return: a n*k matrix with average PD rooted at node n for all choices ranging from 1 to k
    """
    if (k == 0):
        return Decimal(0)

    beta = np.full((len(post_order_nodes), k + 1), Decimal(0))
    alpha = np.full((len(post_order_nodes), k + 1), Decimal(0))
    avgvals =[]
    getcontext().prec =64

    beta[:, 0] = Decimal(0)  # sums_arr_beta
    alpha[:, 0] = Decimal(0)  #avg_arr_alpha

    # calculate D(v,k) for all nodes
    for node in post_order_nodes:
        # base case 2- D(v,k) = 0 for all leaf nodes
        if node.is_terminal():
            beta[node_lookup[node]] =Decimal(0)
            alpha[node_lookup[node]] =Decimal(0)
            # D(v.k) for all k = 0
        else:
            # recursive case
            x = node.clades[0]  # right tree
            y = node.clades[1]  # left tree
            v_x = Decimal(tree.distance(node, x))  # edge (v,x)
            v_y = Decimal(tree.distance(node, y))  # edge (v,y)
            lamb = Decimal(v_x+v_y)
            for p in range(min(k, countterm[node_lookup[node]]) + 1):  # cant select more than |v| nodes at that vertex , thus min(k,|v|)
                # for average the extreme cases always exist that is right and left subtree are always chosen

                sums = Decimal(0)
                for r in range(max(0, p - countterm[node_lookup[y]]), min(countterm[node_lookup[x]], p) + 1):  # goes over range for possible r values such that r+l =k
                    l = p - r  # l = k-r
                    if (l==0 and r==p):
                        sums+= Decimal(beta[node_lookup[x]][p] + choosearr[countterm[node_lookup[x]]][p]*v_x)
                    elif(r==0 and l==p):
                        sums += Decimal(beta[node_lookup[y]][p] + choosearr[countterm[node_lookup[y]]][p]*v_y)
                    else:
                        sums +=Decimal((choosearr[countterm[node_lookup[y]]][l]*(beta[node_lookup[x]][r])) \
                                + (choosearr[countterm[node_lookup[x]]][r]*(beta[node_lookup[y]][l]))\
                                + (choosearr[countterm[node_lookup[y]]][l]*choosearr[countterm[node_lookup[x]]][r]*lamb))

                beta[node_lookup[node]][p] =Decimal(sums) #this is beta

                alpha[node_lookup[node]][p] =Decimal(sums/choosearr[countterm[node_lookup[node]]][p]) #this is alpha

    alpha[node_lookup[tree.root]][0] = Decimal(0)
    beta[node_lookup[tree.root]][0] = Decimal(0)

    for node in post_order_nodes:
        avgvals.append(float(alpha[node_lookup[node]][k]))

    #for node in post_order_nodes:
        #avgvals.append(arr1[node_lookup[node]][k])


    return (alpha[node_lookup[tree.root]], beta[node_lookup[tree.root]],beta,alpha,avgvals)