import numpy as np
from decimal import *

def Getsumsq(tree, k, beta, post_order_nodes, node_lookup, comb, countterm):
    """
    Given a Bio.Phylo tree and the preprocessed data, computes the \sum_PD^2 PD statistics on a tree
    :param tree: input tree
    :param k: user selection for choice of k leaves
    :param beta: a n*k matrix consisting of sum_PD rooted at node n for all choices ranging from 1 to k
    [the following is the preprocessed data on the tree]
    :param post_order_nodes:
    :param node_lookup:
    :param comb:
    :param countterm:
    :return: a n*k matrix with \sum_PD^2 rooted at node n for all choices ranging from 1 to k
        """

    getcontext().prec =64
    v = len(post_order_nodes)
    k = k

    #create array of (|v|k) size with all values = inf
    arr = np.full((v, k + 1), Decimal(0))
    arr1 = np.full((v, k + 1), Decimal(0))

    arr[:, 0] = Decimal(0) #base case 1- for k=0, d(v,k)=0 regardless of vertex
    arr1[:, 0] = Decimal(0)


    #calculate D(v,k) for all nodes
    for node in post_order_nodes:
        if node.is_terminal():
            arr[node_lookup[node]] = Decimal(0)
            arr1[node_lookup[node]] = Decimal(0)

        #base case 2- D(v,k) = 0 for all leaf nodes
            # D(v.k) for all k = 0
        else:
            #recursive case
            x = node.clades[0] #right tree
            y = node.clades[1] #left tree
            v_x = Decimal(tree.distance(node, x)) #edge (v,x)
            v_y = Decimal(tree.distance(node, y)) #edge (v,y)
            lamb = Decimal(v_x+v_y)
            lamb2= Decimal(lamb*lamb)
            for p in range(min(k,countterm[node_lookup[node]])+1): # cant select more than |v| nodes at that vertex , thus min(k,|v|)
                #for average the extreme cases always exist that is right and left subtree are always chosen
               # print ("p =", p)
                sums_x = Decimal(0)
                sums_y = Decimal(0)
                sums = Decimal(0)
                for r in range(max(0, p - countterm[node_lookup[y]]), min(countterm[node_lookup[x]],p) + 1):  # goes over range for possible r values such that r+l =k
                    l = p - r  # l = k-r
                    if (l==0 and r==p):
                        sums += Decimal(arr[node_lookup[x]][p] + comb[countterm[node_lookup[x]]][p] * (v_x) * (v_x) + 2 * (
                            v_x) * (beta[node_lookup[x]][p]))
                    elif(r==0 and l==p):
                        sums += Decimal(arr[node_lookup[y]][p] + comb[countterm[node_lookup[y]]][p] * (v_y) * (v_y) + 2 * (
                            v_y) * (beta[node_lookup[y]][p]))
                    else:
                        sums += Decimal(comb[countterm[node_lookup[y]]][l]*(arr[node_lookup[x]][r]) \
                                + comb[countterm[node_lookup[x]]][r]*(arr[node_lookup[y]][l]) \
                                + comb[countterm[node_lookup[y]]][l] *comb[countterm[node_lookup[x]]][r]*(lamb2) +\
                                2*(beta[node_lookup[x]][r]*beta[node_lookup[y]][l] +\
                                   lamb*(comb[countterm[node_lookup[y]]][l]\
                                         *beta[node_lookup[x]][r]+comb[countterm[node_lookup[x]]][r]*beta[node_lookup[y]][l])))

                #sums += arr[node_lookup[x]][p] + comb[countterm[node_lookup[x]]][p]*(v_x)*(v_x) + 2*(v_x)*(beta[node_lookup[x]][p])
                #sums += arr[node_lookup[y]][p] + comb[countterm[node_lookup[y]]][p] * (v_y) * (v_y) + 2 * (v_y) * (beta[node_lookup[y]][p])


               # print(sum_avg/(choose(node.count_terminals(), p)))
                arr[node_lookup[node]][p] = sums
                arr[node_lookup[node]][0] = 0


    return (arr[node_lookup[tree.root]],arr)