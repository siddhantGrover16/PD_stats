import numpy as np
from decimal import *

def Getdir_var(tree,k, alpha, beta, gamma, post_order_nodes, node_lookup, comb, countterm):
    """
    :param tree: input tree
    :param k: user selction for choice of k leaves
    :param alpha: a n*k matrix consisting of average PD rooted at node n for all choices ranging from 1 to k
    :param beta: a n*k matrix consisting of sum_PD rooted at node n for all choices ranging from 1 to k
    :param gamma: a n*k matrix consisting of \sum_PD^2 rooted at node n for all choices ranging from 1 to k
    [the following is the preprocessed data on the tree]
    :param post_order_nodes:
    :param node_lookup:
    :param comb:
    :param countterm:
    :return:
    """
    getcontext().prec = 64
    v = len(post_order_nodes)
    k = k
    varvals = []

    arr = np.full((v, k + 1), Decimal(0))

    for node in post_order_nodes:
        if node.is_terminal():
            arr[node_lookup[node]] = Decimal(0)

        else:
            for p in range(min(k, countterm[node_lookup[node]]) + 1):
                a = gamma[node_lookup[node]][p]/comb[countterm[node_lookup[node]]][p]
                b = alpha[node_lookup[node]][p]*alpha[node_lookup[node]][p]
                arr[node_lookup[node]][p] = Decimal((gamma[node_lookup[node]][p]/comb[countterm[node_lookup[node]]][p] -alpha[node_lookup[node]][p]*alpha[node_lookup[node]][p]))

    for node in post_order_nodes:
        varvals.append(float(arr[node_lookup[node]][k]))

    return (arr[node_lookup[tree.root]],arr,varvals)







