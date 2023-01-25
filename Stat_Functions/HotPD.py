import numpy as np
from decimal import *

def Getdiv_hot(tree,k, avg, max, post_order_nodes, node_lookup, countterm):
    """
    :param tree: input tree
    :param k: user selction for choice of k leaves
    :param avg: a n*k matrix consisting of average PD rooted at node n for all choices ranging from 1 to k
    :param max: a n*k matrix consisting of max_PD rooted at node n for all choices ranging from 1 to k
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
    hotvals = []

    arr = np.full((v, k + 1), Decimal(0))

    for node in post_order_nodes:
        if node.is_terminal():
            arr[node_lookup[node]] = Decimal(0)

        else:
            for p in range(min(k, countterm[node_lookup[node]]) + 1):
                if max[node_lookup[node]][p] == 0:
                    arr[node_lookup[node]][p] = 0
                else:
                    arr[node_lookup[node]][p] = avg[node_lookup[node]][p]/max[node_lookup[node]][p]

    for node in post_order_nodes:
        hotvals.append(float(arr[node_lookup[node]][k]))

    return (arr[node_lookup[tree.root]],arr,hotvals)