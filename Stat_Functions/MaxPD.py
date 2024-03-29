import numpy as np
from decimal import *

def GetMaxPD(tree, k,post_order_nodes,node_lookup,countterm):
    """
    Given a Bio.Phylo tree and the preprocessed data, computes the MinPD statistics on a tree
    :param tree: input tree
    :param k: user selction for choice of k leaves
    [the following is the preprocessed data on the tree] Note:-max preprocessing doesnt involve computation of pascal triangel
    :param post_order_nodes:
    :param node_lookup:
    :param countterm:
    :return: a n*k matrix with maximum PD rooted at node n for all choices ranging from 1 to k
    """
    getcontext().prec = 64
    if(k==0):
        return Decimal(0)

    v = len(post_order_nodes)
    arr = np.full((v, k + 1), Decimal('Infinity'))
    # taxa_arr = np.full((v,k+1),set(),set)
    arr[:, 0] = Decimal(0)  # base case 1- for k=0, d(v,k)=0 regardless of vertex
    getcontext().prec = 64
    maxvals=[]

    #calculate D(v,k) for all nodes
    for node in post_order_nodes:
        #base case 2- D(v,k) = 0 for all leaf nodes
        if node.is_terminal():
            arr[node_lookup[node]] =Decimal(0)
           # for i in range(k+1):
              #  taxa_arr[node_lookup[node]][1:] = set([node])
            # D(v.k) for all k = 0
        else:
            #recursive case
            x = node.clades[0] #right tree
            y = node.clades[1] #left tree
            v_x =Decimal(tree.distance(node, x)) #edge (v,x)
            v_y =Decimal(tree.distance(node, y)) #edge (v,y)
            for p in range(min(k, node.count_terminals()) + 1): # cant select more than |v| nodes at that vertex , thus min(k,|v|)
                mindist = Decimal('-Infinity')
                #mintaxa = set()
                for r in range(max(0, p - y.count_terminals()), min(x.count_terminals(), p) + 1): #goes over range for possible r values such that r+l =k
                    l = p - r # l = k-r
                    dist = Decimal(arr[node_lookup[x]][r] + arr[node_lookup[y]][l] + v_x * (min(r, 1)) + v_y * (min(l, 1))) #if,r is 0, wont chosose (v,x), same for l=0
                    if dist > mindist:
                        mindist = dist
                arr[node_lookup[node]][p] = mindist # fill d(v,k) with value

    arr[node_lookup[tree.root]][0] = Decimal(0)

    for node in post_order_nodes:
        maxvals.append(float(arr[node_lookup[node]][k]))

    return (arr[node_lookup[tree.root]],tree,arr,maxvals)
