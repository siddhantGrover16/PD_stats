import math
import time
import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
from io import StringIO
import dendropy
from dendropy.simulate import treesim
import copy
import random
from io import StringIO
from skbio import TreeNode
from skbio.diversity.alpha import faith_pd
from Preprocess.treePP import getNexTreeFromFile

tree2 = dendropy.Tree.get(path="t500_.tre", schema="newick")
tree = getNexTreeFromFile("t500.tre")
t2= copy.deepcopy(tree2)

"""
Precalculated the mean below using AvgPd Algorithm
"""
mean50  = 208.625
mean100 = 314.130
mean150 = 392.026
mean200 = 455.272
mean250= 509.146
mean300  = 556.412
mean350 =  598.723
mean400 = 637.132
mean450 = 672.337


random.seed(7)# random seed to recompute results


def mean_sampling(tree2,t2,nruns,mean,k,track):
    """
    The method below samples and computes the pd of each sampled set
    :param tree2:  original tree with 500 taxa (t500_.tre)
    :param t2: deepcopy of tree
    :param mynumber: number of sets sampled
    :param mean: mean value by algorithm
    :param k: user input for subset size
    :param track: keeps track of number of runs
    :return:
    """
    tree2=tree2
    t2=t2
    x = tree2.as_string(schema="newick")
    u_counts = [1] * k
    mytree = TreeNode.read(StringIO(x))
    l2 = []
    total = 0
    for x in t2.leaf_node_iter():
        l2.append(x.taxon.label)

    for j in range(nruns):
        t2 = copy.deepcopy(tree2)
        random.shuffle(l2)
        l2_2 = l2[0:k]
        pd = faith_pd(u_counts, l2_2, mytree)
        #print(j,pd)
        total=total+round(pd,3)

    out = total/nruns
    if (0.995*mean<=out<=1.005 *mean):
    #print("done", out)
     return out
    else:
     track=track+track
     #print(track,out)
     mean_sampling(tree2, t2, nruns+track, mean,k,track)

"""
The method above randomly samples taxa till the average converges to a 0.5% difference from the actual value. 
This process is repeated and the time is evaluated for k ranging from 50 to 450
"""

st = time.time()
mean_sampling(tree2,t2,200,mean50,50,100)
et= time.time()
print(et-st)


st = time.time()
mean_sampling(tree2,t2,200,mean100,100,100)
et= time.time()
print(et-st)


st = time.time()
mean_sampling(tree2,t2,200,mean150,150,100)
et= time.time()
print(et-st)

st = time.time()
mean_sampling(tree2,t2,200,mean200,200,100)
et= time.time()
print(et-st)


st = time.time()
mean_sampling(tree2,t2,200,mean250,250,100)
et= time.time()
print(et-st)


st = time.time()
mean_sampling(tree2,t2,200,mean300,300,100)
et= time.time()
print(et-st)


st = time.time()
mean_sampling(tree2,t2,200,mean350,350,100)
et= time.time()
print(et-st)


st = time.time()
mean_sampling(tree2,t2,200,mean400,400,100)
et= time.time()
print(et-st)



st = time.time()
mean_sampling(tree2,t2,200,mean450,450,100)
et= time.time()
print(et-st)







