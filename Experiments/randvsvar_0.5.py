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


random.seed(1)
t2= copy.deepcopy(tree2)

"""
Precalculated the mean and variance below using AvgPd,SumSq,DirVar Algorithm
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

var50= 78.365
var100=94.014
var150=95.548
var200=91.329
var250=82.947
var300=71.251
var350=56.885
var400=40.221
var450=21.329


def var_sampling(tree2,t2,nruns,mean,var,k,track):
    """
    The method below samples and computes the pd of each sampled set, and then computes the variance manually at the end of
    the runs

   :param tree2:  original tree with 500 taxa (t500_.tre)
    :param t2: deepcopy of tree
    :param mynumber: number of sets sampled
    :param mean: mean value for given k
    :param var: variance for given k
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
    total=0
    onum=0
    for x in t2.leaf_node_iter():
        l2.append(x.taxon.label)

    for j in range(nruns):
        t2 = copy.deepcopy(tree2)
        random.shuffle(l2)
        l2_2 = l2[0:k]
        pd = faith_pd(u_counts, l2_2, mytree)
        diff = pd-mean
        temp = diff*diff
        onum=onum+temp
    out = onum/nruns

    if (0.995*var <= out <= 1.005*var):
        # print("done", out)
        return out
    else:
        track = track + 200
        #print(track,out)
        var_sampling(tree2, t2, nruns + track, mean,var, k, track)



"""
The method above randomly samples taxa till the variance converges to a 0.5% difference from the actual value. 
This process is repeated and the time is evaluated for k ranging from 50 to 450
"""

st = time.time()
var_sampling(tree2,t2,200,mean50,var50,50,100)
et= time.time()
print(et-st)


st = time.time()
var_sampling(tree2,t2,200,mean100,var100,100,100)
et= time.time()
print(et-st)

#
st = time.time()
var_sampling(tree2,t2,200,mean150,var150,150,100)
et= time.time()
print(et-st)


st = time.time()
var_sampling(tree2,t2,200,mean200,var200,200,100)
et= time.time()
print(et-st)

st = time.time()
var_sampling(tree2,t2,200,mean250,var250,250,100)
et= time.time()
print(et-st)


st = time.time()
var_sampling(tree2,t2,200,mean300,var300,300,100)
et= time.time()
print(et-st)


st = time.time()
var_sampling(tree2,t2,200,mean350,var350,350,100)
et= time.time()
print(et-st)


st = time.time()
var_sampling(tree2,t2,200,mean400,var400,400,100)
et= time.time()
print(et-st)



st = time.time()
var_sampling(tree2,t2,200,mean450,var450,450,100)
et= time.time()
print(et-st)







