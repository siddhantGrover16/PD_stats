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
from dendropy.simulate import treesim
from Preprocess.treePP import getTreeFromFile
from Preprocess.PP import preprocessTree
from Stat_Functions.SumSq import Getsumsq
from Stat_Functions.AvgPD import GetSums_Avg
from Stat_Functions.DirVar import Getdir_var
from Preprocess.treePP import getNexTreeFromFile

"""
Dendropy is used to simulate the birth-death trees used for all the experiments done for the project
"""
def GenerateTree(birthrate,deathrate,number_of_Tax):
    """

    :param birthrate: birth rate
    :param deathrate:  death rate
    :param number_of_Tax: desired number of taxa
    :return: tree and the newick and nexus file of tree
    """
    t = treesim.birth_death_tree(birth_rate=birthrate, death_rate=deathrate, ntax=number_of_Tax)
    t.write(path="t_"+str(number_of_Tax)+"nex",schema="nexus")
    t.write(path="t_" + str(number_of_Tax) + "new", schema="newick")
    return t,"t_"+str(number_of_Tax)+"nex"

def mean_sampling(tree2,k,nruns,mean,track):
 """
 The method below samples and computes the pd of each sampled set
 :param tree2:  original tree with desired taxa (t500_.tre)
 :param nruns: number of sets sampled
 :param mean: mean value by algorithm
 :param k: user input for subset size
 :param track: keeps track of number of runs
 :return:
 """
 mean = float(round(mean,4))
 tree2=tree2
 t2= copy.deepcopy(tree2)
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
     total=total+pd

 out = total/nruns
 if (0.999*mean<=out<=1.001 *mean):
  return out
 else:
     nruns+=track
     #print(track)
     mean_sampling(tree2, k, nruns, mean, track)
     #tree2, k, nruns, mean, rate, track

