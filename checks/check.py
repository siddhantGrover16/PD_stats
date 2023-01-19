import math
import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
import dendropy
from Stat_Functions.MaxPD import GetMaxPD
from Stat_Functions.MinPD import GetMinPD
from Stat_Functions.SumSq import Getsumsq
from Stat_Functions.AvgPD import GetSums_Avg
from Stat_Functions.DirVar import Getdir_var
from Preprocess.treePP import binarize_tree
from Preprocess.treePP import getTreeFromString
from Preprocess.treePP import getTree
from Preprocess.treePP import getDendrotree
from Postprocess.GetRatio import GetRatio
from Postprocess.AnnotateTree import AnnotateTree
from Preprocess.PP import preprocessTree
import dendropy

tree = getTree("t100.tre","nexus",True)
o_tree = getDendrotree("t100.tre","nexus")

pp=preprocessTree(tree,3)

minPD=GetMinPD(tree,3,pp[0],pp[1],pp[3])
maxPD=GetMaxPD(tree,3,pp[0],pp[1],pp[3])
avg=GetSums_Avg(tree,3,pp[0],pp[1],pp[2],pp[3])
z = Getsumsq(tree,3,avg[2],pp[0],pp[1],pp[2],pp[3])
var = Getdir_var(tree,3,avg[3],avg[2],z[1],pp[0],pp[1],pp[2],pp[3])

minvals = minPD[3]
maxvals = maxPD[3]
avgvals = avg[4]
varvals = var[2]
res = GetRatio(avgvals,maxvals)

AnnotateTree(o_tree,"a_m",res,10)
AnnotateTree(o_tree,"min",minvals,10)

