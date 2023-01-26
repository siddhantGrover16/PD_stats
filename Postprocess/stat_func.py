from Postprocess.AnnotateTree import AnnotateTree
from Preprocess.PP import preprocessTree
from Preprocess.treePP import getTree, getDendrotree
from Stat_Functions.AvgPD import GetSums_Avg
from Stat_Functions.DirVar import Getdir_var
from Stat_Functions.HotPD import Getdiv_hot
from Stat_Functions.MaxPD import GetMaxPD
from Stat_Functions.MinPD import GetMinPD
from Stat_Functions.SumSq import Getsumsq


def base(tree,k):
    t1 = getTree(tree, "nexus", True)
    o_tree = getDendrotree(tree, "nexus")
    pp = preprocessTree(t1, k)
    return(t1,o_tree,pp,k)

def MinPD(t1,o_tree,pp,k):
    minPD = GetMinPD(t1, k, pp[0], pp[1], pp[3])
    minvals = minPD[3]
    AnnotateTree(o_tree, "min", minvals,k)

def MaxPD(t1,o_tree,pp,k):
    maxPD = GetMaxPD(t1, k, pp[0], pp[1], pp[3])
    maxvals = maxPD[3]
    AnnotateTree(o_tree, "max", maxvals,k)

def AvgPD(t1,o_tree,pp,k):
    avg = GetSums_Avg(t1, k, pp[0], pp[1], pp[2], pp[3])
    avgvals = avg[4]
    AnnotateTree(o_tree, "avg", avgvals,k)

def VarPD(t1,o_tree,pp,k):
    avg = GetSums_Avg(t1, k, pp[0], pp[1], pp[2], pp[3])
    z = Getsumsq(t1, k, avg[2], pp[0], pp[1], pp[2], pp[3])
    var = Getdir_var(t1, k, avg[3], avg[2], z[1], pp[0], pp[1], pp[2], pp[3])
    varvals = var[2]
    AnnotateTree(o_tree, "var", varvals,k)

def HotPD(t1,o_tree,pp,k):
    max = GetMaxPD(t1, k, pp[0], pp[1], pp[3])
    avg = GetSums_Avg(t1, k, pp[0], pp[1], pp[2],pp[3])
    hot = Getdiv_hot(t1, k, avg[3], max[2], pp[0], pp[1], pp[3])
    hotvals = hot[2]
    AnnotateTree(o_tree, "hot", hotvals,k)

def AllPD(t1,o_tree,pp,k):

    MinPD(t1, o_tree, pp, k)
    MaxPD(t1, o_tree, pp, k)
    AvgPD(t1, o_tree, pp, k)
    VarPD(t1, o_tree, pp, k)
    HotPD(t1, o_tree, pp, k)




