import argparse
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


import sys




def Main():

 ap = argparse.ArgumentParser(description='Compute PD Statistics on a given nexus tree')
 stat_g = ap.add_mutually_exclusive_group()
 stat_g.add_argument("-fm", "--fmin", help="Computes Min PD", action="store_true")
 stat_g.add_argument("-fx", "--fmax", help="Computes Max PD", action="store_true")
 stat_g.add_argument("-fa", "--favg", help="Computes Avg PD", action="store_true")
 stat_g.add_argument("-fv", "--fvar", help="Computes Var PD", action="store_true")

 ap.add_argument("filename",help = "nexus file treename",type=str)
 ap.add_argument("kval", help="input k", type=int)

 args = ap.parse_args()
#args.num1,args.num2
 tree = getTree(args.filename, "nexus", True)
 o_tree = getDendrotree(args.filename, "nexus")
 pp = preprocessTree(tree,ap.kval)


 if args.fmin:
   minPD = GetMinPD(tree, args.kval, pp[0], pp[1], pp[3])
   minvals = minPD[3]
   AnnotateTree(o_tree, "min", minvals, 10)

 elif args.fmax:
   maxPD = GetMaxPD(tree, args.kval, pp[0], pp[1], pp[3])
   maxvals = maxPD[3]
   AnnotateTree(o_tree, "max", maxvals, 10)

 elif args.favg:
   avg=GetSums_Avg(tree,args.kval,pp[0],pp[1],pp[2],pp[3])
   avgvals = avg[4]
   AnnotateTree(o_tree, "avg", avgvals, 10)

 elif args.fvar:
   avg = GetSums_Avg(tree, args.kval, pp[0], pp[1], pp[2], pp[3])
   z = Getsumsq(tree, args.kval, avg[2], pp[0], pp[1], pp[2], pp[3])
   var = Getdir_var(tree, args.kval, avg[3], avg[2], z[1], pp[0], pp[1], pp[2], pp[3])
   varvals = var[2]
   AnnotateTree(o_tree, "var", varvals, 10)

 else:
   print("Error:Requires an argument to perform an action")


if __name__ == '__main__':
  Main()











