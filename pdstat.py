import argparse
from Postprocess.stat_func import MinPD
from Postprocess.stat_func import MaxPD
from Postprocess.stat_func import AvgPD
from Postprocess.stat_func import VarPD
from Postprocess.stat_func import HotPD
from Postprocess.stat_func import AllPD
from Postprocess.stat_func import base

def Main():

 ap = argparse.ArgumentParser(description='Compute PD Statistics on a given nexus tree')
 ap.add_argument("-fmin", help="Computes Min PD", action="store_true")
 ap.add_argument("-fmax", help="Computes Max PD", action="store_true")
 ap.add_argument("-favg", help="Computes Avg PD", action="store_true")
 ap.add_argument("-fvar", help="Computes Var PD", action="store_true")
 ap.add_argument("-fhot", help="Computes PD hotspot", action="store_true")
 ap.add_argument("-fall", help="Computes All PD statistics", action="store_true")


 ap.add_argument("filename",help = "nexus file treename",type=str)
 ap.add_argument("kval", help="input k", type=int)

 args = ap.parse_args()

 init=base(args.filename,args.kval)
#args.num1,args.num2
 #tree = getTree(args.filename, "nexus", True)
 #o_tree = getDendrotree(args.filename, "nexus")
 #pp = preprocessTree(tree,args.kval)

 if args.fmin:
     MinPD(init[0],init[1],init[2],init[3])
   #minPD = GetMinPD(tree, args.kval, pp[0], pp[1], pp[3])
   #minvals = minPD[3]
   #AnnotateTree(o_tree, "min", minvals, 10)

 if args.fmax:
     MaxPD(init[0], init[1], init[2], init[3])

   #maxPD = GetMaxPD(tree, args.kval, pp[0], pp[1], pp[3])
   #maxvals = maxPD[3]
   #AnnotateTree(o_tree, "max", maxvals, 10)

 if args.favg:
     AvgPD(init[0], init[1], init[2], init[3])
   #avg=GetSums_Avg(tree,args.kval,pp[0],pp[1],pp[2],pp[3])
   #avgvals = avg[4]
   #AnnotateTree(o_tree, "avg", avgvals, 10)

 if args.fvar:
     VarPD(init[0], init[1], init[2], init[3])
   #avg = GetSums_Avg(tree, args.kval, pp[0], pp[1], pp[2], pp[3])
   #z = Getsumsq(tree, args.kval, avg[2], pp[0], pp[1], pp[2], pp[3])
   #var = Getdir_var(tree, args.kval, avg[3], avg[2], z[1], pp[0], pp[1], pp[2], pp[3])
   #varvals = var[2]
   #AnnotateTree(o_tree, "var", varvals, 10)

 if args.fhot:
     HotPD(init[0], init[1], init[2], init[3])

 if args.fall:
     AllPD(init[0], init[1], init[2], init[3])


if __name__ == '__main__':
  Main()











