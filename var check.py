import math
from MinPD import GetMinPD
from MaxPD import GetMaxPD
from avg import GetnewAvgPD
from sum_sq import Getsumsq
from trees import getTreeFromFile
from preprocessing import preprocessTree
from scipy.special import comb as choose

k=1000
tree = getTreeFromFile("t_11_1000")
x = preprocessTree(tree,k)
minn = GetMinPD(tree,1000,x[0],x[1],x[3])
maxx = GetMaxPD(tree,1000,x[0],x[1],x[3])
y = GetnewAvgPD(tree,k,x[0],x[1],x[2],x[3])
z = Getsumsq(tree,k,y[2],x[0],x[1],x[2],x[3])

def getvar(k):
    print("new way for k = ", k)
    nck = choose(1000, k)
    var = abs((z[k] / nck) - (y[0][k] * y[0][k]))
    #print("min :", minn[k])
    #print("avg : ", y[0][k])
    #print("max :", maxx[k])
    print("var: ", var)
    std8 = math.sqrt(var)
    print("std : ", std8)
    #print("range : ", y[0][k] - std8, " to", y[0][k] + std8)

getvar(25)
getvar(50)
getvar(100)
getvar(200)
getvar(400)
getvar(600)
getvar(800)
getvar(1000)




