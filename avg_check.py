import math
import time
from Bio import Phylo
from io import StringIO
from avg import GetnewAvgPD
from trees import getTreeFromFile
from preprocessing import preprocessTree

tree = getTreeFromFile("t_11_1000")
tic = time.perf_counter()
x=preprocessTree(tree,1000)
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y=GetnewAvgPD(tree,1000,x[0],x[1],x[2],x[3])
print(y[0][200])
print(y[0][400])
print(y[0][600])
print(y[0][800])
print(y[0][1000])
toc = time.perf_counter()
print(toc-tic)
