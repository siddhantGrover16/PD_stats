import time
from MaxPD import GetMaxPD
from trees import getTreeFromFile
from preprocessing import preprocessTree

tree = getTreeFromFile("t_11_1000")
tic = time.perf_counter()
x=preprocessTree(tree,1000)
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y=GetMaxPD(tree,1000,x[0],x[1],x[3])
print(y[200])
print(y[400])
print(y[600])
print(y[800])
print(y[1000])
toc = time.perf_counter()
print(toc-tic)