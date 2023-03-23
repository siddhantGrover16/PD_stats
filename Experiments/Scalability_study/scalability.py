import time

from Stat_Functions.MaxPD import GetMaxPD
from Stat_Functions.MinPD import GetMinPD
from Preprocess.treePP import getTreeFromFile
from Preprocess.PP import preprocessTree
from Stat_Functions.SumSq import Getsumsq
from Stat_Functions.AvgPD import GetSums_Avg

for i in range(1001):
    if i%100 == 0:
        print(i)

t1 = getTreeFromFile("t_11_1000")
t2 = getTreeFromFile("t_10_2000")
t3 = getTreeFromFile("t_9_3000")
t4 = getTreeFromFile("t_8_4000")
t5 = getTreeFromFile("t_7_5000")
t6 = getTreeFromFile("t_6_6000")
t7 = getTreeFromFile("t_5_7000")
t8 = getTreeFromFile("t_4_8000")
t9 = getTreeFromFile("t_3_9000")
t10 = getTreeFromFile("t_2_10000")

print("-------------prepro times--------------")
k=1000
tic = time.perf_counter()
x1=preprocessTree(t1,k)
toc = time.perf_counter()
print(toc-tic)

k=2000
tic = time.perf_counter()
x2=preprocessTree(t2,k)
toc = time.perf_counter()
print(toc-tic)

k=3000
tic = time.perf_counter()
x3=preprocessTree(t3,k)
toc = time.perf_counter()
print(toc-tic)

k=4000
tic = time.perf_counter()
x4=preprocessTree(t4,k)
toc = time.perf_counter()
print(toc-tic)

k=5000
tic = time.perf_counter()
x5=preprocessTree(t5,k)
toc = time.perf_counter()
print(toc-tic)

k=6000
tic = time.perf_counter()
x6=preprocessTree(t6,k)
toc = time.perf_counter()
print(toc-tic)

k=7000
tic = time.perf_counter()
x7=preprocessTree(t7,k)
toc = time.perf_counter()
print(toc-tic)

k=8000
tic = time.perf_counter()
x8=preprocessTree(t8,k)
toc = time.perf_counter()
print(toc-tic)

k=9000
tic = time.perf_counter()
x9=preprocessTree(t9,k)
toc = time.perf_counter()
print(toc-tic)

k=10000
tic = time.perf_counter()
x10=preprocessTree(t10,k)
toc = time.perf_counter()
print(toc-tic)




print("-------------avg times--------------")
tic = time.perf_counter()
y1=GetSums_Avg(t1,1000,x1[0],x1[1],x1[2],x1[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y2=GetSums_Avg(t2,2000,x2[0],x2[1],x2[2],x2[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y3=GetSums_Avg(t3,3000,x3[0],x3[1],x3[2],x3[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y4=GetSums_Avg(t4,4000,x4[0],x4[1],x4[2],x4[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y5=GetSums_Avg(t5,5000,x5[0],x5[1],x5[2],x5[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y6=GetSums_Avg(t6,6000,x6[0],x6[1],x6[2],x6[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y7=GetSums_Avg(t7,7000,x7[0],x7[1],x7[2],x7[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y8=GetSums_Avg(t8,8000,x8[0],x8[1],x8[2],x8[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y9=GetSums_Avg(t9,9000,x9[0],x9[1],x9[2],x9[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
y10=GetSums_Avg(t10,10000,x10[0],x10[1],x10[2],x10[3])
toc = time.perf_counter()
print(toc-tic)

print("-------------sumsq times--------------")#the computation of sum square diretcly facilitates the variance and all the ther operations are on precomputed data(timed as well) which are done in constant time
tic = time.perf_counter()
Getsumsq(t1,1000,y1[2],x1[0],x1[1],x1[2],x1[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t2,2000,y2[2],x2[0],x2[1],x2[2],x2[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t3,3000,y3[2],x3[0],x3[1],x3[2],x3[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t4,4000,y4[2],x4[0],x4[1],x4[2],x4[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t5,5000,y5[2],x5[0],x5[1],x5[2],x5[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t6,6000,y6[2],x6[0],x6[1],x6[2],x6[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t7,7000,y7[2],x7[0],x7[1],x7[2],x7[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t8,8000,y8[2],x8[0],x8[1],x8[2],x8[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t9,9000,y9[2],x9[0],x9[1],x9[2],x9[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
Getsumsq(t10,10000,y10[2],x10[0],x10[1],x10[2],x10[3])
toc = time.perf_counter()
print(toc-tic)

print("-------------min times--------------")
tic = time.perf_counter()
GetMinPD(t1,1000,x1[0],x1[1],x1[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t2,2000,x2[0],x2[1],x2[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t3,3000,x3[0],x3[1],x3[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t4,4000,x4[0],x4[1],x4[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t5,5000,x5[0],x5[1],x5[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t6,6000,x6[0],x6[1],x6[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t7,7000,x7[0],x7[1],x7[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t8,8000,x8[0],x8[1],x8[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t9,9000,x9[0],x9[1],x9[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMinPD(t10,10000,x10[0],x10[1],x10[3])
toc = time.perf_counter()
print(toc-tic)

print("-------------max times--------------")
tic = time.perf_counter()
GetMaxPD(t1,1000,x1[0],x1[1],x1[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t2,2000,x2[0],x2[1],x2[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t3,3000,x3[0],x3[1],x3[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t4,4000,x4[0],x4[1],x4[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t5,5000,x5[0],x5[1],x5[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t6,6000,x6[0],x6[1],x6[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t7,7000,x7[0],x7[1],x7[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t8,8000,x8[0],x8[1],x8[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t9,9000,x9[0],x9[1],x9[3])
toc = time.perf_counter()
print(toc-tic)

tic = time.perf_counter()
GetMaxPD(t10,10000,x10[0],x10[1],x10[3])
toc = time.perf_counter()
print(toc-tic)












