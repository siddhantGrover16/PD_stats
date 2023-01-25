# PD_stats
A project designed with algorithms to compute PD Statisttics- Minimum, Maximum, Average, Variance on a nexus tree

# Requirements
  * Python 3.7
 
# Packages Used
  * dendropy
  * Bio
  * decimal
  * scipy
  * numpy
  * argparse
  
  
the statistic arguements are as follows
-fmin:- min
-fmax:- max
-favg:- avg
-fvar:- var
-fhot:- hotspot measure
 
To run this tool on a nexus tree please follow the steps below
1) place nexus treefile in the same directory as the project
2) run pdstat.py with the following arguemenets:- treename, (int)k, min/max/avg/var
  * example : pdstat.py treename.tre 10 -fm 
    The example above runs pdstat.py on the tree treefile where 10 is the user input and wants to finds the minPD 
3) There will be a file named "annotated_tree_nex" with the clades annotated with the minPD at clade

A test file "t100.tre" is also placed in the directory for usage

Scalability of Algorithms:-
The above algorithms were run on different trees with the number of leaves ranging from 1000 to 10000 taxa with k = n
![image](https://user-images.githubusercontent.com/46168937/213595654-48da5734-dcf1-460d-b7e7-1f0c94bc804b.png)


