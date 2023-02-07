# Phylogenetic Diversity Statistics for All Clades in a Phylogeny
A project designed with algorithms to compute Phylogenetic Diversity (PD) Statistics - Minimum, Maximum, Average, Variance, and HotSpots on a phylogenetic tree.

Grover, S., Markin, A., Anderson, T.K., and Eulenstein, O. (in review). Phylogenetic Diversity Statistics for All Clades in a Phylogeny.

Abstract: The classic quantitative measure of phylogenetic diversity, PD, has been used to address problems in conservation biology, microbial ecology, and evolutionary biology. PD is the minimum total length of the branches in a phylogeny required to cover a specified set of taxa on the phylogeny. A general goal in the application of PD has been identifying taxa that maximize PD on a given phylogeny; this has been mirrored in active research to develop efficient algorithms for the problem. Other descriptive statistics, such as the minimum PD, average PD, and standard deviation of PD, can provide invaluable insight into the distribution of PD across a phylogeny. However, there has been limited or no research on computing these statistics, especially when required for each clade in a phylogeny, enabling direct comparisons of PD between clades. We introduce efficient algorithms for computing PD and the associated descriptive statistics for a given phylogeny and each of its clades. In simulation studies, we demonstrate the ability of our algorithms to analyze large-scale phylogenies with applications in ecology and evolutionary biology

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
-fall:- computes all the statistics above
 
To run this tool on a nexus tree please follow the steps below:
1) place nexus binary treefile in the same directory as the project
2) run pdstat.py with the following arguements: -treename, (int)k, min/max/avg/var
  * example : pdstat.py treename.tre 10 -fmax 
    The example above runs pdstat.py on the tree treefile and finds the maxPD for k=10 taxa in the tree.
    
    Note:- the user can run multiple functions on the the treefile, eg. pdstat.py treename.tre 10 -fmax -fmin
3) The output of the above command will be a file named "annotated_tree_nex" with the clades annotated with the maxPD at clades with k=10. This tree file can be opened with FigTree, and the nodes annotated using the "node labels" or "node shapes" menu.

Input test files with different numbers of input taxa are also present in the directory (t50.tre,t100.tre,t200.tre,t300.tre,t400.tre).

Scalability of Algorithms:-
The above algorithms were run on different trees with the number of leaves ranging from 1000 to 10000 taxa with k = n
![image](https://user-images.githubusercontent.com/46168937/213595654-48da5734-dcf1-460d-b7e7-1f0c94bc804b.png)


