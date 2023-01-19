import math
import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
import time
import dendropy
from Stat_Functions.MaxPD import GetMaxPD
from Stat_Functions.MinPD import GetMinPD
from Stat_Functions.SumSq import Getsumsq
from Stat_Functions.AvgPD import GetSums_Avg
from Stat_Functions.DirVar import Getdir_var
from Preprocess.treePP import binarize_tree
from Preprocess.treePP import getTreeFromString
from Preprocess.PP import preprocessTree
import dendropy

def AnnotateTree(treed,annot,annot_list,num_chil):
    """
    :param treed: input dendropy tree
    :param annot: the annotation string - "min", "max", "avg"--etc
    :param annot_list: the list of annotations
    :param num_child: restrict the annotations on nodes with atleast num_child
    :return: an annotated nexus tree
    """
    i=0
    for node in treed.postorder_node_iter():
        node.annotations[annot] = annot_list[i]
        i=i+1

    for node in treed.postorder_node_iter():
        if len(node.leaf_nodes()) < num_chil:
            node.annotations.clear()

    for node in treed.postorder_node_iter():
        if len(node.leaf_nodes()) < num_chil:
            node.annotations["node_size"] = 0.0
        else:
            node.annotations["node_size"] = 11

    treed.write(
        path="annotated_tree.nex",
        schema="nexus",
    )



