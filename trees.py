from Bio import Phylo
from io import StringIO

def getTreeFromFile(file):
    tree = Phylo.read(file,"newick")
    return tree

def getTreeFromString(string):
    tree = Phylo.read(StringIO(string), "newick")
    return tree