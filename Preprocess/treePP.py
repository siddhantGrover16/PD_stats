import dendropy
from Bio import Phylo
from io import StringIO
from dendropy import Tree, Node

def getTree(path,schema,isbinary):
    """
    :param path: takes in string of tree file absolute path.
    :param schema:  takes input format of the tree resticted to “fasta”, “newick”, “nexus”, “nexml”, or “phylip”.
    :param isbinary: takes in boolean value indicating if the input is a binary tree
    :return: a tree
    """
    dendro = getDendrotree(path,schema)
    if (isbinary == False):
        bin_tree=binarize_tree(dendro,0)
    else:
        bin_tree =dendro

    tree = getTreeFromString(bin_tree.as_string(schema="newick"))
    return tree

def getDendrotree(path,schema):
    tree = dendropy.Tree.get(path=path, schema=schema)
    return tree

def getTreeFromString(string):
    tree = Phylo.read(StringIO(string),"newick")
    return tree

def getNwkTreeFromFile(file):
    """
    Gets Newick Tree from file
    :param newick file: Newick tree with edge lengths.
    """
    tree = Phylo.read(file,"newick")
    return tree

def getNexTreeFromFile(file):
    """
    Gets Nexus Tree from file
    :param nexus file: Nexus tree with edge lengths.
    """
    tree = Phylo.read(file, "nexus")
    return tree

def binarize_tree(tree: Tree, edge_length=0):
    """
    Adds/removes nodes from the tree to make it fully binary (added edges will have length 'edge_length')
    :param tree: Dendropy tree to be made bifurcating.
    """

    # First suppress unifurcations.
    tree.suppress_unifurcations()

    # Now binarize multifurcations.
    for node in tree.postorder_node_iter():
        assert isinstance(node, Node)
        if node.child_nodes() and len(node.child_nodes()) > 2:
            num_children = len(node.child_nodes())
            children = node.child_nodes()
            interim_node = node
            # Creates a caterpillar structure with children on the left of the trunk:
            for child_ind in range(len(children) - 2):
                new_node = Node(edge_length=edge_length)
                interim_node.set_child_nodes([children[child_ind], new_node])
                interim_node = new_node
            interim_node.set_child_nodes(children[num_children - 2:])



