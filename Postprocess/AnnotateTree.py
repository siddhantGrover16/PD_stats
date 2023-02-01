
def AnnotateTree(treed,annot,annot_list,k):
    """
    :param treed: input dendropy tree
    :param annot: the annotation string - "min", "max", "avg"--etc
    :param annot_list: the list of annotations
    :return: an annotated nexus tree
    """
    i = 0
    for node in treed.postorder_node_iter():
        node.annotations[annot] = annot_list[i]
        #node.annotations["i_n"] = 1
        if (annot_list[i] == float('inf') or annot_list[i] == float('-inf') or annot_list[i] == float('0')):
            node.annotations.clear()
           # node.annotations["i_n"]= 0
        i = i + 1
    for node in treed.postorder_node_iter():
        if (len(node.leaf_nodes())<k+1):
            node.annotations.clear()
          #  node.annotations["i_n"] = 0

    treed.write(
        path="annotated_tree.nex",
        schema="nexus",
    )



