
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
        i = i + 1
    for node in treed.postorder_node_iter():
        if (len(node.leaf_nodes())<k+1):
            node.annotations.clear()

    treed.write(
        path="annotated_tree.nex",
        schema="nexus",
    )



