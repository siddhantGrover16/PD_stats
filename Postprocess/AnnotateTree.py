
def AnnotateTree(treed,annot,annot_list):
    """
    :param treed: input dendropy tree
    :param annot: the annotation string - "min", "max", "avg"--etc
    :param annot_list: the list of annotations
    :return: an annotated nexus tree
    """
    i = 0
    for node in treed.postorder_node_iter():
        if (annot_list[i] == float('inf') or annot_list[i] == float(0)):
            node.annotations.clear()

        else:
            node.annotations[annot] = annot_list[i]

        i = i + 1

    treed.write(
        path="annotated_tree.nex",
        schema="nexus",
    )



