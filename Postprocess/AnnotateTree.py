
def AnnotateTree(treed,annot,annot_list,k):
    f = open("PDSTAT.txt", "a")
    f.write("Node Label, k, AnnotationSet\n")
    """
    :param treed: input dendropy tree
    :param annot: the annotation string - "min", "max", "avg"--etc
    :param annot_list: the list of annotations
    :return: an annotated nexus tree along with a file that contains the node id, k, annotation set
    
    """
    i = 0
    for node in treed.postorder_node_iter():
        node.annotations[annot] = annot_list[i]
        if node.label == None:
            node.label = "Post"+str(i)
        #node.annotations["i_n"] = 1
        if (annot_list[i] == float('inf') or annot_list[i] == float('-inf') or annot_list[i] == float('0')):
            node.annotations.clear()
           # node.annotations["i_n"]= 0
        i = i + 1


    for node in treed.postorder_node_iter():
        if (len(node.leaf_nodes())<k+1):
            node.annotations.clear()
        if (len(node.leaf_nodes())>=k+1):
            f.write(node.label + ", "+ str(k)+ ", " + str(node.annotations)+ "\n")

    f.close()
          #  node.annotations["i_n"] = 0

    treed.write(
        path="annotated_tree.nex",
        schema="nexus",
    )




