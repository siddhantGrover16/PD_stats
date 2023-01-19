
def GetRatio(list1,list2):
    """
    :param list1: numerator list
    :param list2: denominator list
    :return: a list with num/den values at every given index
    """
    res = [0 for element in range(len(list1))]
    i = 0
    while i < len(list1):
        if list2[i] != 0:
            res[i] = list1[i] / list2[i]
        i += 1

    return res
