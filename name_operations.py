"""
This script contains a set of functions useful for conversions between different name formats and notations different functions require.
It is only imported in relationship_operations.py and sm_csm_operations.py. The user does not need to import it to the applicatory script.
It also contains the functions used for finding the intersections among lists of node names and node states.

Functions:

turn_to_name(a): Takes a dictionary of nodes and their states in plant pollinator (pl-po) notation, converts the dictionary to a stinrg,
and returns the string.

turn_string_to_list_n(rel_node,mapping): Takes a string (a node name in network of functional relationships) in plant pollinator (pl-po)
notation and returns the list of corresponding node states in number index notation.

def turn_string_to_dict_pl_po(string): Takes a string (a node name in network of functional relationships) in plant pollinator (pl-po)
notation and returns a dictionary of the nodes with their states in the same notation.

find_nodes_in_this_motif(primary, mapping): Takes a dictionary of nodes and their states in either notations and returns a list of node
states in number index notation.

intersection(list1, list2): Takes two lists of node states, and returns True if the two have mutual node states, and False otherwise.

intersection_v2(list1, list2): Takes two lists of node states, and returns True if the two have mutual nodes, and False otherwise.

intersection_items(list1, list2): Takes two lists, and returns the list of mutual items between the two lists.

intersection_negation(list1,list2): Takes two lists of node states, and returns True if one of the lists contains the negation of one node
state in the other list, and False otherwise.

Author: Fatemeh Sadat Fatemi Nasrollahi unless otherwise noted.
Date: May 2021
Python Version: 3.7
"""




def turn_to_name(a):


    """
    Converts a dictionary of nodes and their states in plant pollinator (pl-po) notation to a string.

    Keyword arguments:
        a -- a dictionary which is a motif/motif group such as {pl_1p:0, po_1p:0, pl_2p:0, po_2p:0}
    Returns:
        x -- A list containing a string. In the above example it returns ['pl_1p=0;po_1p=0;pl_2p=0;po_2p=0;']
    """


    a=a[0]
    x = []
    s = ''
    for item in a:
        s += '%s=%s;' % (item, str(a[item]))
    x.append(s)
    return x




def turn_string_to_list_n(rel_node,mapping):
    """
    Returns a list of node states in 'number index' notation.

    Keyword arguments:
        rel_node -- a string such as 'pl_0p=0;po_0p=0;pl_1p=0;po_1p=0;' in plant pollinator (pl-po) notation
        mapping -- a dictionary mapping from pl-po notation to 'number index' notation

    Returns:
        nodes_in_this_motif_node_number -- a list of node states in 'number index' notation  such as [~n1n, ~n2n, ~n3n, ~n4n]
    """


    nodes_in_this_motif_node_number = []

    node_states_in_the_rel_node = rel_node.split(';')[:-1]
    for i in range(len(node_states_in_the_rel_node)):
        node, state = node_states_in_the_rel_node[i].split('=')
        if state == '1':
            nodes_in_this_motif_node_number.append(mapping[node])
        else:
            nodes_in_this_motif_node_number.append(mapping['~' + node])

    return nodes_in_this_motif_node_number




def turn_string_to_dict_pl_po(string):
    """
    Receives a string and returns a dictionary of nodes and their states in plant pollinator (pl-po) notation.

    Keyword arguments:
        string -- a string such as 'pl_0p=0;po_0p=0;pl_1p=0;po_1p=0;' in pl-po notation

    Returns:
        nodes_in_this_motif -- a dictionary such as {pl_op:0, po_0p:0, pl_1p:0, po_1p:0}
    """


    nodes_in_this_motif = {}

    node_states_in_the_rel_node = string.split(';')[:-1]
    for i in range(len(node_states_in_the_rel_node)):
        node, state = node_states_in_the_rel_node[i].split('=')

        nodes_in_this_motif.update({node:int(state)})

    return nodes_in_this_motif




def find_nodes_in_this_motif(primary, mapping):
    """
    Receives a dictionary of nodes and their states in pl-po notation or 'number index' notation and returns
    a list of node states in 'number index' notation

    Keyword arguments:
        primary -- a dictionary of nodes and their states in the pl-po notation such as {pl_op:0, po_0p:0, pl_1p:0, po_1p:0} or
        in 'number index' notation such as {n1n:0, n2n:0, n3n:0, n4n:0}
        mapping -- a dictionary mapping from pl-po notation to 'number index' notation

    Returns:
        nodes_in_this_motif_node_number -- a list of node states in 'number index' notation  such as [~n1n, ~n2n, ~n3n, ~n4n]
    """


    nodes_in_this_motif_node_number = []
    for node in primary:
        if 'n' not in node:
            state = primary[node]
            if state == 1:
                nodes_in_this_motif_node_number.append(mapping[node])
            else:
                nodes_in_this_motif_node_number.append(mapping['~' + node])
        else:
            state = primary[node]
            if state == 1:
                nodes_in_this_motif_node_number.append(node)
            else:
                nodes_in_this_motif_node_number.append('~' + node)

    return nodes_in_this_motif_node_number




def intersection(list1, list2):
    """
    Receives two lists of node states and returns a boolean value showing if the two have mutual node states or not.

    Keyword arguments:
        list1 and list2 -- two lists containing node states

    Returns:
        Result -- a boolean value that is True if the two lists intersect, and False if they do not.
    """

    result = False


    for x in list1:
        for y in list2:

            # if one common
            if x == y:
                result = True
                return result

    return result




def intersection_v2(list1, list2):
    """
    Receives two lists of node states and returns a boolean value showing if the two have mutual nodes or not.

    Keyword arguments:
        list1 and list2 -- two lists containing node states

    Returns:
        Result -- a boolean value that is True if the two lists have a shared node, and False otherwise.
    """


    result = False
    new_list1 = [x.replace('~', '') for x in list1]
    new_list2 = [x.replace('~', '') for x in list2]
    # traverse in the 1st list
    for x in new_list1:

        # traverse in the 2nd list

        for y in new_list2:
            if x == y:
                result = True
                return result

    return result




def intersection_items(list1, list2):
    """
    Receives two lists and returns the items that are mutual in the two.

    Keyword arguments:
        list1 and list2 -- two lists containing node states

    Returns:
        result -- the list of mutual items that are in both lists
    """


    result=[]


    for x in list1:
        for y in list2:

            # if one common
            if x == y:
                if x not in result:
                    result.append(x)

    return result




def intersection_negation(list1,list2):
    """
    Receives two lists and returns a boolean value showing if the negation of a node state in one list is in the other list.

    Keyword arguments:
        list1 and list2 -- two lists containing node states

    Returns:
        True if the negation of a node in one list is in the other and False otherwise.
    """


    for item in list1:
        if '~' not in item:
            if '~' + item in list2:
                return True
        else:
            if item.replace('~', '') in list2:
                return True
    return False