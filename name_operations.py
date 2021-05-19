




def turn_to_name(a):



    #converts a dictionary of nodes and their states in pl-po notation to a string
    #input: a dictionary which is a motif/motif group such as {pl_1p:0, po_1p:0, pl_2p:0, po_2p:0}
    #output: A list containing a string. In the above example it returns 'pl_1p=0;po_1p=0;pl_2p=0;po_2p=0;'



    a=a[0]
    x = []
    s = ''
    for item in a:
        s += '%s=%s;' % (item, str(a[item]))
    x.append(s)
    return x




def number_of_negative_edges_text(lines):



    #returns the number of negative edges in a Boolean model
    #input:
    #lines: the text file of Boolean model read with readlines() function
    #output: the number of negative edges in the model



    total_count=0
    for i in range (len(lines)):
        count = lines[i].count('not')
        total_count=total_count+count
    return total_count






def turn_string_to_list_n(rel_node,mapping):

    #returns a list of node states in 'number index' notation

    #inputs:
    #rel_node: a string such as 'pl_0p=0;po_0p=0;pl_1p=0;po_1p=0;' in pl-po notation
    #mapping: a dictionary mapping from pl-po notation to 'number index' notation

    #output: a list of node states in 'number index' notation  such as [~n1n, ~n2n, ~n3n, ~n4n]

    nodes_in_this_motif = []
    nodes_in_this_motif_node_number = []

    node_states_in_the_rel_node = rel_node.split(';')[:-1]
    for i in range(len(node_states_in_the_rel_node)):
        node, state = node_states_in_the_rel_node[i].split('=')

        if state == '1':
            nodes_in_this_motif_node_number.append(mapping[node])
            nodes_in_this_motif.append(node)
        else:
            nodes_in_this_motif_node_number.append(mapping['~' + node])
            nodes_in_this_motif.append('~' + node)

    return nodes_in_this_motif_node_number




def turn_string_to_dict_pl_po(string):

    #receives a string and returns a dictionary of nodes and their states in the pl-po notation

    #inputs:
    #string: a string such as 'pl_0p=0;po_0p=0;pl_1p=0;po_1p=0;' in pl-po notation
    #output: a dictionary such as {pl_op:0, po_0p:0, pl_1p:0, po_1p:0}


    nodes_in_this_motif = {}


    node_states_in_the_rel_node = string.split(';')[:-1]
    for i in range(len(node_states_in_the_rel_node)):
        node, state = node_states_in_the_rel_node[i].split('=')

        nodes_in_this_motif.update({node:int(state)})

    return nodes_in_this_motif






def find_nodes_in_this_motif(primary, mapping):

    # returns a list of node states in 'number index' notation

    # inputs:
    # primary: a dictionary of nodes and their states in the pl-po notation such as {pl_op:0, po_0p:0, pl_1p:0, po_1p:0} or
    # in 'number index' notation such as {n1n:0, n2n:0, n3n:0, n4n:0}
    # mapping: a dictionary mapping from pl-po notation to 'number index' notation

    # output: a list of node states in 'number index' notation  such as [~n1n, ~n2n, ~n3n, ~n4n]

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

    #returns a boolean value showing if two lists have mutual items or not.
    #input: two lists list1 and list2
    #output: True if the two lists intersect, and False if they do not.



    result = False


    for x in list1:
        for y in list2:

            # if one common
            if x == y:
                result = True
                return result

    return result





def intersection_v2(list1, list2):

    #returns a boolean value showing if two lists have mutual nodes or not.
    #input: two lists list1 and list2
    #output: True if the two lists intersect in a node (not node state), and False if they do not.

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

    # returns the items that are mutual in list1 and list2
    # input: two lists list1 and list2
    # output: the mutual items that are in both lists


    result=[]


    for x in list1:
        for y in list2:

            # if one common
            if x == y:
                if x not in result:
                    result.append(x)


    return result




def intersection_negation(list1,list2):

    #returns a boolean value showing if the negation of a node state in one list is in another list
    #inputs: two lists list1 and list2
    #output: True if the negation of a node in one list is in the other and False if not.

    for item in list1:
        if '~' not in item:
            if '~' + item in list2:
                return True
        else:
            if item.replace('~', '') in list2:
                return True
    return False











