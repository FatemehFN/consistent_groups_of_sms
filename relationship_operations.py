"""
This script contains a set of functions useful for constructing the network of functional relationships between the stable motifs
and motif groups. It also contains the functions that find the functional relationships such as mutual exclusivity and LDOI among
the motifs and motif groups

Functions:

construct_relationships_network(lines,rules,write_cycle_graph): Takes the Boolean functions read from a text file and returns
network of functional relationships in a networkx DiGraph object.

consistent_cycles(G_expanded): Takes expanded network and returns the consistent cycles of it.

find_DOIs(G_expanded,nodes_in_this_motif_node_number): Takes expanded network and a list of node states and returns the logical
domain of influence of the list of node states.

find_DOIs_residue(G_expanded,nodes_in_this_motif_node_number): Takes the expanded network and a list of node states and returns
the conflicting points during the LDOI search.

check_for_mutual_exclusivity_in_a_comb(comb,G_rel): Takes a list of names of motifs and motif groups and returns True if the
list is not consistent, and False otherwise.

check_for_mutual_exclusivity(comb1, comb2, G_rel): Takes two lists of motifs and motif groups and returns True if the
the two are not consistent, and False otherwise.

DOI_check(comb,G_rel): Takes a list of motifs and motif groups and returns True if it contains two motifs/ motif groups that are in the
LDOI of each other, and False otherwise.

how_many_motifs(string,list_of_names): Takes s string and the list of names of stable motifs and conditionally stable motifs and counts
how many stable motifs and conditionally stable motifs are within the string.

is_subset(new_dict,support_dict,mapping): Takes a dictionary of a newly found combination of motifs and motif group that satisfy the
conditions of a conditionally stable motif, the dictionary of the found supports and the mapping. Returns True if the new combination
contains any of the found supports, and False otherwise.

repetition(list,mapping): Takes a list of names of stable motifs and motif groups and returns True if the list contains repition
and False otherwise

merge_mutual_DOIs(G_rel): Takes the network of functional relationships, and merges the nodes that have LDOI relationships. Returns the
network of functional relationships after the merging.

Author: Fatemeh Sadat Fatemi Nasrollahi.
Date: May 2021
Python Version: 3.7
"""


import networkx as nx
import BooleanDOI_DOI as BDOI
from collections import ChainMap
import cool_bool_tools as cbt
import sm_csm_operations as SMO
import model_operations as MO
import name_operations as NO
import BooleanDOI_processing as BDOIp
import FVS
import copy




def construct_relationships_network(lines,rules,write_cycle_graph):


    """
    Constructs the cycle graph and saves it in a digraph object in which nodes are stable motifs and edges are functional
    relationships between them.

    Keyword arguments:
        lines -- Boolean functions read from the text file using function readlines()
        rules -- Boolean functions read from the text file using function read()

    Returns:
        G_rel -- the network of functional relationships
    """


    G_rel=nx.DiGraph()
    len_rel_node = []
    primary_sm=[]
    DOIs = []
    nodes_in_all_motifs = []


    #build the network from the lines
    Gread,readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    FVS_size = len(FVS.FVS(Gread))



    #build the expanded network
    prefix, suffix = 'n', ''
    G_expanded=BDOIp.Get_expanded_network(Gread,prefix=prefix,suffix=suffix)

    #calculate the mapping from string nodename to index
    mapping={}              #nodename to number index
    inverse_mapping={}      #number index to nodename
    for i,node in enumerate(readnodes):
      index=prefix+str(i)+suffix
      mapping[node]=index
      inverse_mapping[index]=node
      mapping['~'+node]='~'+index
      inverse_mapping['~'+index]='~'+node


    #add composite node to node mapping
    for node in G_expanded.nodes():
        if node not in mapping:
            components=node.split('_')
            composite_node='_'.join([inverse_mapping[x] for x in components])
            mapping[composite_node]=node
            inverse_mapping[node]=composite_node



    #find stable motifs
    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')

    sm_list=copy.deepcopy(SMO.stable_motifs(rules))

    print('stable motifs')
    print(sm_list)

    #add the SM nodes to G_rel
    for item in sm_list:

        name=NO.turn_to_name([item])
        G_rel.add_nodes_from(name)
        len_rel_node.append({name[0]:1})
        primary_sm.append([item,{'number of sms':1}])
        nodes_in_all_motifs.append(NO.find_nodes_in_this_motif(item,mapping))
        DOIs.append(find_DOIs(G_expanded,NO.find_nodes_in_this_motif(item,mapping)))

    len_rel_node = dict(ChainMap(*len_rel_node))
    nx.set_node_attributes(G_rel, 'number of sms', len_rel_node)


    #add the G_rel edges
    done=[]
    for q in range (len(nodes_in_all_motifs)):
            for m in range (len(nodes_in_all_motifs)):
                if m!=q and NO.intersection_negation(nodes_in_all_motifs[q],nodes_in_all_motifs[m]):
                    if [q,m] not in done and [m,q] not in done:

                        G_rel.add_edge(G_rel.nodes()[m], G_rel.nodes()[q], relationship='mutual exclusivity')
                        G_rel.add_edge(G_rel.nodes()[q], G_rel.nodes()[m], relationship='mutual exclusivity')
                        done.append([m,q])
                        done.append([q,m])
                if m!=q and set(nodes_in_all_motifs[q]).issubset(set(nodes_in_all_motifs[m]+DOIs[m])):
                    G_rel.add_edge(G_rel.nodes()[m],G_rel.nodes()[q],  relationship='DOI')
                if m!=q and set(nodes_in_all_motifs[m]).issubset(set(nodes_in_all_motifs[q]+DOIs[q])):
                    G_rel.add_edge(G_rel.nodes()[q],G_rel.nodes()[m],  relationship='DOI')



    list_of_names = G_rel.nodes()


    #consistent cycles of expanded network
    formatted_M_list=consistent_cycles(G_expanded)


    #conditionally stable motifs
    number_of_neg_edges = MO.number_of_negative_edges_text(lines)

    if number_of_neg_edges==0:
        csms= SMO.csm_finder_positive_edges(formatted_M_list, mapping,write_cycle_graph=write_cycle_graph)
    else:
        csms= SMO.csm_finder_general_function(formatted_M_list, mapping,G_expanded,write_cycle_graph=write_cycle_graph)




    #complete the list of names
    for item in csms:
        buffer={}
        for c in item[0]:
            buffer.update({inverse_mapping[c]: item[0][c]})

        name= NO.turn_to_name([buffer])
        list_of_names+=name



    #find supports of conditionally stable motifs
    found_index = []
    found_counter = 0
    found_counter_array = []
    while (found_counter != len(csms)):

        if found_counter != 0 and len(found_counter_array) >= 2:
            if found_counter_array[-1] == found_counter_array[-2]:
                break

        for ind in range(len(csms)):

            if ind not in found_index:
                print('\n')
                print('csm')
                print(csms[ind])
                support_dict = SMO.find_supports(csms[ind], mapping, G_expanded, G_rel,FVS_size,list_of_names, number_of_neg_edges)

                print('supports')
                print(support_dict)
                if support_dict != []:
                    found_counter = found_counter + 1
                    found_index.append(ind)
                    buffer = {}

                    for c in csms[ind][0]:  # c is the key
                        buffer.update({inverse_mapping[c]: csms[ind][0][c]})

                    for pp in support_dict:

                        succ = [[{**pp[0], **buffer}, {'number of sms': pp[1]}]]
                        succ1 = {**pp[0], **buffer}

                        primary_sm = primary_sm + succ
                        name= NO.turn_to_name([succ1])
                        if name[0] not in G_rel.nodes():
                            G_rel.add_nodes_from(name)
                        len_rel_node.update({name[0]: pp[1]})
                        nx.set_node_attributes(G_rel, 'number of sms', len_rel_node)
                        n = NO.find_nodes_in_this_motif(succ1, mapping)
                        doi_n = find_DOIs(G_expanded, n)

                        for m in G_rel.nodes():
                            if m != name[0]:
                                node_number_m = NO.find_nodes_in_this_motif(NO.turn_string_to_dict_pl_po(m),mapping)
                                if NO.intersection_negation(n, node_number_m):
                                    G_rel.add_edge(m, name[0], relationship='mutual exclusivity')
                                    G_rel.add_edge(name[0], m, relationship='mutual exclusivity')

                                if set(n).issubset(set(node_number_m + find_DOIs(G_expanded, node_number_m))):
                                    G_rel.add_edge( m,name[0], relationship='DOI')
                                if set(node_number_m).issubset(set(n + doi_n)):
                                    G_rel.add_edge(name[0],m,  relationship='DOI')



        found_counter_array.append(found_counter)

    alpha=[] # alpha -- the number of unstabilized nodes after substituting the specific motif/motif group
    for rel_node in G_rel.nodes():


        dict_nodes=NO.turn_string_to_dict_pl_po(rel_node)
        alpha_single=MO.number_of_unstabilized_nodes(dict_nodes,rules,len(readnodes))
        alpha.append({rel_node:alpha_single})

    alpha=dict(ChainMap(*alpha))
    nx.set_node_attributes(G_rel, 'alpha', alpha)

    return G_rel,list_of_names




def consistent_cycles(G_expanded):


    """
    Returns a list of consistent cycles of the expanded network in 'number index' notation.

    Keyword arguments:
        G_expanded -- the expanded network as Gang Yang's function Get_expanded_network() produces

    Returns:
        formatted_M_list -- a list of consistent cycles of the expanded network in 'number index' notation
        example -- [[{'n0': 0, 'n1': 0}, {'~n2', '~n0'}], [{'n3': 0, 'n2': 0}, {'~n2', '~n0'}]]
        the first member of each entity is a dictionary containing the virtual nodes inside the cycle and the second
        member is the composite nodes broken down.
    """


    Ms = []
    M_list = []
    for i in list(nx.simple_cycles(G_expanded)):
        if cbt.are_subsets_consistent(i, i):
            Ms.append(i)

    for m in Ms:
        conditions = []
        for n in m:
            if '_' in n:
                for cc in n.split('_'):

                    if cc not in conditions:
                        conditions.append(cc)
        M_list.append([m, conditions])

    formatted_M_list = []
    for item in M_list:
        buffer = {}
        for i in item[0]:
            if '_' not in i:
                if '~' in i:
                    buffer[i.replace('~', '')] = 0
                else:
                    buffer[i] = 1
        if [buffer, set(item[1])] not in formatted_M_list:
            formatted_M_list.append([buffer, set(item[1])])



    return formatted_M_list




def find_DOIs(G_expanded,nodes_in_this_motif_node_number):


    """
    Returns the LDOI of a list of node states in 'number index' notation such as [~n0n, ~n1n].

    Keyword arguments:
        G_expanded -- the expanded network as Gang Yang's function Get_expanded_network() produces
        nodes_in_this_motif_node_number -- a list of node states in 'number index' notation such as [~n0n, ~n1n]

    Returns:
        list(DOI_of_a_m_node_numbers[0]) -- the LDOI of the list of node states given in nodes_in_this_motif_node_number
        in ' number index' notation. For example [~n2n, ~n3n, n4n, ~n3n_n4n]
    """


    DOI_of_a_m_node_numbers = (BDOI.truncated_node_of_influence_BFS(G_expanded, nodes_in_this_motif_node_number))
    return list(DOI_of_a_m_node_numbers[0])




def find_DOIs_residue(G_expanded,nodes_in_this_motif_node_number):


    """
    Returns the conflicting node states to the source during the search for LDOI of a list of node states.

    Keyword arguments:
        G_expanded -- the expanded network as Gang Yang's function Get_expanded_network() produces
        nodes_in_this_motif_node_number: a list of node states in 'number index' notation such as [~n0n, ~n1n]

    Returns:
        list(DOI_of_a_m_node_numbers[2]) -- a list of conflicting node states found during LDOI search. For example [n0n]
    """


    DOI_of_a_m_node_numbers = (BDOI.truncated_node_of_influence_BFS(G_expanded, nodes_in_this_motif_node_number))
    return list(DOI_of_a_m_node_numbers[2])




def check_for_mutual_exclusivity_in_a_comb(comb,G_rel):


    """
    Returns a boolean value showing if a combination of stable motifs and motif groups is self-consistent or not
    to find this, it uses the info stored in G_rel edges

    Keyword arguments:
        comb -- a list of stable motifs and/or motif groups. Just the names in strings without data=True
        G_rel -- network of functional relationships between the stable motifs and motif groups in the form of a digraph

    Returns:
        True if the list contains two motifs/motif groups that are not consistent, and False otherwise
    """


    for i in range (len(comb)):
        for j in range (len(comb)):
            if i!=j:
                if G_rel.has_edge(comb[i],comb[j]):
                    if G_rel [comb[i]][comb[j]]['relationship']=='mutual exclusivity':
                        return True
    return False




def check_for_mutual_exclusivity(comb1, comb2, G_rel):


    """
    Returns a boolean value showing if two combinations of stable motifs and motif groups are consistent or not
    to find this, it uses the info stored in G_rel edges.

    Keyword arguments:
        comb1, comb2 -- two lists of stable motifs and/or motif groups with data=True
        G_rel -- network of functional relationships between the stable motifs and motif groups in the form of a digraph

    Returns:
        True if the two lists contains two motifs/motif groups that are not consistent, and False otherwise
    """


    for i in range(len(comb1)):
        for j in range(len(comb2)):
            if G_rel.has_edge(comb1[i][0], comb2[j][0]):
                if G_rel[comb1[i][0]][comb2[j][0]]['relationship'] == 'mutual exclusivity':
                    return True
    return False




def DOI_check(comb,G_rel):


    """
    Returns a boolean value showing if a combination of stable motifs and motif groups has two items that are in LDOI of each other
    to find this, it uses the info stored in G_rel edges.

    Keyword arguments:
        comb -- a list of stable motifs and/or motif groups
        G_rel -- network of functional relationships between the stable motifs and motif groups in the form of a digraph

    Returns:
        True if the list contains two motifs/motif groups that are in LDOI of each other , and False otherwise
    """


    for i in range (len(comb)):
        for j in range (len(comb)):
            #print(i)
            #print(j)
            if i!=j:
                if G_rel.has_edge(comb[i],comb[j]):
                    if G_rel [comb[i]][comb[j]]['relationship']=='DOI':
                        if G_rel.has_edge(comb[j], comb[i]):
                            if G_rel [comb[j]][comb[i]]['relationship']=='DOI':
                                return True
    return False




def how_many_motifs(string,list_of_names):


    """
    Returns how many stable motifs and conditionally stable motifs are included within one string

    Keyword arguments:
        string -- a string of nodes with their states in pl-po notation such as 'pl_0p=0;po_0p=0;pl_1p=0;po_1p=0;'
        list_of_names -- a list of strings, each of which is a stable motif or a conditionally stable motif such as ['pl_0p=0;po_0p=0;', 'pl_1p=0;po_1p=0;']

    Returns:
        count -- the number showing how many motifs are included in a string. In the above example the function returns 2.
    """


    count=0

    until_here=[]
    for item in list_of_names:
        bool = []

        list=item.split(';')[:-1]
        for l in list:
            if l in string:

                bool.append(True)
            else:
                bool.append(False)

        if all([x for x in bool]):
            count+=1
            until_here+=list

        if set(string.split(';')[:-1]).issubset(set(until_here)) or set(string.split(';')[:-1])==set(until_here):
            break


    return count




def is_subset(new_dict,support_dict,mapping):


    """
    Returns a boolean value showing if any member of the list of previously found supports of a conditionally stable motif is
    a subset of a newly found one
    for example if the combination of motifs 1, 2 is a support, then the combination of motifs 1, 2, and 3 is not a support.

    Keyword arguments:
        new_dict -- the new combination that supports the conditionally stable motif in the form of a dictionary
        support_dict -- a list of previously found supports of a conditionally stable motif
        mapping -- a dictionary mapping from pl-po notation to 'number index' notation

    Returns:
        True if the newly found combination contains any of the previously found ones, and False otherwise
    """


    for s in support_dict:
        if set(NO.find_nodes_in_this_motif(s[0],mapping)).issubset(set(NO.find_nodes_in_this_motif(new_dict,mapping))):
            return True

    return False




def repetition(list,mapping):


    """
    Returns a boolean value showing if there is repition among members of a list of motifs/motif groups. repition is allowed between motif groups
    because two conditionally stable motifs could have the same supports, but it's not allowed between a motif group and a stable motif.

    Keyword arguments:
        list -- a list of motifs/motif groups such as ['pl_0p=0;po_0p=0;', 'pl_1p=0;po_1p=0;']
        mapping -- a dictionary mapping from pl-po notation to 'number index' notation
    Returns:
        True if the list contains two motifs/motif groups that are repetitive , and False otherwise
    """


    for item1 in list:
        for item2 in list:
            if item1!=item2:
                if NO.intersection(NO.turn_string_to_list_n(item1[0],mapping),NO.turn_string_to_list_n(item2[0],mapping)):
                    if item1[1]['number of sms']==1 and item2[1]['number of sms']!=1 or item1[1]['number of sms']==1 and item2[1]['number of sms']==1:
                        return True
    return False




def merge_mutual_DOIs(G_rel):


    """
    Receives network of functional relationships and merges the nodes that have LDOI relationships into one of them keeping the edges
    of those that are merged going to and coming from the kept one.

    Keyword arguments:
        G_rel -- network of functional relationships in the form of a networkx digraph object

    Returns:
        G_rel -- network of functional relationships G_rel with nodes that have LDOI relationships merged into one
    """


    def mutual_DOI_detection(G_rel):
        for node1 in G_rel.nodes(data=True):
            for node2 in G_rel.nodes(data=True):
                if node1[0]!=node2[0]:
                    if G_rel.has_edge(node1[0], node2[0]) and G_rel.has_edge(node2[0], node1[0]):
                        if G_rel[node1[0]][node2[0]]['relationship'] == 'DOI' and G_rel[node2[0]][node1[0]]['relationship'] == 'DOI':
                            return [True,[node1,node2]]
        return [False,[]]



    while True:
        [bool_mutual_DOI, items] = mutual_DOI_detection(G_rel)
        if bool_mutual_DOI==False:
            break
        else:
            G_rel=nx.contracted_nodes(G_rel,items[0][0],items[1][0])


    return G_rel