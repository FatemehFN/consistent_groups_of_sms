
import PyBoolNet
import name_operations as NO
import relationship_operations as RO
import networkx as nx
import itertools
import BooleanDOI_processing as BDOIp
import model_operations as MO


def map(cycles):


    #generates the cycles_mapping that is a dictionary assigning a 'c'+number to each consistent cycle of the expanded network
    #example: {'c0': [{'n2': 0, 'n3': 0}, {'~n2', '~n0'}], 'c1': [{'n2': 0, 'n1': 0, 'n0': 0, 'n3': 0}, {'~n2', '~n0'}]}
    #input:
    #consistent cycles of the expanded network as the function consistent_cycles() of relationship_operations produces.
    #output: a map from c1,c2,... to the consistent cycles

    cycles_mapping = {}
    for i in range(len(cycles)):
        name = 'c' + str(i)
        cycles_mapping[name] = cycles[i]

    return cycles_mapping





def stable_motifs(r):


    #returns the stable motifs of a boolean model saved in r. r should be compatible with PyBoolNet format

    #input: Boolean functions compatible with PyBoolNet format
    #output: stable motifs of the Boolean model. Note that in plant pollinator networks we don't have self loops, hence we are
    # only looking for maximal trap spaces or stable motifs with length larger than 1. For other models, this line should be commented.

    primes = PyBoolNet.FileExchange.bnet2primes(r)
    PyBoolNet.PrimeImplicants.percolate_and_remove_constants(primes)
    maxts=PyBoolNet.AspSolver.trap_spaces(primes, "max")
    maxts=[x for x in maxts if len(x) > 1] #The line explained above
    return maxts







def order(csms,mapping):


    #pre-checks for the possibility of one csm containing the condition of another and sorts them accordingly
    #if csm1 contains any of the conditions of csm2, it is possible that csm1 is part of the support of csm2, so the support of csm1 should be found first.

    #inputs:
    #csms: list of conditionally stable motifs
    #mapping: a dictionary mapping from pl-po notation to 'number index' notation

    #output: the ordered list of conditionally stable motifs

    import math
    ignore = []
    swap = []
    done = []
    while True:
        bool = False
        bool1 = False
        Bool2 = False

        for i1 in range(len(csms)):

            for j1 in range(len(csms)):

                if i1 != j1 and [csms[i1], csms[j1]] not in swap and [csms[i1], csms[j1]] not in ignore and [csms[i1],csms[j1]] not in done:
                    if NO.intersection(list(csms[j1][1]), NO.find_nodes_in_this_motif(csms[i1][0], mapping)) and i1 > j1 and \
                            NO.intersection(NO.find_nodes_in_this_motif(csms[j1][0], mapping),NO.find_nodes_in_this_motif(csms[i1][0], mapping)) == False and \
                            NO.intersection(list(csms[i1][1]),NO.find_nodes_in_this_motif(csms[j1][0], mapping)) and j1 > i1 and \
                            NO.intersection(NO.find_nodes_in_this_motif(csms[i1][0], mapping),NO.find_nodes_in_this_motif(csms[j1][0], mapping)) == False:
                        ignore.append([csms[i1], csms[j1]])
                        ignore.append([csms[j1], csms[i1]])
                        bool = True
                        break



                    elif NO.intersection(list(csms[i1][1]), NO.find_nodes_in_this_motif(csms[j1][0], mapping)) and j1 > i1 and \
                            NO.intersection(NO.find_nodes_in_this_motif(csms[i1][0], mapping),
                                         NO.find_nodes_in_this_motif(csms[j1][0], mapping)) == False:
                        bool1 = True
                        bool = True

                        break

                    elif NO.intersection(list(csms[j1][1]), NO.find_nodes_in_this_motif(csms[i1][0], mapping)) and i1 > j1 and \
                            NO.intersection(NO.find_nodes_in_this_motif(csms[j1][0], mapping),
                                         NO.find_nodes_in_this_motif(csms[i1][0], mapping)) == False:
                        bool2 = True
                        bool = True
                        break

                    if bool == False:
                        done.append([csms[i1], csms[j1]])
                        done.append([csms[j1], csms[i1]])

            if bool == True:
                break

        if (len(swap) + len(ignore) + len(done)) / 2 == math.factorial(len(csms)) // math.factorial(
                (len(csms) - 2)) // math.factorial(2):
            break



        else:
            #print('here')
            #print(csms[i1])
            #print(csms[j1])
            #print(list(csms[i1][1]))
            #print(find_nodes_in_this_motif(csms[j1][0], mapping))

            if bool1 == True:
                csms[i1], csms[j1] = csms[j1], csms[i1]
                swap.append([csms[j1], csms[i1]])
                swap.append([csms[i1], csms[j1]])
            elif bool2 == True:
                csms[j1], csms[i1] = csms[i1], csms[j1]
                swap.append([csms[j1], csms[i1]])
                swap.append([csms[i1], csms[j1]])



    return csms


def self_consitence_check(comp,cycles_mapping,mapping):

    #returns True if an SCC in the cycle graph is consistent, and False otherwise
    #inputs:
    #comp: an SCC of cycle graph like {'c0', 'c1', 'c2'}
    #cycles_mapping: a dictionary assigning a 'c'+number to each consistent cycle of the expanded network. the function map() generates this
    #mapping: a dictionary mapping from pl-po notation to 'number index' notation


    list = []
    for c in comp:  # c:'c15'
        list += NO.find_nodes_in_this_motif(cycles_mapping[c][0], mapping)
    for l in list:
        if '~' + l in list:
            # print(l)
            return False
    return True









def csm_finder_positive_edges(cycles, mapping,write_cycle_graph):



    # the function that finds the conditionally stable motifs of a negative edge free network based on cycle graph constructed only using the inactive cycles
    # inputs:
    # cycles: consistent cycles of the expanded network as the function consistent_cycles() of relationship_operations produces.
    # mapping: a dictionary mapping from pl-po notation to 'number index' notation
    #write_cycle_graph: if it's True, cycle graph will be written in a gml file

    # output:
    # the ordered list of inactive csms





    cycle_graph = nx.DiGraph()

    cycles_mapping = {}
    for i in range(len(cycles)):
        name = 'c' + str(i)
        cycles_mapping[name] = cycles[i]

    key_list = list(cycles_mapping.keys())  # c1 c2
    val_list = list(cycles_mapping.values())  # list of csms


    # add the edges between conditional cycles
    conditional_cycles = [x for x in cycles if len(x[1]) != 0]
    checked=[]
    for i in range(len(conditional_cycles)):

        for j in range(len(conditional_cycles)):

            if i != j and [i,j] not in checked and [j,i] not in checked:

                checked.append([i,j])
                checked.append([j, i])
                n = list(conditional_cycles[i][0].values())
                m = list(conditional_cycles[j][0].values())
                if 1 not in n and 1 not in m:  # all inactive
                    position_i = val_list.index(conditional_cycles[i])
                    position_j = val_list.index(conditional_cycles[j])

                    if NO.intersection_negation(list(conditional_cycles[i][1]),list(conditional_cycles[j][1]))==False\
                        and NO.intersection_negation(list(conditional_cycles[i][1]),NO.find_nodes_in_this_motif(conditional_cycles[j][0], mapping))==False \
                        and NO.intersection_negation(list(conditional_cycles[j][1]),NO.find_nodes_in_this_motif(conditional_cycles[i][0], mapping))==False:




                        if NO.intersection(list(conditional_cycles[i][1]),NO.find_nodes_in_this_motif(conditional_cycles[j][0], mapping)) == True:


                            if cycle_graph.has_edge(key_list[position_j], key_list[position_i]) == False:
                                cycle_graph.add_edge(key_list[position_j], key_list[position_i])

                        if NO.intersection(list(conditional_cycles[j][1]), NO.find_nodes_in_this_motif(conditional_cycles[i][0],mapping)) == True:


                            if cycle_graph.has_edge(key_list[position_i], key_list[position_j]) == False:
                                cycle_graph.add_edge(key_list[position_i], key_list[position_j])





    if write_cycle_graph==True:
        nx.write_gml(cycle_graph, 'cycle_graph' + '.gml')





    cycle_graph_scc = [x for x in list(nx.strongly_connected_components(cycle_graph)) if len(x) > 1]




    csms = []
    SCCs_of_cycle_graph_including_sms = []







    # this loop gives the maximal sccs that are CSMs
    for item in cycle_graph_scc:
        if self_consitence_check(item,cycles_mapping,mapping) == True:
            buffer = {}
            buffer_cond = set()
            remaining_conds = []
            for c in item:

                buffer.update(cycles_mapping[c][0])
                buffer_cond.update(cycles_mapping[c][1])
            for cond in buffer_cond:
                if cond.replace('~', '') not in buffer:
                    remaining_conds.append(cond)

            SCCs_of_cycle_graph_including_sms.append([buffer, set(remaining_conds)])
            if len(remaining_conds) != 0:
                csms.append([buffer, set(remaining_conds)])



    # finding the cycles that did not participate in cycle graph SCC
    single_cycles_not_in_any_SCC = []

    for item in conditional_cycles:



            position_item = val_list.index(item)

            if item not in single_cycles_not_in_any_SCC:

                bool = False  # default: item is not in any of the SCCs
                for scc in cycle_graph_scc:
                    if key_list[position_item] in scc:
                        bool = True
                        break
                for scc_nodes in SCCs_of_cycle_graph_including_sms:

                    if NO.intersection(NO.find_nodes_in_this_motif(item[0], mapping),
                                    NO.find_nodes_in_this_motif(scc_nodes[0],mapping)):

                        bool = True
                        break
                if bool == False:
                    single_cycles_not_in_any_SCC.append(item)




    # fixing the conditions for this set
    for item in single_cycles_not_in_any_SCC:
        rc = []
        for c in item[1]:
            if c.replace('~', '') not in item[0].keys():
                rc.append(c)
        item[1] = set(rc)


    if csms!=[] and len(csms)>=2:
        csms = order(csms, mapping)
    csms = csms + single_cycles_not_in_any_SCC

    return csms






def csm_finder_general_function(cycles, mapping, G_expanded,write_cycle_graph):

    # the function that finds the conditionally stable motifs of any network based on cycle graph constructed using all consistent cycles
    # of the expanded network.
    # This function looks for smaller SCSs within the maximal SCCs in the cycle graph
    # This function constructs the cycle_graph_ME, a digragh object that has the mutual exclusivity information between the consistent cycles. It uses
    # the information in this object to break the maximal inconsistent SCC into consistent pieces.
    #function break_the_SCC() is embedded within this function to break the inconsistent SCCs of cycle graph into maximal consistent SCCs.


    # inputs:
    # cycles: consistent cycles of the expanded network as the function consistent_cycles() of relationship_operations produces.
    # mapping: a dictionary mapping from pl-po notation to 'number index' notation
    #write_cycle_graph: if it's True, cycle graph and cycle_graph_ME will be written in two separate gml files

    # output:
    # the ordered list of csms


    fishy=[]

    cycle_graph = nx.DiGraph()
    cycle_graph_ME=nx.DiGraph()

    cycles_mapping = {}
    for i in range(len(cycles)):
        name = 'c' + str(i)
        cycles_mapping[name] = cycles[i]

    key_list = list(cycles_mapping.keys())  # c1 c2
    val_list = list(cycles_mapping.values())  # list of csms


    # add the edges between conditional cycles
    conditional_cycles = [x for x in cycles if len(x[1]) != 0]
    checked=[]
    for i in range(len(conditional_cycles)):

        for j in range(len(conditional_cycles)):

            if i != j and [i,j] not in checked and [j,i] not in checked:

                    checked.append([i,j])
                    checked.append([j, i])
                    position_i = val_list.index(conditional_cycles[i])
                    position_j = val_list.index(conditional_cycles[j])

                    if NO.intersection_negation(list(conditional_cycles[i][1]),list(conditional_cycles[j][1]))==False\
                        and NO.intersection_negation(list(conditional_cycles[i][1]),NO.find_nodes_in_this_motif(conditional_cycles[j][0], mapping))==False \
                        and NO.intersection_negation(list(conditional_cycles[j][1]),NO.find_nodes_in_this_motif(conditional_cycles[i][0], mapping))==False:




                        if NO.intersection(list(conditional_cycles[i][1]),NO.find_nodes_in_this_motif(conditional_cycles[j][0], mapping)) == True:


                            if cycle_graph.has_edge(key_list[position_j], key_list[position_i]) == False:
                                cycle_graph.add_edge(key_list[position_j], key_list[position_i])

                        if NO.intersection(list(conditional_cycles[j][1]), NO.find_nodes_in_this_motif(conditional_cycles[i][0],mapping)) == True:


                            if cycle_graph.has_edge(key_list[position_i], key_list[position_j]) == False:
                                cycle_graph.add_edge(key_list[position_i], key_list[position_j])

                    else:
                        if cycle_graph_ME.has_edge(key_list[position_j], key_list[position_i]) == False:
                                cycle_graph_ME.add_edge(key_list[position_j], key_list[position_i])
                                cycle_graph_ME.add_edge(key_list[position_i], key_list[position_j])



    if write_cycle_graph==True:
        nx.write_gml(cycle_graph, 'cycle_graph' + '.gml')
        nx.write_gml(cycle_graph_ME, 'cycle_graph_ME' + '.gml')







    cycle_graph_scc = [x for x in list(nx.strongly_connected_components(cycle_graph)) if len(x) > 1]




    csms = []
    SCCs_of_cycle_graph_including_sms = []





    def break_the_SCC(item): #item: {'c15', 'c14'}
        csms_b=[]
        checked=[]
        for c in item:
            edges=[E for E in (cycle_graph_ME.edges()) if E[1]==c]
            forbidden=[x[0] for x in edges]
            if forbidden!=[]:
                allowed=[c1 for c1 in list(item) if c1 not in forbidden]
                non_ME_SCCs=list(nx.strongly_connected_components(cycle_graph.subgraph(allowed)))
                for a in non_ME_SCCs:
                    if a not in checked:

                        buffer = {}
                        buffer_cond = set()
                        remaining_conds = []
                        for c2 in a:
                            buffer.update(cycles_mapping[c2][0])
                            buffer_cond.update(cycles_mapping[c2][1])
                        for cond in buffer_cond:
                            if cond.replace('~', '') not in buffer:
                                remaining_conds.append(cond)
                        if [buffer, set(remaining_conds)] not in SCCs_of_cycle_graph_including_sms:
                            SCCs_of_cycle_graph_including_sms.append([buffer, set(remaining_conds)])
                        if len(remaining_conds) != 0 and [buffer, set(remaining_conds)] not in csms_b:
                            csms_b.append([buffer, set(remaining_conds)])
                        checked.append(a)


        return csms_b



    # this loop gives the maximal sccs that are CSMs
    for item in cycle_graph_scc:
        if self_consitence_check(item,cycles_mapping,mapping) == True:
            buffer = {}
            buffer_cond = set()
            remaining_conds = []
            for c in item:

                buffer.update(cycles_mapping[c][0])
                buffer_cond.update(cycles_mapping[c][1])
            for cond in buffer_cond:
                if cond.replace('~', '') not in buffer:
                    remaining_conds.append(cond)

            SCCs_of_cycle_graph_including_sms.append([buffer, set(remaining_conds)])
            if len(remaining_conds) != 0:
                csms.append([buffer, set(remaining_conds)])

        else:

            broken_csms=break_the_SCC(item)
            if broken_csms!=[]:
                csms+=break_the_SCC(item)




    # finding the cycles that did not participate in cycle graph SCC
    single_cycles_not_in_any_SCC = []

    for item in conditional_cycles:



            position_item = val_list.index(item)

            if item not in single_cycles_not_in_any_SCC:

                bool = False  # default: item is not in any of the SCCs
                for scc in cycle_graph_scc:
                    if key_list[position_item] in scc:
                        bool = True
                        break
                for scc_nodes in SCCs_of_cycle_graph_including_sms:

                    if NO.intersection(NO.find_nodes_in_this_motif(item[0], mapping),
                                    NO.find_nodes_in_this_motif(scc_nodes[0],mapping)):

                        bool = True
                        break
                if bool == False:
                    single_cycles_not_in_any_SCC.append(item)




    # fixing the conditions for this set
    for item in single_cycles_not_in_any_SCC:
        rc = []
        for c in item[1]:
            if c.replace('~', '') not in item[0].keys():
                rc.append(c)
        item[1] = set(rc)


    # finding state=1 motifs
    one_motifs = []
    one_motifs_DOI = []

    # finding the unconditional cycles (SMs)
    unconditional_cycles = [x for x in cycles if len(x[1]) == 0]

    for item in unconditional_cycles:


        if 0 not in list(item[0].values()):
            one_motifs.append(NO.find_nodes_in_this_motif(item[0], mapping))
            one_motifs_DOI.append(RO.find_DOIs(G_expanded, NO.find_nodes_in_this_motif(item[0], mapping)))


    # finding the bridge
    bridge_list = []
    for doi in one_motifs_DOI:
        for scc in SCCs_of_cycle_graph_including_sms:
            if NO.intersection(NO.find_nodes_in_this_motif(scc[0], mapping), doi):
                bridge_list.append(NO.intersection_items(NO.find_nodes_in_this_motif(scc[0], mapping), doi))


    for b in range(len(bridge_list)):
        bridge_list[b] = list(set(bridge_list[b]))

    bridge = list(set(itertools.chain.from_iterable(bridge_list)))

    # the all=1 stable motifs that have the bridge in their DOI
    one_motif_bridge = []
    one_motif_bridge_DOI = []
    for i in range(len(one_motifs_DOI)):
        if NO.intersection(one_motifs[i] + one_motifs_DOI[i], bridge):
            one_motif_bridge.append(one_motifs[i])
            one_motif_bridge_DOI.append(one_motifs_DOI[i])


    bridge_list_inactive = []
    for br in bridge_list:
        for item in br:
            if '~' in item:
                bridge_list_inactive.append(br)
                break


    if len(bridge) != 0:


        for item in bridge_list_inactive:
            forbidden_cycles = []
            for c in cycles_mapping:
                if NO.intersection(NO.find_nodes_in_this_motif(cycles_mapping[c][0], mapping),
                                item) and c not in forbidden_cycles:
                    forbidden_cycles.append(c)


            # find the intermediate SCCs
            intermediate = []
            for item in cycle_graph_scc:
                combs = []
                allowed_cycles = []
                for c in item:
                    if c not in forbidden_cycles:
                        allowed_cycles.append(c)


                if len(item) == len(
                        allowed_cycles):  # means all the cycles inside this scc are allowed and this scc does not have the bridge in it
                    continue

                intermediate += list(nx.strongly_connected_components(cycle_graph.subgraph(allowed_cycles)))

            all_intermediate_SCSs = []

            for item in intermediate:
                buffer = {}
                buffer_cond = set()
                remaining_conds = []
                for c in item:

                    buffer.update(cycles_mapping[c][0])
                    buffer_cond.update(cycles_mapping[c][1])
                for cond in buffer_cond:
                    if cond.replace('~', '') not in buffer:
                        remaining_conds.append(cond)

                all_intermediate_SCSs.append([buffer, set(remaining_conds)])



            for SCS in all_intermediate_SCSs:

                if len(SCS[1]) != 0 and SCS not in csms and SCS not in fishy:

                    fishy.append(SCS)




    csms = csms + fishy
    if csms!=[] and len(csms)>=2:
        csms = order(csms, mapping)
    csms = csms + single_cycles_not_in_any_SCC
    #print('csms')
    #print(csms)
    return csms























def find_supports(csm, mapping, G_expanded, G_rel,FVS_size,list_of_names, number_of_neg_edges):


    # returns the list of supports of each CSM
    # inputs:
    # csm: conditionally stable motif
    # mapping: a dictionary mapping from pl-po notation to 'number index' notation
    # G_expanded: the expanded network as Gang Yang's function Get_expanded_network() produces
    # G_rel: network of functional relationships between the stable motifs and motif groups in the form of a digraph
    # FVS_size: feedback vertex set size
    # list_of_names: a list of strings, each of which is a stable motif or a conditionally stable motif such as ['pl_0p=0;po_0p=0;', 'pl_1p=0;po_1p=0;']
    # number_of_neg_edges: number of negative edges in the network

    G_rel_nodes=[]
    if number_of_neg_edges==0:
        for g in G_rel.nodes(data=True):
            if '0' in g[0]:
                G_rel_nodes.append(g)


    else:
        G_rel_nodes=G_rel.nodes(data=True)




    l = [item for item in range(0, len(G_rel_nodes))]
    support_dict=[]
    combo_index=[]
    for i in range (1,FVS_size):







        combo_index+=list(itertools.combinations(l, i))

        for item in combo_index:

            string = ''
            for n in item:
                string += G_rel_nodes[n][0]




            if RO.check_for_mutual_exclusivity_in_a_comb([G_rel_nodes[x][0] for x in item],G_rel)==False:

                if RO.repetition([G_rel_nodes[x] for x in item],mapping)==False and RO.DOI_check([G_rel_nodes[x][0] for x in item],G_rel)==False:



                    if RO.how_many_motifs(string, list_of_names) == i:




                        #find the node combo complete list and its DOI
                        nodes_n_all=NO.turn_string_to_list_n(string,mapping)
                        DOI=RO.find_DOIs(G_expanded,nodes_n_all)



                        #check if the DOI U node combo includes the conditions of the csm
                        if csm[1].issubset(set(DOI + nodes_n_all)): #and \




                            if NO.intersection_v2(list(csm[0]), nodes_n_all) == False and NO.intersection_negation(list(csm[1]), RO.find_DOIs_residue(G_expanded, nodes_n_all)) == False:


                                if RO.is_subset(NO.turn_string_to_dict_pl_po(string),support_dict,mapping)==False:

                                    support_dict.append([NO.turn_string_to_dict_pl_po(string), i + 1])





    return support_dict










def cycle_graph_virtual_node_based(lines,write_cycle_graph=False):

    #constructs cycle graph using the virtual nodes. For each virtual node in the expanded network adds an edge from
    # the cycles that contain that virtual node inside them to the cycles that have that virtual node in their composite nodes and
    # are consistent


    #inputs:
    #lines: Boolean functions read from the text file using function readlines()
    #write_cycle_graph: if it's True, cycle graph and cycle_graph_ME will be written in two separate gml files

    #output: writes cycle graph in a gml file.




    # build the network from the lines
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)


    # build the expanded network
    prefix, suffix = 'n', ''
    G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)





    # calculate the mapping from string nodename to index
    mapping = {}
    inverse_mapping = {}
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node





    # add composite node to node mapping
    for node in G_expanded.nodes():
        if node not in mapping:
            components = node.split('_')
            composite_node = '_'.join([inverse_mapping[x] for x in components])
            mapping[composite_node] = node
            inverse_mapping[node] = composite_node





    cycles = RO.consistent_cycles(G_expanded)
    conditional_cycles=[x for x in cycles if len(x[1]) != 0]







    cycle_graph = nx.DiGraph()




    cycles_mapping = {}
    for i in range(len(cycles)):
        name = 'c' + str(i)
        cycles_mapping[name] = cycles[i]

    print('cycle mapping')
    for item in cycles_mapping:
        print(item,cycles_mapping[item])


    key_list = list(cycles_mapping.keys())  # c1 c2
    val_list = list(cycles_mapping.values())  # list of csms







    virtual_nodes=[x for x in G_expanded.nodes() if '_' not in x]


    for i in virtual_nodes:


        source_cycles=[y for y in conditional_cycles if i in NO.find_nodes_in_this_motif(y[0],mapping)]
        target_cycles=[z for z in conditional_cycles if i in list(z[1])]
        for sc in source_cycles:
            for tc in target_cycles:
                if tc!=sc:
                    if NO.intersection_negation(list(sc[1]), list(tc[1])) == False \
                            and NO.intersection_negation(list(sc[1]), NO.find_nodes_in_this_motif(tc[0],mapping)) == False \
                            and NO.intersection_negation(list(tc[1]),NO.find_nodes_in_this_motif(sc[0],mapping)) == False:



                        position_sc = val_list.index(sc)
                        position_tc = val_list.index(tc)


                        if cycle_graph.has_edge(key_list[position_sc], key_list[position_tc]) == False:



                            cycle_graph.add_edge(key_list[position_sc], key_list[position_tc])


    if write_cycle_graph==True:
        nx.write_gml(cycle_graph, 'cycle_graph.gml')




    return cycle_graph









