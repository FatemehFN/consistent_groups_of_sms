

import itertools
import relationship_operations as RO
import BooleanDOI_processing as BDOIp
import FVS




def is_new_attractor(SCcom,attrs,G_rel):

    #returns True if the consistent combination of motifs/motif groups excludes all previously found attractor, and False otherwise.
    #inputs:
    #SCcom: the combination of motifs/motif groups
    #attrs: the list of attractors
    #G_rel: network of functional relationships between the stable motifs and motif groups in the form of a digraph

    bool_check = []
    for A in attrs:
        bool_check.append(RO.check_for_mutual_exclusivity(SCcom, A, G_rel))
    if all(i for i in bool_check):
        return True
    return False





def consistent_groups_of_sms(G_rel, list_of_names,lines):

    #calculates the number of q-attrs/minimal trap spaces based on the relationships between the motifs/motif groups
    #inputs:
    # G_rel: network of functional relationships between the motifs/motif groups. In this network nodes are motifs and motif groups and edges are the functional
    # relationships between the motifs/motif groups
    # FVS_size: size of the feedback vertex set
    # list_of_names: The list of the names of all stable motifs and conditionally stable motifs. This is for counting purposes.

    #output: number of consistent groups of motifs/motif groups, which is the same as number of q-attrs/minimal trap spaces


    #Extras:
    #alpha: the number of unstabilized nodes after substituting the specific motif/motif group
    #G_rel nodes format: [[{pl_1p:0, po_1p:0},{'number of sms':1, 'alpha': 14}],[{pl_2p:0, po_2p:0},{'number of sms':1, 'alpha': 0}],...]
    #G_rel edges format: [('pl_1p=0;po_1p=0;', 'pl_2p=0;po_2p=0;', {'relationship': 'DOI'}),...]


    G_rel = RO.merge_mutual_DOIs(G_rel)
    print('motifs/ motif groups in network of functional relationships')
    print(G_rel.nodes(data=True))



    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    FVS_size = len(FVS.FVS(Gread))


    nodes = G_rel.nodes(data=True)[::-1]
    attrs = []
    alpha_zero_nodes = [node for node in nodes if node[1]['alpha'] == 0]
    alpha_non_zero_nodes = [node1 for node1 in nodes if node1[1]['alpha'] != 0]


    for item in alpha_zero_nodes:
        if item not in attrs:
            if len(attrs) == 0:
                attrs.append([item])
            else:
                bool = False
                for A in attrs:
                    if G_rel[A[0][0]][item[0]]['relationship'] != 'mutual exclusivity':
                        bool = True
                if bool == False:
                    attrs.append([item])


    #build consistent groups of non-zero alpha motifs/motif groups
    if len(alpha_non_zero_nodes) > FVS_size:
        r = FVS_size
    else:
        r = len(alpha_non_zero_nodes)

    combs = []
    for i in range(2, r + 1):


        combo = list(itertools.combinations(alpha_non_zero_nodes, i))

        for item in combo:
            string = ''
            for n in item:
                string += n[0]
            if RO.how_many_motifs(string, list_of_names) <= r:
                combs.append(item)


    self_consistent_combs = []
    for comb in combs:
        if RO.check_for_mutual_exclusivity_in_a_comb([x[0] for x in comb],G_rel)==False:
            self_consistent_combs.append(comb)


    self_consistent_combs.sort(key=len)
    self_consistent_combs = self_consistent_combs[::-1]

    for SCcom in self_consistent_combs:
        if is_new_attractor(SCcom,attrs,G_rel)==True:
            attrs.append(SCcom)



    # To add the alpha non zeros that cannot be combined with anything else (oscillations)
    for node in alpha_non_zero_nodes:
        check_var = True
        for comparison_item1 in self_consistent_combs:
            if node in comparison_item1:
                check_var = False

        if check_var == True:
            all_bools = []
            for comparison_item in nodes:
                bool = True
                if node != comparison_item:
                    if G_rel.has_edge(node[0], comparison_item[0]):
                        if G_rel[node[0]][comparison_item[0]]['relationship'] != 'mutual exclusivity':
                            bool = False
                    else:
                        bool = False
                all_bools.append(bool)


            if all(z for z in all_bools):
                attrs.append([node])

    return attrs,len(attrs)