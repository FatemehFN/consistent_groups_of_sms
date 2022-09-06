import relationship_operations as RO
import name_operations as NO
import BooleanDOI_processing as BDOIp
import itertools
import networkx as nx
import PyBoolNet
import random
import structure_operations as SO
import numpy as np
import model_operations as MO

def redundancy(list,meta_list):

    for item in meta_list:
        if set(item).issubset(set(list)):
            return True

    return False





def find_driver(motif,max_size,G_expanded,mapping,driver_previous_motif=[]):


    string_name_motif=motif
    #print(string_name_motif)

    drivers = []
    node_list_n=NO.turn_string_to_list_n(string_name_motif,mapping=mapping)
    #print(node_list_n)
    node_dict_pl_po=NO.turn_string_to_dict_pl_po(string_name_motif)



    def part_of_the_motif(node, node_list_n):

        for n1 in node_list_n:
            for n2 in node_list_n:
                if n1!=n2 and G_expanded.has_edge(node,n1) and G_expanded.has_edge(n2,node):
                    return True

        return False







    if 0 not in node_dict_pl_po.values():
        #for n in node_list_n:
            #drivers.append([n])
        drivers=[node_list_n[0]]

    else:
        all_composite_nodes=[x for x in G_expanded.nodes() if '_' in x]
        relevant_composite_nodes=[]
        for node in all_composite_nodes:
            #print('here5')
            if part_of_the_motif(node,node_list_n)==True:
                    relevant_composite_nodes.append(node)

        #print(relevant_composite_nodes)
        if relevant_composite_nodes==[]:
            drivers = [node_list_n[0]]
            #print('drivers')
            #print(drivers)

        else:
            target_nodes_of_composite_nodes=[]
            for cn in relevant_composite_nodes:
                target_nodes_of_composite_nodes+=[x for x in node_list_n if G_expanded.has_edge(cn,x)]

            #print('targets of composite nodes')
            #rint(target_nodes_of_composite_nodes)
            combo=[]
            for i in range (1,len(target_nodes_of_composite_nodes)+1):
                #print('here4')
                combo += list(itertools.combinations(target_nodes_of_composite_nodes, i))
            #print(combo)
            for item in combo:
                #print(combo)
                #print(set(node_list_n).issubset(RO.find_DOIs(G_expanded,list(item))+list(item)))
                #print(redundancy(list(item),drivers))
                if set(node_list_n).issubset(set(RO.find_DOIs(G_expanded,list(item)+driver_previous_motif)+list(item))) \
                        and redundancy(list(item),drivers)==False:
                    #drivers.append(list(item))
                    drivers=list(item)
                    #print(drivers)
                    break


    #print('drivers')
    #print(drivers)
    return drivers




def find_previous_motifs(G_rel_node,G_rel):



    length=G_rel_node[1]['number of sms']
    pre=[]
    driver_p=[]
    for node in G_rel.nodes(data=True)[::-1]:
        if node[0] in G_rel_node[0] and node[1]['number of sms']<length:
            pre.append(node[0])
            driver_p+=node[1]['driver']
            #length=node[1]['number of sms']

    #print('driver pre')
    #print(driver_p)
    return pre,list(set(driver_p))






def stabilized_nodes(nodes_to_find_downstream_of,rules):
    """
    Returns the stabilized nodes after stabilizing a specific motif/motif group.

    Keyword arguments:
        nodes_to_find_downstream_of -- the motif dictionary in pl-po notation, e.g., {pl_1p:0, po_1p:0}
        rules -- Boolean rules

    Returns:
        constants -- stabilized nodes after stabilizing a specific motif/motif group
    """


    node_substitutions = {}
    for n, v in nodes_to_find_downstream_of.items():
        node_substitutions[n] = str(v)
    new_lines = []
    for line in rules.strip().split('\n'):
        n, rule = line.split(',')
        for node, value in node_substitutions.items():
            rule = rule.replace(node, str(value))
        new_lines.append(n + ',\t' + rule)
    rules_new = '\n'.join(new_lines)

    primes = PyBoolNet.FileExchange.bnet2primes(rules_new)
    constants = PyBoolNet.PrimeImplicants.percolate_and_remove_constants(primes)

    #return primes
    return constants


def one_way_DOI(motif,list,G_rel):

    for item in list:
        if item!=motif and G_rel.has_edge(item[0], motif[0]):
                if G_rel[item[0]][motif[0]]['relationship']=='DOI':
                    return True

    return False




def find_control_sets(G_rel,rules):



    lines=rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''
    G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)



    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node


    length_list=[]
    for item in G_rel.nodes(data=True):
        length_list.append(item[1]['number of sms'])

    max_number_of_motifs_in_a_group=max(length_list)


    for i in range (1,max_number_of_motifs_in_a_group+1):
        #print('here')
        for G_rel_node in G_rel.nodes(data=True):

            if G_rel_node[1]['number of sms']==i and i==1:
                driver=find_driver(G_rel_node[0], max_size=2, G_expanded=G_expanded, mapping=mapping)
                nx.set_node_attributes(G_rel, 'driver',{G_rel_node[0]: driver})

            elif G_rel_node[1]['number of sms']==i and i!=1:

                previous_motif,driver_previous_motif=find_previous_motifs(G_rel_node,G_rel)
                current_group = G_rel_node[0]
                for item in previous_motif:
                    if item in previous_motif:
                        current_group=current_group.replace(item,"")

                driver = find_driver(current_group, max_size=2, G_expanded=G_expanded,
                                     mapping=mapping,driver_previous_motif=driver_previous_motif)
                full_driver=driver_previous_motif+driver
                nx.set_node_attributes(G_rel, 'driver' ,{G_rel_node[0]: full_driver})



    #print(G_rel.nodes(data=True))
    return G_rel





def find_control_solutions(attractors,G_rel,target,rules):


    lines=rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node


    sols=[]
    for attr in attractors:
        full_gp_string=''
        for item in attr:
            full_gp_string+=item[0]

        dict_gp=NO.turn_string_to_dict_pl_po(full_gp_string)



        constants=stabilized_nodes(dict_gp, rules)

        if set(NO.find_nodes_in_this_motif(target,mapping)).\
                issubset(set((NO.find_nodes_in_this_motif({**dict_gp, **constants},mapping)))):
            sols.append(attr)




    controls=[]
    for s in sols: #each sol is a consistent group of motifs/motif groups
        control_s=[]
        for item in s: #each item is a motif/motif group
            if one_way_DOI(item,s,G_rel)==False:
                control_s+=item[1]['driver']

        controls.append(control_s) #this is in number index



    #print('\ncontrol set for the target:')
    controls_pl_po=[]
    for item in controls:
        temp=[]
        for n in item:
            temp+=[inverse_mapping[n]]
        controls_pl_po.append(temp)
    #print(controls_pl_po)

    return controls



def find_max_rich_community(mints):


    count_ones=[]
    for item in mints:
        count_ones.append(list(item.values()).count(1))



    index_max_rich=count_ones.index(max(count_ones))

    return(mints[index_max_rich])




def find_min_rich_community(mints):

    count_zeros=[]
    for item in mints:
        count_zeros.append(list(item.values()).count(0))


    index_max_rich=count_zeros.index(max(count_zeros))

    return(mints[index_max_rich])



def calculate_damage_percentage(max_rich_attr,new_attractor,size):


    count_ones_original=list(max_rich_attr.values()).count(1)
    count_ones_damaged=list(new_attractor.values()).count(1)

    return (abs((count_ones_damaged+size)-count_ones_original)/count_ones_original)*100




def count_zeros(attr):

    return  list(attr.values()).count(0)





def in_motifs(node,node_p, G_rel):


    for item in G_rel.nodes(data=True):
        if node_p in item[0] or '~'+node in item[1]['driver']:
            return True
        else:
            continue

    return False







def extinction(G_rel,max_rich_attr,candidates_for_extinction,rules,mode):


    lines = rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node



    #in the case centrality is needed

    dict_centrality = nx.betweenness_centrality(Gread)
    #print(dict_centrality)






    #for item in mapping:
        #print('%s : %s' % (mapping[item], item))


    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')



    size_of_interest=4



    if mode=='maximum':
        motif_candidates_for_maximum_damage=[]
        for node in G_rel.nodes(data=True):
            if '=1' not in node[0] and len(node[1]['driver'])==size_of_interest:
                for ex in node[1]['driver']:
                    if inverse_mapping[ex].replace('~','') in candidates_for_extinction:
                        motif_candidates_for_maximum_damage.append(node)
                        break

        all_per = []
        centrality=[]
        for item in motif_candidates_for_maximum_damage:
            temp={}
            for x in item[1]['driver']:
                temp.update({inverse_mapping[x].replace('~',''): 0})
                centrality.append(dict_centrality[int(x.replace('n','').replace('~',''))])




            constants = stabilized_nodes(temp, rules)
            new_attractor = {**max_rich_attr, **{**constants,**temp}}

            damage_percentage = calculate_damage_percentage(max_rich_attr, new_attractor,size_of_interest)

            print(item[1]['driver'])
            print(damage_percentage)

            all_per.append(damage_percentage)

        if len(motif_candidates_for_maximum_damage)==0:
            return 'No',None
        else:
            return sum(all_per) / len(all_per),sum(centrality)/len(centrality)




    elif mode=='minimum':
        node_candidates_for_minimum_damage=[]
        G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)
        for n in G_expanded.nodes():
            if '_' not in n:
                if len(G_expanded.edges(n))==0 and inverse_mapping[n].replace('~','') in candidates_for_extinction:
                    node_candidates_for_minimum_damage.append(n)

        if node_candidates_for_minimum_damage!=[]:
            all_per=[]
            for node in node_candidates_for_minimum_damage:
                temp={inverse_mapping[node].replace('~',''):0}

                constants=stabilized_nodes(temp,rules)
                new_attractor = {**max_rich_attr, **{**constants,**temp}}  #check to see if this works
                damage_percentage=calculate_damage_percentage(max_rich_attr,new_attractor,size_of_interest)
                all_per.append(damage_percentage)


            return sum(all_per)/len(all_per)
        else:
            return 0



    elif mode=='random_non_stable_motif':
        node_candidates_for_random_damage=[]
        G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)
        for n in G_expanded.nodes():
            if '_' not in n:
                if inverse_mapping[n] in candidates_for_extinction:
                    if in_motifs(n, inverse_mapping[n],G_rel)==False:
                        node_candidates_for_random_damage.append(n)


        if len(node_candidates_for_random_damage)>=size_of_interest:
            combo = list(itertools.combinations(node_candidates_for_random_damage, size_of_interest))

            all_per=[]
            centrality=[]
            for item in combo:
                temp={}
                for node in item:
                    temp.update({inverse_mapping[node].replace('~',''):0}) #i don't think this is right
                    centrality.append(dict_centrality[int(node.replace('n', '').replace('~', ''))])




                constants=stabilized_nodes(temp,rules)
                new_attractor = {**max_rich_attr, **{**constants,**temp}}  #check to see if this works
                damage_percentage=calculate_damage_percentage(max_rich_attr,new_attractor,size_of_interest)
                all_per.append(damage_percentage)


            return sum(all_per)/len(all_per),sum(centrality)/len(centrality)
        else:
            return 'No',None




    else: #partial

        motif_candidates_for_partial_damage=[]
        for node in G_rel.nodes(data=True):
            if '0' in node[0] and len(node[1]['driver'])==size_of_interest:
                for ex in node[1]['driver']:
                    if inverse_mapping[ex].replace('~','') in candidates_for_extinction:
                        motif_candidates_for_partial_damage.append(node)
                        break



        G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)



        all_per = []
        for item in motif_candidates_for_partial_damage:
            temp = {}
            list_inactive_in_control=[x for x in item[1]['driver'] if '~' in x]
            select_some=random.sample(list_inactive_in_control,int(len(list_inactive_in_control)/2))





            for n in G_expanded.nodes():
                if '_' not in n:
                    if inverse_mapping[n] in candidates_for_extinction:
                        if in_motifs(n, inverse_mapping[n], G_rel) == False and len(select_some)<size_of_interest:
                            select_some+=[n]






            print('select some')
            print(select_some)

            for x in select_some:  # this may need if len(select_some)!=0:
                temp.update({inverse_mapping[x].replace('~',''): 0})


            #print(temp)


            constants = stabilized_nodes(temp, rules)
            new_attractor = {**max_rich_attr, **{**constants, **temp}}  # check to see if this works
            # print(temp)
            #print(new_attractor)


            damage_percentage = calculate_damage_percentage(max_rich_attr, new_attractor,size_of_interest)

            all_per.append(damage_percentage)

        if len(motif_candidates_for_partial_damage)==0:
            return 'No', None
        else:
            return sum(all_per) / len(all_per), None









def extinction_after_block(G_rel,max_rich_attr,candidates_for_extinction,rules,mode,block_size,block,size_of_extinction):


    lines = rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node


    #for item in mapping:
        #print('%s : %s' % (mapping[item], item))


    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')



    drivers_of_inactive=[]
    for g_rel_node in G_rel.nodes(data=True):
        if '=0' in g_rel_node[0] and len(g_rel_node[1]['driver'])<=block_size:
            print(g_rel_node)
            drivers_of_inactive+=g_rel_node[1]['driver']

    drivers_of_inactive=list(set(drivers_of_inactive))
    print('drivers of inactive')
    print(drivers_of_inactive)

    if mode=='random':
        node_candidates_for_random_damage=[]
        G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)
        for n in G_expanded.nodes():
            if '_' not in n and block==True:
                if inverse_mapping[n] in candidates_for_extinction and '~'+n not in drivers_of_inactive:
                        node_candidates_for_random_damage.append(n)

            elif '_' not in n and block==False:
                if inverse_mapping[n] in candidates_for_extinction:
                        node_candidates_for_random_damage.append(n)

        print('node_candidates_for_random_damage')
        print(node_candidates_for_random_damage)
        if len(node_candidates_for_random_damage)>=size_of_extinction:
            combo = list(itertools.combinations(node_candidates_for_random_damage, size_of_extinction))

            all_per=[]
            for item in combo:
                temp={}
                for node in item:
                    temp.update({inverse_mapping[node].replace('~',''):0})




                constants=stabilized_nodes(temp,rules)
                new_attractor = {**max_rich_attr, **{**constants,**temp}}  #check to see if this works
                damage_percentage=calculate_damage_percentage(max_rich_attr,new_attractor,size_of_extinction)
                all_per.append(damage_percentage)


            return sum(all_per)/len(all_per)
        else:
            return 'No'



def extinction_block_collapsed(G_rel,max_rich_attr,candidates_for_extinction,to_be_blocked,rules,block):


    lines = rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node


    #for item in mapping:
        #print('%s : %s' % (mapping[item], item))


    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')




    node_candidates_for_random_damage=[]
    G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)
    for n in G_expanded.nodes():
        if '_' not in n and block==True:
            if inverse_mapping[n] in candidates_for_extinction and '~'+n not in to_be_blocked[0]: #actually this needs a function to go through all solutions
                    node_candidates_for_random_damage.append(n)

        elif '_' not in n and block==False:
            if inverse_mapping[n] in candidates_for_extinction:
                    node_candidates_for_random_damage.append(n)

    print('node_candidates_for_random_damage')
    print(node_candidates_for_random_damage)
    if len(node_candidates_for_random_damage)>=len(to_be_blocked):
        combo = list(itertools.combinations(node_candidates_for_random_damage, len(to_be_blocked)))

        all_per=[]
        for item in combo:
            temp={}
            for node in item:
                temp.update({inverse_mapping[node].replace('~',''):0})




            constants=stabilized_nodes(temp,rules)
            new_attractor = {**max_rich_attr, **{**constants,**temp}}  #check to see if this works
            damage_percentage=calculate_damage_percentage(max_rich_attr,new_attractor,len(to_be_blocked))
            all_per.append(damage_percentage)


        return sum(all_per)/len(all_per)

    else: return 'No'




def extinction_block_largest_sm(G_rel,max_rich_attr,candidates_for_extinction,rules,block):


    lines = rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node


    #for item in mapping:
        #print('%s : %s' % (mapping[item], item))


    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')


    max=0
    for g_rel_node in G_rel.nodes(data=True):
        if '=0' in g_rel_node[0] and g_rel_node[0].count('=0')>max:
            largest_string,driver_largest=g_rel_node[0],g_rel_node[1]['driver']
            max=g_rel_node[0].count('=0')


    #print('largest driver')
    #print(driver_largest)
    largest_n=NO.turn_string_to_list_n(largest_string,mapping)



    node_candidates_for_random_damage=[]

    G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)
    for n in G_expanded.nodes():
        if '_' not in n and block==True:
            if inverse_mapping[n] in candidates_for_extinction and '~'+n not in driver_largest:
                node_candidates_for_random_damage.append(n)

        elif '_' not in n and block==False:
            if inverse_mapping[n] in candidates_for_extinction:
                    node_candidates_for_random_damage.append(n)

    #print('largest inactive sm')
    #print(largest_n)

    print('largest inactive sm driver')
    print(driver_largest)


    temp={}
    for node in driver_largest:
        temp.update({inverse_mapping[node].replace('~', ''): 0})

    constants_most = stabilized_nodes(temp, rules)
    new_attractor_most = {**max_rich_attr, **{**constants_most, **temp}}  # check to see if this works
    damage_percentage_most = calculate_damage_percentage(max_rich_attr, new_attractor_most, len(driver_largest))
    print('damage percentage most')
    print(damage_percentage_most)



    print('node_candidates_for_random_damage')
    print(node_candidates_for_random_damage)
    if len(node_candidates_for_random_damage)>=len(driver_largest):
        combo = list(itertools.combinations(node_candidates_for_random_damage, len(driver_largest)))

        all_per=[]
        bool_another_motif=[]
        for item in combo:
            if block==True:
                if set(largest_n).issubset(set(RO.find_DOIs(G_expanded, ['~' + x for x in item]))):

                    continue
                else: #to block the constituent motifs too
                    bool=False
                    for g_rel_node in G_rel.nodes(data=True):
                        if g_rel_node[0] in largest_string and \
                                set(NO.turn_string_to_list_n(g_rel_node[0],mapping)).issubset(set(RO.find_DOIs(G_expanded, ['~' + x for x in item]))):
                            bool=True
                            break
                    if bool==True:

                        continue


            temp={}
            for node in item:
                temp.update({inverse_mapping[node].replace('~',''):0})



            #print('\n')
            #print(temp)

            constants=stabilized_nodes(temp,rules)
            new_attractor = {**max_rich_attr, **{**constants,**temp}}  #check to see if this works
            damage_percentage=calculate_damage_percentage(max_rich_attr,new_attractor,len(driver_largest))
            if damage_percentage!=0:
                all_per.append(damage_percentage)




            if damage_percentage >= 0.8 *damage_percentage_most: #to find if there is an inactive
                # motif that leads to 0.9 of the effect of the largest inactive motif.

                for g_rel_node in G_rel.nodes(data=True):
                    if '=0' in g_rel_node[0] and set(NO.turn_string_to_list_n(g_rel_node[0], mapping)).issubset(
                                set(RO.find_DOIs(G_expanded, ['~' + x for x in item]))):
                        bool_another_motif.append(True)
                        print('the next most damage motif drivers')
                        print(item)
                        print(damage_percentage)
                        break





        if len(all_per)!=0:
            return sum(all_per)/len(all_per),bool_another_motif

    else: return 'No'

    return 'No'




def extinction_block_all_inactive(G_rel,candidates_for_extinction,rules):


    lines = rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node


    for item in mapping:
        print('%s : %s' % (mapping[item], item))


    node_candidates_for_damage=[]

    G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)
    for species in candidates_for_extinction:
        node_candidates_for_damage.append(mapping[species])

    print('node candidates for damage')
    print(node_candidates_for_damage)
    #find the largest driver 

    max_size_driver=0
    for n1 in  G_rel.nodes(data=True):
        if len(n1[1]['driver'])>max_size_driver:
            max_size_driver=len(n1[1]['driver'])





    inative_motifs=[x for x in G_rel.nodes(data=True) if '=0' in x[0]]
    to_be_blocked_all=[]
    combo=[]
    for i in range (1,max_size_driver+1):
        combo += list(itertools.combinations(node_candidates_for_damage, i))




    for g_rel_node in inative_motifs:
        if '=0' in g_rel_node[0]:
            to_be_blocked=[]
            for item in combo:
                if set(NO.turn_string_to_list_n(g_rel_node[0],mapping)).issubset(set(RO.find_DOIs(G_expanded, ['~' + x for x in item]))):
                    #print('motif')
                    #print(set(NO.turn_string_to_list_n(g_rel_node[0],mapping)))
                    #print('to be blocked')
                    #print(to_be_blocked)
                    bool_loop=False
                    for c in to_be_blocked:
                        if set(c).issubset(set(item)):

                            #print('here')
                            #print(c)
                            #print(item)
                            bool_loop=True
                            break

                    if bool_loop==False:
                        to_be_blocked.append(item)

            to_be_blocked_all.append(list(set(to_be_blocked)))



    print('to be blocked')
    print(list(set(list(itertools.chain(*to_be_blocked_all)))))

    final_block=[]
    for item in list(set(list(itertools.chain(*to_be_blocked_all)))):
        for i in item:
            if i not in final_block:
                final_block.append(i)
    #print(list(set(to_be_blocked)))

    print('final block')
    print(final_block)
    return final_block







def extinction_block_parent_size(G_rel):



    G_rel = RO.merge_mutual_DOIs(G_rel)


    parent_inactive=[]
    for g_rel_node in G_rel.nodes(data=True):
        if '=0' in g_rel_node[0] and len(g_rel_node[1].keys())==4:
            print('here')
            print(g_rel_node)
            parent_inactive=g_rel_node


    if parent_inactive!=[]:
        return len(parent_inactive[1]['driver'])

    else:
        return 'No'




def central_extinction(max_rich_attr,candidates_for_extinction,rules,size_of_extinction):


    lines = rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node






    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')


    combo = list(itertools.combinations(candidates_for_extinction, size_of_extinction))


    all_per = []

    for item in combo:
        temp = {}
        for node in item: temp.update({node: 0})



        constants = stabilized_nodes(temp, rules)
        new_attractor = {**max_rich_attr, **{**constants, **temp}}


        damage_percentage = calculate_damage_percentage(max_rich_attr, new_attractor,size_of_extinction)


        all_per.append(damage_percentage)

    return sum(all_per) / len(all_per)



def calculate_restoration_percentage(most_damaged_community,max_rich_attr, new_attractor,size):

    count_ones_damaged = list(most_damaged_community.values()).count(1)
    count_ones_restored = list(new_attractor.values()).count(1)
    count_ones_max_rich = list(max_rich_attr.values()).count(1)

    return (abs((count_ones_restored - size) - count_ones_damaged) / (count_ones_max_rich-count_ones_damaged)) * 100








def restoration(G_rel,max_rich_attr,candidates_for_extinction,rules,size_of_restoration,mode):





    lines = rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node


    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')


    #combo = list(itertools.combinations(candidates_for_extinction, size_of_restoration))




    #damage happens
    motif_candidates_for_sm_damage = []
    for node in G_rel.nodes(data=True):
        if '=1' not in node[0] and len(node[1]['driver']) == size_of_restoration:
            for ex in node[1]['driver']:
                if inverse_mapping[ex].replace('~', '') in candidates_for_extinction:
                    motif_candidates_for_sm_damage.append(node)
                    break

    if len(motif_candidates_for_sm_damage) == 0:
        return 'No'

    else:

        max = 0
        most_damaged_community = {}
        for item in motif_candidates_for_sm_damage:
            temp = {}
            for x in item[1]['driver']:
                temp.update({inverse_mapping[x].replace('~', ''): 0})

            constants = stabilized_nodes(temp, rules)
            new_attractor = {**max_rich_attr, **{**constants, **temp}}
            if max < count_zeros(new_attractor):
                max = count_zeros(new_attractor)
                most_damaged_community = new_attractor
                # damage_percentage = calculate_damage_percentage(max_rich_attr, new_attractor,size_of_restoration)


    #---------------------

    if mode=='sm':


        #restoration begins
        candidates_for_restoration = [] #the ones that die AFTER extinction AKA reintroduction
        for key in max_rich_attr.keys():
            if max_rich_attr[key]==1 and most_damaged_community[key]==0:
                candidates_for_restoration.append(key)

        print('\n')
        print('most damaged community')
        print(most_damaged_community)
        print('candidates for restoration')
        print(candidates_for_restoration)

        motif_candidates_for_sm_restoration = []

        for node in G_rel.nodes(data=True):
            if '=0' not in node[0] and len(node[1]['driver']) == 1:
                for ex in node[1]['driver']:
                    if inverse_mapping[ex].replace('~', '') in candidates_for_restoration:
                        motif_candidates_for_sm_restoration.append(node)
                        break



        print('motif candidates for restoration')
        print(motif_candidates_for_sm_restoration)
        if len(motif_candidates_for_sm_restoration)<size_of_restoration:
            return 'No'

        else:

            combo_to_be_explored = list(itertools.combinations(motif_candidates_for_sm_restoration, size_of_restoration))

            all_per=[]
            for item in combo_to_be_explored:
                temp={}
                for active_motif in item:
                    for x in active_motif[1]['driver']:
                        temp.update({inverse_mapping[x].replace('~',''): 1})
                #print('temp')
                #print(temp)
                constants = stabilized_nodes(temp, rules)
                new_attractor = {**most_damaged_community, **{**constants,**temp}}

                restoration_percentage = calculate_restoration_percentage(most_damaged_community,
                                                                          max_rich_attr, new_attractor,size_of_restoration)

                all_per.append(restoration_percentage)


            return sum(all_per) / len(all_per)




    elif mode=='random':

        node_candidates_for_random_restoration=[]
        G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)

        for key in max_rich_attr.keys():
            if max_rich_attr[key] == 1 and most_damaged_community[key] == 0:
                if in_motifs(mapping[key], key,G_rel)==False:
                    node_candidates_for_random_restoration.append(mapping[key])




        for node in G_expanded:
            if '~' not in node and '_' not in node:
                if inverse_mapping[node] not in max_rich_attr.keys() and \
                        inverse_mapping[node] not in most_damaged_community.keys(): #oscillatories
                    node_candidates_for_random_restoration.append(node)



        if len(node_candidates_for_random_restoration)>=size_of_restoration:
            combo = list(itertools.combinations(node_candidates_for_random_restoration, size_of_restoration))

            all_per=[]
            #centrality=[]
            for item in combo:
                temp={}
                for node in item:
                    temp.update({inverse_mapping[node].replace('~',''):1}) #i don't think this is right
                    #centrality.append(dict_centrality[int(node.replace('n', '').replace('~', ''))])



                constants=stabilized_nodes(temp,rules)
                new_attractor = {**most_damaged_community, **{**constants,**temp}}  #check to see if this works
                restoration_percentage = calculate_restoration_percentage(most_damaged_community,
                                                                          max_rich_attr, new_attractor,size_of_restoration)

                all_per.append(restoration_percentage)

            return sum(all_per)/len(all_per)
        else:
            return 'No'





def restoration_FD(G_rel,max_rich_attr,candidates_for_extinction,rules,size_of_restoration,cluster):





    lines = rules.splitlines(True)
    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)
    prefix, suffix = 'n', ''

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    for i, node in enumerate(readnodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node


    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')





    #----------damage happens------------------
    motif_candidates_for_sm_damage = []
    for node in G_rel.nodes(data=True):
        if '=1' not in node[0] and len(node[1]['driver']) == size_of_restoration:
            for ex in node[1]['driver']:
                if inverse_mapping[ex].replace('~', '') in candidates_for_extinction:
                    motif_candidates_for_sm_damage.append(node)
                    break

    if len(motif_candidates_for_sm_damage) == 0:
        return 'No'

    else:

        max = 0
        most_damaged_community = {}
        for item in motif_candidates_for_sm_damage:
            temp = {}
            for x in item[1]['driver']:
                temp.update({inverse_mapping[x].replace('~', ''): 0})

            constants = stabilized_nodes(temp, rules)
            new_attractor = {**max_rich_attr, **{**constants, **temp}}
            if max < count_zeros(new_attractor):
                max = count_zeros(new_attractor)
                most_damaged_community = new_attractor
                # damage_percentage = calculate_damage_percentage(max_rich_attr, new_attractor,size_of_restoration)

    #----------------------------------


    #find out what species survived

    survived_the_damage = []  # the ones that die AFTER extinction AKA reintroduction
    for key in most_damaged_community.keys():
        if most_damaged_community[key] == 1:
            survived_the_damage.append(key)


    #print('most damaged community')
    #print(most_damaged_community)



    #restoration begins
    candidates_for_restoration = [] #the ones that die AFTER extinction AKA reintroduction
    for key in max_rich_attr.keys():
        if max_rich_attr[key]==1 and most_damaged_community[key]==0:
            candidates_for_restoration.append(key)

    #print('candidates for restoration')
    #print(candidates_for_restoration)


    bi_adj_matrix, cluster_list = SO.form_bi_adj_matrix(Gread, inverse_mapping, cluster=cluster)


    if len(cluster_list)<2: #there cannot be clusters
        return 'No'




    # map the survived to the index in the linkage format
    if cluster=='pollinator':
        survived_index=[]
        for k in [x for x in survived_the_damage if 'po' in x]:
            if Gread.out_degree(int(mapping[k].replace('n','')))!=0:
                survived_index.append(cluster_list.index(int(mapping[k].replace('n',''))))

    else:
        survived_index = []
        for k in [x for x in survived_the_damage if 'pl' in x]:
            if Gread.out_degree(mapping[k].replace('n', '')) != 0:
                survived_index.append(cluster_list.index(int(mapping[k].replace('n', ''))))


    #print('survived index in linkage numbering')
    #print(survived_index)

    dict_added_FD = SO.restoration_functional_diversity(bi_adj_matrix, cluster_list, survived=survived_index,
                                                        size_of_restoration=size_of_restoration)

    cut_off=sorted(dict_added_FD.values())[-1*(int(len(dict_added_FD)*0.1)+1)]
    to_come_back=[x for x in dict_added_FD.keys() if dict_added_FD[x]>=cut_off]


    all_per=[]
    for item in to_come_back: #to come back is in linkage index


            temp={}
            bool=[]
            for x in item:
                if inverse_mapping['n'+str(cluster_list[x])] in candidates_for_restoration:
                    bool.append(True)
                else:
                    bool.append(False)
            if all(bool):
                for y in item:
                    temp.update({inverse_mapping['n'+str(cluster_list[y])]: 1})


            if temp=={}:
                continue
            #print('temp')
            #print(temp)
            constants = stabilized_nodes(temp, rules)
            new_attractor = {**most_damaged_community, **{**constants,**temp}}
            #print('new attractor')
            #print(new_attractor)
            restoration_percentage = calculate_restoration_percentage(most_damaged_community,
                                                                      max_rich_attr, new_attractor,size_of_restoration)

            all_per.append(restoration_percentage)


    if all_per!=[]:
        return sum(all_per) / len(all_per)
    else:
        return 'No'




