
import PyBoolNet
import networkx as nx

import pandas as pd
import math
import numpy as np
import itertools
from itertools import permutations
import random





def mapping(x):
    return x + 'p'





def disjunctive_prime_form_text_file(G, network_number, plant, pollinator, IC_number=0,
                                     IC=[], write_IC=False, perturb=False, wp=0,request_for_p_n_edges=False):

    #writes the Boolean functions of a pl-po threshold model in disjunctive prime form in a text file
    #if request_for_p_n_edges==True it returns the number of remaining positive and negative edges after converting to disjunctive
    #prime form, otherwise doesn't return anything.
    #this is currently hard coded for positive weight of 4 and negative weight of -1 for the Campbell et al. model.

    #inputs:
    #G: the original digraph including all prositive and negative edges. This object is taken from Colin's ensembles
    #network_number: the index of a network within each batch. This index starts at 0 and goes to 999 for each batch.
    #plant: the number of plants in the network
    #pollinator: the number of pollinators in the network
    #write_IC: if True, the function writes initial conditions at the beginning of the text file
    #IC: initial conditions in the form of a list of 0s and 1s showing what the state of each species is at t=0. This is only implemented if one
    #wants to try synchronous update or study perturbations.
    #perturb: if True, the perturbed version of the Boolean functions are written within the text file
    #wp: weight of the perturbation. for example if set to 0.3 in a network of 100 nodes, 30 nodes will be perturbed to their active states.


    G = nx.relabel_nodes(G, mapping)
    nodes = list(G.nodes(data=True))
    edges = list(G.edges(data=True))

    length = len(nodes)
    name_of_nodes = [item[0] for item in nodes]
    false_nodes, name_of_remained_nodes = find_false_nodes(G)
    perturbed_to_be_true = name_of_remained_nodes[0:int((plant + pollinator) * wp)]


    address='disjunctive_prime_form_' + str(plant) + '_' + str(pollinator) + '_' + str(network_number)
    if write_IC==True:
        address+='_IC'+str(IC_number)
    f = open(address+'.txt', 'w')











    if write_IC == True:
        for i in range(len(IC)):
            if IC[i] == 0:
                f.write(name_of_nodes[i] + '=False' + '\n')
            else:
                f.write(name_of_nodes[i] + '=True' + '\n')

        f.write('\n')
    counter_p = 1

    negative_edges_in_this_net=0
    positive_edges_in_this_net = 0



    for j in range(length):
        positive_regulators = []
        negative_regulators = []
        regulators = []
        counter = 0
        sample_node = name_of_nodes[j]
        for i in range(len(edges)):
            if edges[i][1] == sample_node:
                regulators.append(edges[i])

        regulator_signs = [item[2] for item in regulators]
        length_regulators = len(regulators)
        for k in range(length_regulators):
            if regulator_signs[k].get('data') == -1.0:
                negative_regulators.append(regulators[k])
            else:
                check = 0
                for r in range(len(false_nodes)):
                    if regulators[k][0] == false_nodes[r]:
                        check = 1

                if check == 0:
                    positive_regulators.append(regulators[k])


        positive_edges_in_this_net += len(positive_regulators)
        if len(negative_regulators)>3 and len(positive_regulators)!=0:
            negative_edges_in_this_net+=len(negative_regulators)





        name_of_positive_reg = [item[0] for item in positive_regulators]
        name_of_negative_reg = [item[0] for item in negative_regulators]
        char = name_of_nodes[j] + '*=('
        if perturb == True and counter_p <= int((plant + pollinator) * wp):
            bool_p = False
            for a in range(len(perturbed_to_be_true)):
                if name_of_nodes[j] == perturbed_to_be_true[a]:
                    bool_p = True

            if bool_p == True:
                char = name_of_nodes[j] + '*=True'

                counter_p = counter_p + 1
        else:
            for active_pos_reg in range(1, len(positive_regulators) + 1):
                inactive_neg_reg = len(negative_regulators) - (4 * active_pos_reg - 1)
                if inactive_neg_reg <= 0:
                    counter = counter + 1
                if inactive_neg_reg <= 0:
                    if active_pos_reg != 1:
                        if counter > 1:
                            break
                if char != name_of_nodes[j] + '*=(':
                    char = char + ' or ('
                if inactive_neg_reg <= 0 and active_pos_reg == 1:
                    for n in range(len(name_of_positive_reg)):
                        char = char + name_of_positive_reg[n]
                        if n != len(name_of_positive_reg) - 1:
                            char = char + ' or '

                    char = char + ')'
                    break
                else:
                    bool = False
                    if inactive_neg_reg > 0:
                        combo_neg = list(set(itertools.combinations(name_of_negative_reg, inactive_neg_reg)))
                        bool = True
                    combo_pos = list(set(itertools.combinations(name_of_positive_reg, active_pos_reg)))
                    for u in range(len(combo_pos)):
                        for n in range(len(combo_pos[u])):
                            char = char + combo_pos[u][n]
                            if n != len(combo_pos[u]) - 1:
                                char = char + ' and '
                            if u != len(combo_pos) - 1 and n == len(combo_pos[u]) - 1:
                                char = char + ' or '

                    char = char + ')'
                    if bool == True:
                        char = char + ' and ('
                        for a in range(len(combo_neg)):
                            for m in range(len(combo_neg[a])):
                                char = char + 'not ' + combo_neg[a][m]
                                if m != len(combo_neg[a]) - 1:
                                    char = char + ' and '
                                if a != len(combo_neg) - 1 and m == len(combo_neg[a]) - 1:
                                    char = char + ' or '

                        char = char + ')'

        if len(positive_regulators) == 0:
            char = name_of_nodes[j] + '*=' + 'False'

        f.write(char + '\n')

    f.close()

    if request_for_p_n_edges==True:
        return positive_edges_in_this_net,negative_edges_in_this_net





def find_false_nodes(G):

    #finds the species that cannot establish under any circumstances and become extinct in all attractors. Figure 2 on the paper describes these nodes.
    #inputs:
    #G: the original digraph including all prositive and negative edges. This object is taken from Colin's ensembles
    #output:
    #false_nodes: the species that cannot establish in a list



    nodes = list(G.nodes(data=True))
    edges = list(G.edges(data=True))
    length = len(nodes)
    name_of_nodes = [item[0] for item in nodes]
    false_nodes = []
    df = pd.DataFrame(nodes)
    for s in range(length):
        positive_regulators = []
        sample_node = name_of_nodes[s]
        for i in range(len(edges)):
            if edges[i][1] == sample_node:
                sign = edges[i][2].get('data')
                if sign == 1:
                    positive_regulators.append(edges[i])

        len_p = len(positive_regulators)
        if len_p == 0:
            false_nodes.append(sample_node)
            df = df.drop(s)

    memory = 0
    while True:
        nodes = df.values.tolist()
        df = pd.DataFrame(nodes)
        name_of_remained_nodes = [item[0] for item in nodes]
        number_of_remained_nodes = len(name_of_remained_nodes)
        if memory == number_of_remained_nodes:
            break
        else:
            memory = number_of_remained_nodes
            for j in range(number_of_remained_nodes):
                positive_regulators = []
                sample_node = name_of_remained_nodes[j]
                for i in range(len(edges)):
                    if edges[i][1] == sample_node:
                        sign = edges[i][2].get('data')
                        if sign == 1:
                            counter = 0
                            for a in range(len(false_nodes)):
                                if edges[i][0] == false_nodes[a]:
                                    counter = 1

                            if counter == 0:
                                positive_regulators.append(edges[i])

                len_p = len(positive_regulators)
                if len_p == 0:
                    df = df.drop(j)
                    false_nodes.append(name_of_remained_nodes[j])

    return (false_nodes, name_of_remained_nodes)





def simplification_text_file(G,network_number,plant,pollinator,IC=[],IC_number=0,
                             write_IC=False,more=False,less=True,perturb=False,wp=0):


    #simplifies and writes the Boolean functions of a pl-po threshold model in disjunctive prime form in a text file. simplification method is
    # explained in the paper. It samples negative edges instead of keeping all of them and preserves the probability of being active in the traget node.
    #both directions of the inequality discussed in the paper are implemented in this function. 'less' stands for 'keeping less negative edges' in the
    # final Boolean function when p_b>=p_t (the case we decided to go with), and 'more' stands for 'keeping more negative edges' in the
    # final Boolean function when p_b<=p_t
    #this is currently hard coded for positive weight of 4 and negative weight of -1 for the Campbell et al. model.


    #inputs:
    #G: the original digraph including all prositive and negative edges. This object is taken from Colin's ensembles
    #network_number: the index of a network within each batch. This index starts at 0 and goes to 999 for each batch.
    #plant: the number of plants in the network
    #pollinator: the number of pollinators in the network
    #write_IC: if True, the function writes initial conditions at the beginning of the text file
    #IC: initial conditions in the form of a list of 0s and 1s showing what the state of each species is at t=0.
    #This is only implemented if one wants to try synchronous update or perturbations.
    #perturb: if True, the perturbed version of the Boolean functions are written within the text file
    #wp: weight of the perturbation. for example if set to 0.3 in a network of 100 nodes, 30 nodes will be perturbed to their active states.


    def mapping(x):
        return x + 'p'

    max=0
    G = nx.relabel_nodes(G, mapping)
    nodes = list(G.nodes(data=True))
    edges = list(G.edges(data=True))
    length = len(nodes)
    name_of_nodes = [item[0] for item in nodes]



    string='simplified_'+str(plant) + '_' + str(pollinator) + '_' + str(network_number)
    if write_IC==True:
        string+= '_IC'+str(IC_number)





    if more==True and less==True:
        string_more=string+'_more.txt'
        string_less=string+'_less.txt'
        f1 = open(string_more, 'w')
        f2 = open(string_less, 'w')
    elif more==True and less==False:
        string=string+'_more.txt'
        f1 = open(string, 'w')
    elif more==False and less==True:
        string=string+'_less.txt'
        f2 = open(string, 'w')
    false_nodes = []
    df = pd.DataFrame(nodes)

    #f2.write('#BOOLEAN RULES'+'\n') #if one wants to run Jorge's java stable motif code on the output of this function, this line is necessary.



    if write_IC==True:
        for i in range (len(IC)):
            if IC[i]==0:
                if more==True:
                    f1.write(name_of_nodes[i] + '=False'+'\n')
                if less==True:
                    f2.write(name_of_nodes[i] + '=False'+'\n')
            else:
                if more==True:
                    f1.write(name_of_nodes[i] + '=True'+'\n')
                if less==True:
                    f2.write(name_of_nodes[i] + '=True'+'\n')
        if more==True:
            f1.write('\n')
        if less==True:
            f2.write('\n')

    for s in range(length):  # immediate false nodes
        positive_regulators = []
        negative_regulators=[]
        sample_node = name_of_nodes[s]
        for i in range(len(edges)):
            if edges[i][1] == sample_node:
                sign = edges[i][2].get('data')
                if sign == 1:
                    positive_regulators.append(edges[i])
                else:
                    negative_regulators.append(edges[i])
        len_p = len(positive_regulators)
        if len_p == 0:
            false_nodes.append(sample_node)
            df = df.drop(s)

        if len(positive_regulators+negative_regulators)>max:
            max=len(positive_regulators+negative_regulators)

    memory = 0
    while True:  # the rest of the nodes that become false
        nodes = df.values.tolist()
        df = pd.DataFrame(nodes)
        name_of_remained_nodes = [item[0] for item in nodes]
        number_of_remained_nodes = len(name_of_remained_nodes)
        if memory == number_of_remained_nodes:
            break
        else:
            memory = number_of_remained_nodes
            for j in range(number_of_remained_nodes):
                positive_regulators = []
                sample_node = name_of_remained_nodes[j]
                for i in range(len(edges)):
                    if edges[i][1] == sample_node:
                        sign = edges[i][2].get('data')
                        if sign == 1:
                            counter = 0
                            for a in range(len(false_nodes)):
                                if edges[i][0] == false_nodes[a]:
                                    counter = 1
                            if counter == 0:
                                positive_regulators.append(edges[i])
                len_p = len(positive_regulators)
                if len_p == 0:
                    df = df.drop(j)
                    false_nodes.append(name_of_remained_nodes[j])


    print('species that cannot establish')
    print(false_nodes)
    perturbed_to_be_true=name_of_remained_nodes[0:int((plant+pollinator)*wp)]
    nodes = list(G.nodes(data=True))
    length = len(nodes)
    counter_p=1
    for j in range(length):
        regulators = []
        positive_regulators=[]
        all_positive_regulators=[]
        negative_regulators=[]
        sample_node = name_of_nodes[j]
        for i in range(len(edges)):
            if edges[i][1] == sample_node:
                sign= edges[i][2].get('data')
                if sign==1:
                    all_positive_regulators.append(edges[i])
                    counter=0
                    for a in range (len(false_nodes)):
                        if edges[i][0]==false_nodes[a]:
                            counter=1
                    if counter==0:
                        positive_regulators.append(edges[i])
                else:
                    negative_regulators.append(edges[i])

        #print('\n')
        #print(sample_node)
        #print(all_positive_regulators)
        #print(positive_regulators)
        #print(negative_regulators)

        name_of_negative_regulators = [item [0] for item in negative_regulators]
        name_of_positive_regulators = [item[0] for item in positive_regulators]

        if len(positive_regulators) != 0:
            regulators=positive_regulators
        if len(negative_regulators) != 0:
            regulators=regulators+negative_regulators

        len_p=len(positive_regulators)
        len_n=len(negative_regulators)

        char = sample_node + '*='
        if len_p == 0:
            char = char + 'False'
            if more==True:
                f1.write(char+'\n')
            if less==True:
                f2.write(char+'\n')
            false_nodes.append(sample_node)
        else:
            if perturb==True and counter_p<=int((plant+pollinator)*wp):
                bool_p=False
                for z in range (len(perturbed_to_be_true)):
                    if sample_node==perturbed_to_be_true[z]:
                        bool_p=True
                if bool_p==True:
                    counter_p = counter_p + 1
                    if more==True:
                        f1.write(char+'True'+'\n')
                    if less==True:
                        f2.write(char+'True'+'\n')
            else:
                summation = 0
                for k in range(1, len_p + 1):
                    k_of_p = math.factorial(len_p) // math.factorial((len_p - k)) // math.factorial(k)
                    sum = 0
                    if len_n >= 4 * k - 1:
                        limit = 4 * k
                    else:
                        limit = len_n + 1
                    for c in range(limit):
                        sum = sum + (math.factorial(len_n) // math.factorial(len_n - c) // math.factorial(c))
                    sum = sum * k_of_p
                    summation = summation + sum

                sol = len_n - (np.log(summation / ((2 ** len_p) - 1)) / np.log(2))


                x_more = 0
                while True:
                    if x_more >= sol:
                        break
                    else:
                        x_more = x_more + 1

                x_less=len_n
                while True:
                    if x_less <= sol:
                        break
                    else:
                        x_less = x_less - 1

                char=char+'('
                for i1 in range (len_p):
                    if i1!=len_p-1:
                        char=char+name_of_positive_regulators[i1]+ ' or '
                    else:
                        char = char + name_of_positive_regulators[i1] + ') '

                neg_more = random.sample(name_of_negative_regulators, x_more)
                neg_less=neg_more
                char_more=char
                char_less=char
                if x_more!=0 and more==True:

                    char_more=char_more+'and '
                    for i2 in range (len(neg_more)):
                        char_more = char_more + ' not ' + neg_more[i2]
                        if i2 != int(x_more) - 1:
                            char_more = char_more + ' and '
                if x_less!=0 and less==True:
                    if x_more != x_less:
                        diff = x_more - x_less
                        for j in range(diff):
                            neg_less.pop(-1)
                    char_less = char_less + 'and '
                    for i2 in range(len(neg_less)):
                        char_less = char_less + ' not ' + neg_less[i2]
                        if i2 != int(x_less) - 1:
                            char_less = char_less + ' and '


                if len_p==0:
                    if more==True:
                        f1.write(char+'\n')
                    if less==True:
                        f2.write(char+'\n')
                else:
                    if more==True:
                        f1.write(char_more+'\n')
                    if less==True:
                        f2.write(char_less+'\n')

    if more==True:
        f1.close()
    if less==True:
        f2.close()



def number_of_false_nodes(lines):

    #counts the number of nodes that are initially False
    #input: Boolean rules in the format of lines
    #output: returns the number of False nodes in the original Boolean model


    total_count = 0
    for i in range(len(lines)):
        count = lines[i].count('False')
        total_count = total_count + count
    return total_count




def write_g_expanded(G_expanded,model_name):

    #writes expanded network in a gml file that can be opened with yED.
    #inputs:
    #the expanded network as Gang Yang's function Get_expanded_network() produces
    #model_name: name of the model in a string


    nx.write_gml(G_expanded, 'expanded_network_' + model_name+ '.gml')




def number_of_unstabilized_nodes(nodes_to_find_downstream_of,rules,number_of_nodes): #ye argumente akhar dashte watch out


    #counts the number of unstabilized nodes after stabilizing a specific motif/motif group
    #inputs:
    # nodes_to_find_downstream_of: the motif dictionary in pl-po notation {pl_1p:0, po_1p:0}
    # rules: Boolean rules
    # number_of_nodes: total number of nodes in the network
    #output: the number of unstabilized nodes after stabilizing a specific motif/motif group




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

    return (number_of_nodes - len(constants))




def number_of_negative_edges_text(lines):

    #returns the number of negative edges in the Boolean model
    #input: lines which is the Boolean model text file read with the function readlines()
    #output: the number of negative edges
    total_count=0
    for i in range (len(lines)):
        count = lines[i].count('not')
        total_count=total_count+count
    return total_count



