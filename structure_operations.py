import networkx as nx
from numpy import linalg as LA
import math
import numpy as np
import control_operations as CO
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt
from iteration_utilities import deepflatten






def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper






def node_nestedness(G):


    """
    The original network nodes should be in numbers like 1,2,3 (the output of Colin's form_network())
    Calculates the node nestedness based on Samuel Jonhson's 2013 paper
    """


    ad_matrix=nx.adjacency_matrix(G,nodelist=[x for x in range(0,len(G.nodes()))])
    ad_matrix_sq=LA.matrix_power(ad_matrix.todense().transpose(), 2)
    in_degree_dict=G.in_degree(G.nodes())


    dict_node_nestedness={}
    for i in range(len(G.nodes())):
        sum = 0
        if in_degree_dict[i] != 0:
            for j in range (len(G.nodes())):
                if j!=i:
                    #print('in the loop')
                    #print(i)
                    #print(j)
                    #print(sum)
                    #print(ad_matrix_sq[i,j])
                    #print(in_degree_dict[i])
                    #print(in_degree_dict[j])
                    if in_degree_dict[j]!=0:
                        sum+= (ad_matrix_sq[i,j]/(in_degree_dict[i]*in_degree_dict[j]))


        dict_node_nestedness.update({i:sum/(len(G.nodes()))})


    return dict_node_nestedness



def node_mus_rank_plpo(G,inverse_mapping,false_nodes):


    # -----separating plants and pollinators
    prefix='n'
    plants=[]
    pollinators=[]


    for node in G.nodes():

        #MUR rev
        if 'pl' in inverse_mapping[prefix+str(node)] and inverse_mapping[prefix+str(node)] not in false_nodes:
            plants.append(node)
        elif 'po' in inverse_mapping[prefix+str(node)] and inverse_mapping[prefix+str(node)] not in false_nodes:
            pollinators.append(node)


    #print('pollinators')
    #print(pollinators)



    def converged(iterations): #needs to be fixed

        bool_1=0
        bool_2=0
        for i in plants+pollinators:
            if iterations[-1][i]==iterations[-2][i]:
                bool_1+=1
            if iterations[-1][i]==iterations[-3][i]:
                bool_2+=1

        #print('here')
        #print(len(plants))
        #print(len(pollinators))
        #print(bool_1)
        #print(bool_2)

        if bool_1==len(plants)+len(pollinators) or bool_2==len(plants)+len(pollinators):
            return True
        else:
            return False





    #calculating importance and vulnerability
    iterations=[]
    IC={}
    for item in plants+pollinators:
        IC.update({item:1})

    iterations.append(IC)


    while (1):
        this_step = {}


        if len(iterations)>=3 and converged(iterations)==True:
                break


        #calculating importances
        temp_po = {}
        for A in pollinators: #in 1,2,3
            I=0
            #for edge in G.out_edges(A):
            for edge in G.edges(A):
                #print('A')
                #print(A)
                #print('edge')
                #print(edge)
                I+=iterations[-1][edge[1]]
                #print('I')
                #print(I)


            temp_po.update({A:I})






        av=sum(temp_po.values())/len(temp_po.values())

        for case in temp_po.keys():
            this_step.update({case:truncate(temp_po[case]/av,3)})


        if 0 in this_step.values(): break


        #calculating the vulnerabilities
        temp_pl = {}
        for P in plants:
            denominator=0
            #print('here')
            #print(P)
            #print(G.edges(P))
            #for edge1 in G.in_edges(P):
            for edge1 in G.edges(P):


                denominator += 1/(iterations[-1][edge1[1]])

            V=1/denominator
            temp_pl.update({P:V})

        av = sum(temp_pl.values()) / len(temp_pl.values())
        for case1 in temp_pl.keys():
            this_step.update({case1: truncate(temp_pl[case1] / av,3)})



        iterations.append(this_step)
        #print('iteration this step')
        #print(this_step)

    importances={}
    for item in iterations[-1].keys():
        if item in pollinators:
            importances.update({item:iterations[-1][item]})

    #print(iterations[-1])
    return importances






def robustness(G,rules,inverse_mapping,false_nodes):




    rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
        ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES', '')



    plants=[]
    pollinators=[]
    prefix='n'

    false_plants=[]
    false_pollinators=[]
    for item in false_nodes:
        if 'pl' in item: false_plants.append(item)
        else: false_pollinators.append(item)




    for node in G.nodes():


        if 'pl' in inverse_mapping[prefix+str(node)] and inverse_mapping[prefix+str(node)] not in false_nodes:
            plants.append(node)
        elif 'po' in inverse_mapping[prefix+str(node)] and inverse_mapping[prefix+str(node)] not in false_nodes:
            pollinators.append(node)




    out_degree_pollinators_dict = G.out_degree(pollinators)
    out_degree_po=dict(reversed(list({k: v for k, v in sorted(out_degree_pollinators_dict.items(), key=lambda item: item[1])}.items())))
    #print(out_degree_po)


    extinct={}
    ATC_data=[]
    for key in out_degree_po.keys():
        extinct.update({inverse_mapping[prefix+str(key)]:0})
        #print(extinct)
        constants = CO.stabilized_nodes(extinct, rules)
        #print(constants)
        following_extinction={**{**constants, **extinct}}
        secondary_extinct_plants=[]
        for key in following_extinction.keys():
            if 'pl' in key and key not in false_plants:
                #print('here')
                secondary_extinct_plants.append(key)

        #print(secondary_extinct_plants)
        fraction_plants_to_survive=(len(plants)-len(false_plants)-len(secondary_extinct_plants))/(len(plants)-len(false_plants))
        #print(fraction_plants_to_survive)
        if fraction_plants_to_survive>=0:
            ATC_data.append([(len(extinct.keys())/(len(pollinators)-len(false_pollinators))),fraction_plants_to_survive])
        if fraction_plants_to_survive==0:
            break

        if fraction_plants_to_survive<0:
            break

    print(ATC_data)






def form_bi_adj_matrix(G,inverse_mapping,cluster):



    plants=[]
    pollinators=[]

    #print(G.nodes())

    for n in G.nodes():
        if 'pl' in inverse_mapping['n'+str(n)] and G.in_degree(n)!=0:
            plants.append(n)
        elif 'po' in inverse_mapping['n'+str(n)] and G.in_degree(n)!=0 and G.out_degree(n)!=0:
            pollinators.append(n)

    #print(len(plants))
    #print(len(pollinators))

    if cluster=='pollinator':

        adj_matrix=np.zeros((len(plants),len(pollinators)))
        #print(adj_matrix)

        for i in range(len(plants)):
            for j in range(len(pollinators)):
                if (pollinators[j],plants[i]) in G.edges():
                    adj_matrix[i][j]=1


        print([inverse_mapping['n'+str(x)] for x in plants])
        print([inverse_mapping['n' + str(x)] for x in pollinators])
        #print([sum(x) for x in adj_matrix])

        return adj_matrix.transpose() ,pollinators
    else:

        adj_matrix = np.zeros((len(pollinators), len(plants)))
        # print(adj_matrix)

        for i in range(len(pollinators)):
            for j in range(len(plants)):
                if (plants[j], pollinators[i]) in G.edges():
                    adj_matrix[i][j] = 1

        #print([inverse_mapping['n' + str(x)] for x in plants])
        #print([inverse_mapping['n' + str(x)] for x in pollinators])
        # print([sum(x) for x in adj_matrix])

        return adj_matrix.transpose(), plants



def restoration_functional_diversity(bi_adj_matrix,cluster_list,survived,size_of_restoration=1):


    #print(bi_adj_matrix)
    Z = hierarchy.linkage(bi_adj_matrix, 'single', metric='euclidean')

    #print(Z)
    clusters = []
    n = len(Z) + 1
    len_Z_fix = len(Z) + 1
    for item in Z:

        if item[0] >= len_Z_fix:
            # print(item[0])
            for c in clusters:
                if c[2] == item[0]:
                    entry_1 = list(deepflatten([c[0]] + [c[1]]))
            # print(entry_1)
        else:
            entry_1 = item[0]

        if item[1] >= len_Z_fix:
            # print(item[1])
            for c1 in clusters:
                if c1[2] == item[1]:
                    entry_2 = list(deepflatten([c1[0]] + [c1[1]]))
            # print('entry2')
            # print(entry_2)
        else:
            entry_2 = item[1]

        if item[0] < len_Z_fix and item[1] < len_Z_fix:
            entry_1 = item[0]
            entry_2 = item[1]

        clusters.append([entry_1, entry_2, n, item[2]])
        n += 1

    #for jm in clusters:
        #print(jm)

    plt.figure()



    #levels=sorted(list(set([x[-2] for x in Z])))


    in_between=0
    for z in range(len_Z_fix, len_Z_fix+len(Z)):
        for item in Z:
            #print('\n')
            #print(z)
            #print(item)
            if z in item[0:2] and Z[z-len_Z_fix][-2]!=item[-2]:
                #print('here')
                in_between+=1
                break

    #print('in between')
    #print(in_between)
    #print(len_Z_fix + in_between)


    distribution_matrix = np.zeros((len_Z_fix + in_between, len_Z_fix))

    #print(distribution_matrix)

    paths_label = list(map(chr, range(97, 97 + len(distribution_matrix))))


    paths_length = {}
    i = 0

    for cl in clusters:
        if type(cl[0]) != list and type(cl[1]) != list:
            paths_length.update({paths_label[i]: cl[-1]})
            distribution_matrix[i][int(cl[0])] = 1
            i += 1
            paths_length.update({paths_label[i]: cl[-1]})
            distribution_matrix[i][int(cl[1])] = 1
            i += 1


        elif type(cl[0]) != list and type(cl[1]) == list:
            paths_length.update({paths_label[i]: cl[-1]})
            distribution_matrix[i][int(cl[0])] = 1
            i += 1




    for z in range(len_Z_fix, len_Z_fix+len(Z)): #z=10,11,12
        for item in Z: #item: [ 3, 5, 1, 2 ]
            if z in item[0:2] and Z[z-len_Z_fix][-2]!=item[-2]:
                first_layer = Z[z-len_Z_fix][-2]
                second_layer=item[-2]
                paths_length.update({paths_label[i]: second_layer - first_layer})



                for f in Z[z-len_Z_fix][0:2]:
                    if f < len_Z_fix:
                        distribution_matrix[i][int(f)] = 1
                    else:
                        for index in list(deepflatten(clusters[int(f)-len_Z_fix][0:2])):
                            distribution_matrix[i][int(index)] = 1

                i+=1
                break

    paths_length_matrix = np.zeros(len(paths_length))
    w = 0
    for key in paths_length.keys():
        paths_length_matrix[w] = paths_length[key]
        w += 1

    #print('path length matrix')
    #print(paths_length_matrix)
    #print('path dict')
    #print(paths_length)

    cut_matrix=np.zeros(len(distribution_matrix))
    #cut_matrix=np.add(distribution_matrix[:, survived[0]], distribution_matrix[:, survived[1]])
    for y in survived:
        cut_matrix = np.add(cut_matrix, distribution_matrix[:, y])

    #print('cut matrix')
    #print(cut_matrix)

    a_0_1 = np.zeros(len(paths_length))

    w = 0
    for w in range(len(paths_length_matrix)):
        if cut_matrix[w] != 0:
            a_0_1[w] = 1

    w = 0
    sum = 0
    for w in range(len(paths_length_matrix)):
        sum += paths_length_matrix[w] * a_0_1[w]


    print('FD before restoration')
    print(sum)

    import itertools


    to_be_explored=[x for x in range(len(cluster_list)) if x not in survived]
    combo_to_be_explored=list(itertools.combinations(to_be_explored, size_of_restoration))

    dict_added_FD={}
    for c in combo_to_be_explored:
        cut_restoration=cut_matrix
        for c_n in c:
            cut_restoration = np.add(cut_restoration, distribution_matrix[:, c_n])


        a_0_1_restore = np.zeros(len(paths_length))

        w = 0
        for w in range(len(paths_length_matrix)):
            if cut_restoration[w] != 0:
                a_0_1_restore[w] = 1

        w = 0
        sum_restore = 0
        for w in range(len(paths_length_matrix)):
            sum_restore += paths_length_matrix[w] * a_0_1_restore[w]

        dict_added_FD.update({c:sum_restore-sum})



    print(dict_added_FD)
    #dn = hierarchy.dendrogram(Z)
    #plt.show()
    return dict_added_FD

