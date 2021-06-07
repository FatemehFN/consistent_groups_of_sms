"""
This script shows the workflow for generating the expanded network and cycle graph of the Boolean model in Figure 4 and Figure 5.
Depending on which file_name the user chooses, the corresponding expanded network and cycle graph are generated and saved in gml files.
First the Boolean model is read from a text file and the network Gread is constructed accordingly. Then this network is fed into the
functions Get_expanded_network() and write_g_expanded to get and write the expanded network respectively. Lastly the Boolean model in lines
is given to function cycle_graph_virtual_node_based() to generate and save the cycle graph.
"""


import model_operations as MO
import BooleanDOI_processing as BDOIp
import sm_csm_operations as SMO




file_name='figure_4_Boolean_functions'
#file_name='figure_5_Boolean_functions'

#read the text file of the Boolean model
f=open('Figures_data/'+file_name+'.txt')
lines=f.readlines()
f.close()

# build the network from the lines
Gread, readnodes = BDOIp.form_network(lines, sorted_nodename=False)

# build the expanded network
prefix, suffix = 'n', ''
G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)

# write the expanded network in a gml file that can be opened with YED
MO.write_g_expanded(G_expanded,model_name=file_name)

# calculate the mapping from original node names to number index
mapping = {}  # nodename to number index
for i, node in enumerate(readnodes):
    index = prefix + str(i) + suffix
    mapping[node] = index
    mapping['~' + node] = '~' + index

print('Labels in the expanded network:')
for item in mapping:
    print('%s : %s'% (mapping[item],item))

#build cycle graph using virtual nodes of expanded network and write it in a gml file that can be opened with YED
cycle_graph=SMO.cycle_graph_virtual_node_based(lines,write_cycle_graph=True)
