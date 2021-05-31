import model_operations as MO
import BooleanDOI_processing as BDOIp
import sm_csm_operations as SMO



file_name='figure_12_Boolean_functions'



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


# calculate the mapping from original node names to index numbers
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
