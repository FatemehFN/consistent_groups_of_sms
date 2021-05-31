import model_operations as MO
import BooleanDOI_processing as BDOIp




#file_name='figure_13_Boolean_functions'
file_name='figure_14_a_Boolean_functions'
#file_name='figure_14_b_Boolean_functions'


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



