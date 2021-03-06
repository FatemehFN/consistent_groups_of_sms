"""
This script shows the workflow for generating and saving the expanded network of the Boolean model in Figure 1 of the paper. First
the Boolean functions are read from a text file and the corresponding network Gread is built using the function form_network() from the
module BooleanDOI_processing. The function form_network() was written by Colin Campbell.
Then the network is fed into the function Get_expanded_network(), and the resulting expanded network is saved in G_expanded.
G_expanded is then given to write_g_expanded() to be written and saved in a gml file.
The user can also choose 'figure_13_Boolean_functions', 'figure_14_a_Boolean_functions', or 'figure_14_b_Boolean_functions' from the
Figures_data directory to get the expanded network of the Boolean models in Figure 13 or Figure 14 in the Appendix.
"""


import model_operations as MO
import BooleanDOI_processing as BDOIp




file_name='figure_1_Boolean_functions'

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



