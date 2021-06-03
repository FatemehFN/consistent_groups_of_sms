"""
This script shows the workflow for generating and writing the cycle graph for any imported Boolean model. First the model is read from
a text file, then it is fed into the function cycle_graph_virtual_node_based() to generate the cycle graph. if write_cycle_graph=True, the
cycle graph is written in a gml file.
"""


import sm_csm_operations as SMO




plant=20
pollinator=10
network_number=987
file_name=str(plant)+'_'+str(pollinator)+'_'+str(network_number)

#read the text file of the Boolean model
f=open('plant_pollinator_models/'+file_name+'_less.txt')
lines=f.readlines()
f.close()

f1=open('plant_pollinator_models/'+file_name+'_less.txt')
rules=f1.read()
f1.close()

#build cycle graph using virtual nodes of expanded network
cycle_graph=SMO.cycle_graph_virtual_node_based(lines,write_cycle_graph=True)



