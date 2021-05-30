

import sm_csm_operations as SMO



plant=30
pollinator=25
network_number=266
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



