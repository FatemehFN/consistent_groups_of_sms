import attrs_operations as AO
import model_operations as MO
import relationship_operations as RO
import FVS
import BooleanDOI_processing as BDOIp
import PyBoolNet


plant=35
pollinator=35
network_number=937
file_name=str(plant)+'_'+str(pollinator)+'_'+str(network_number)



#read the text file of the Boolean model
f=open('plant_pollinator_models/'+file_name+'_less.txt')
lines=f.readlines()
f.close()

f1=open('plant_pollinator_models/'+file_name+'_less.txt')
rules=f1.read()
f1.close()


#number of negative edges in the model
print('number of negative edges in the model')
print(MO.number_of_negative_edges_text(lines))




# number of minimal trap spaces from consistent groups of sms
Gread,readnodes = BDOIp.form_network(lines, sorted_nodename=False)
G_rel,list_of_names=RO.construct_relationships_network(lines,rules,write_cycle_graph=False)
FVS_size = len(FVS.FVS(Gread))
attrs,len_attrs=AO.consistent_groups_of_sms(G_rel,FVS_size,list_of_names)
print('number of minimal trap spaces from consistent groups of sms:')
print(len_attrs)







# number of minimal trap spaces from PyBoolNet
rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
    ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES','')
primes = PyBoolNet.FileExchange.bnet2primes(rules)
mints = PyBoolNet.AspSolver.trap_spaces(primes, "min", MaxOutput=2000)
print('number of minimal trap spaces from PyBoolNet:')
print(len(mints))