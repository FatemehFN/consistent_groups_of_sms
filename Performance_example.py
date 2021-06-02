import PyBoolNet
import StableMotifs as sm
import relationship_operations as RO
import attrs_operations as AO
import time



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




# number of minimal trap spaces from consistent groups of sms
start_cg=time.time()
G_rel,list_of_names=RO.construct_relationships_network(lines,rules,write_cycle_graph=False)
attrs,len_attrs=AO.consistent_groups_of_sms(G_rel,list_of_names,lines)
print('number of minimal trap spaces from consistent groups of sms:')
print(len_attrs)
end_cg=time.time()




# number of minimal trap spaces from PyBoolNet
start_mints=time.time()
rules = rules.replace(' *=', ',\t').replace('*=', ',\t').replace('not ', '!').replace(' and ', ' & ').replace(
    ' or ', ' | ').replace('False', '0').replace('#BOOLEAN RULES','')
primes = PyBoolNet.FileExchange.bnet2primes(rules)
mints = PyBoolNet.AspSolver.trap_spaces(primes, "min", MaxOutput=2000)
print('number of minimal trap spaces from PyBoolNet:')
print(len(mints))
end_mints=time.time()




# number of attractors from Stable Motifs 2021
start_21=time.time()
sm.Format.pretty_print_prime_rules({k: primes[k] for k in sorted(primes)})
max_simulate_size = 0
ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
sum_attrs=0
for a in ar.attractors:
    sum_attrs+=1
print('number of minimal trap spaces from Stable Motifs 2021:')
print(sum_attrs)
end_21=time.time()



print('performance summary:')
print('Consistent groups of stable motifs finished the task in %s seconds' % str(end_cg-start_cg))
print('PyBoolNet finished the task in %s seconds' % str(end_mints-start_mints))
print('StableMotifs2021 finished the task in %s seconds' % str(end_21-start_21))

