"""
This scripts shows the workflow for generating the disjunctive prime form and the simplified Boolean functions of the threshold
functions of plant-pollinator interaction networks. First the network is read from file, and then it is fed into the function
disjunctive_prime_form_text_file() to generate the disjunctive prime form of the Boolean functions in a text file. It can also be fed into
the function simplification_text_file() to generate the simplified Boolean model written in a text file.

Author: Fatemeh Sadat Fatemi Nasrollahi unless otherwise noted.
Date: May 2021
Python Version: 3.7
"""


import pickle
import numpy as np
import model_operations as MO




plant=20
pollinator=10
network_number=987

f = open('graphlist'+str(plant)+'_'+str(pollinator)+'.data', 'rb')
graphlist = pickle.load(f, encoding='latin1')
f.close()

#initial conditions and perturbation parameters if needed
IC=np.ones((plant+pollinator,), dtype=int)
IC=list(IC)
IC_number=100
write_IC=False
more=False
less=True
wp=0.1
perturb=False

G = graphlist[network_number]

print('the edges in the network')
for edge in G.edges(data=True):
    print(edge)

#generates the Boolean functions in disjunctive prime form in a text file
MO.disjunctive_prime_form_text_file(G,network_number,plant,pollinator)

#example of adding initial conditions and perturbations to the text file
#MO.disjunctive_prime_form_text_file(G,network_number,plant,pollinator,write_IC=write_IC,
# request_for_p_n_edges=False,perturb=perturb,IC=IC,wp=wp,IC_number=IC_number)



#generates the simplified Boolean model according to the method described in the paper and writes the update functions in a text file
MO.simplification_text_file(G,network_number,plant,pollinator)

#example of adding initial conditions and perturbations to the text file
#MO.simplification_text_file(G,network_number,plant,pollinator,IC=IC,IC_number=IC_number,write_IC=write_IC,
# more=more,less=less,perturb=perturb,wp=wp)