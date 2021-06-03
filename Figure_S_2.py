"""
This script shows the workflow for generating the network of functional relationships in Figure 5. The Boolean model is
read from the text file and given to function construct_relationships_network() to construct the network of functional relationships among
stable motifs and motif groups. One can see the nodes and edges of this DiGraph in the output.
"""


import relationship_operations as RO


file_name='figure_5_Boolean_functions'

#read the text file of the Boolean model
f=open('Figures_data/'+file_name+'.txt')
lines=f.readlines()
f.close()


f1=open('Figures_data/'+file_name+'.txt')
rules=f1.read()
f1.close()


G_rel,list_of_names=RO.construct_relationships_network(lines,rules,write_cycle_graph=False)

print('\nnodes in network of functional relationships')
for node in G_rel.nodes(data=True):
    print(node)

print('\nedges in network of functional relationships')
for edge in G_rel.edges(data=True):
    print(edge)


