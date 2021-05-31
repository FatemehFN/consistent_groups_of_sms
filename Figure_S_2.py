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


