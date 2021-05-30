import networkx as nx
import model_operations as MO



#construction of the network example
G=nx.DiGraph()
node_list=['pl_1','pl_2','po_1','po_2','po_3','po_4','po_5']
G.add_nodes_from(node_list)
edge_list=[('pl_1','po_1'),('pl_1','po_2'),('pl_2','po_2'),('pl_2','po_3'),('pl_2','po_4'),('pl_2','po_5'),('po_2','pl_2'),
           ('po_1','pl_1'),('po_2','pl_1'),('po_3','pl_2'),('po_4','pl_2'),('po_5','pl_2')]

G.add_edges_from(edge_list)

edge_sign={}
for edge in G.edges()[0:6]:
    edge_sign.update({edge:1})


edge_sign.update({G.edges()[6]:-1})
edge_sign.update({G.edges()[7]:1})
for edge in G.edges()[8:]:
    edge_sign.update({edge:-1})


nx.set_edge_attributes(G,'data',edge_sign)



#the edges in the original network
print('The edges in the original network')
for item in G.edges(data=True):
    print(item)


#generates the Boolean functions in disjunctive prime form in a text file
MO.disjunctive_prime_form_text_file(G,network_number=1,plant=2,pollinator=5)


