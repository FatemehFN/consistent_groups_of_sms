# consistent_groups_of_sms

This repository is the python implementation of the attractor identification method described in the follwoing work:

Fatemeh Sadat Fatemi Nasrollahi, Jorge Gómez Tejeda Zañudo, Colin Campbell, and Réka Albert. The relationships among generalized positive feedback loops determine possible community outcomes in plant-pollinator interaction networks.

In this work we develop an attractor identification method in plant-pollinator interaction networks governed by threshold functions (see [Campbell et al. 2011](https://www.pnas.org/content/108/1/197.short) for more details on the model). 

## Features

This repository has a set of modules that are useful for:
* Generating the Boolean disjunctive prime form from the network of the original species pool and the threshold functions. This is currently hard coded for the threshold functions of Campbell et al. model (with the weight of 4 for positive regulators and -1 for negative regulators).
* Generating the simplified Boolean model from the network of the original species pool and the threshold functions. This is currently hard coded for the threshold functions of Campbell et al. model (with the weight of 4 for positive regulators and -1 for negative regulators).
* Build and save the expanded network for the Boolean models.
* Find the stable generalized positive feedback loops (stable motifs) of Boolean models.
* Construct and save the cycle graph in which nodes are the consistent conditional cycles of the expanded network and edges represent partial condition satisfaction between the nodes. For example an edge from cycle C<sub>1</sub> to cycle C<sub>2</sub> means that a virtual node in C<sub>1</sub> satisfies a condition in C<sub>2</sub>. 
* Find the conditonally stable motifs and their supports in plant-pollinator Boolean models. 
* Construct the network of functional relationships according to the identified stable motifs, conditionally stable motifs and their supports. In this network nodes are either stable motifs or motif groups (each motif group consists of one conditionally stable motif and one of its supports), and edges represent logical determination and mutual exclusivity between the nodes. Logical determination edge from node M<sub>1</sub> to node M<sub>2</sub> shows that if M<sub>1</sub> stabilizes, M<sub>2</sub> automatically stabilizes. Mutual exclusivity between nodes M<sub>1</sub> and M<sub>2</sub> shows that they both contain the same node, but in different states. 
* Find maximal consistent groups of stable motifs and motif groups that mutually exclude each other. These groups lead to attractors of the Boolean model. 

## <a name="attr"></a> Attractor identification workflow

First the simplified Boolean model is read from a text file. These Boolean functions should be in booleannet format such as:
```
Node1*=Node2 or Node3
Node2*=Node1 and not Node3
Node3*=Node2
```
The text before *=* is the node name, the * shows the future state of the node, and the text after *=* is the Boolean regulatory function of the node. 
Note that:
1. The name of one node should not be a part of the name of another node.
2. The functions do not need to be in disjunctive prime form. The code has functions that perform such conversions. 
3. For the node names use only alphanumeric characters (A-Z,a-z), numbers (0-9) and *_*. The reserved words for the program, which shouldn't be used for node names, are: *True*, *False*, *true*, *false*, *0*, *1*, *and*, *or*, and *not*.  

The model then is analyzed to find the stable motifs and the conditionally stable motifs in plant-pollinator Boolean models in the work of Campbell et al. 2011. Then it finds three relationships between such stable motifs: dependence, mutual exclusivity, and logical determination. Based on these relationships, it finds self-consistent groups of these motifs that mutually exclude each other. These groups lead to attractors of the Boolean models. In plant-pollinator interaction networks the attractors correspond to stable community outcomes. 

## Operation files

There are five modules each of which focuses on a group of operations that are similar:

1. **Model operations (model_operations.py):** It contains the functions that generate the Boolean functions of the plant-pollinator models of Campbell et al. work from the original networks. Each original network is a DiGraph object consisting of all the species in the original species pool as nodes and the interactions among them as edges. The data of each edge consists of its sign and determines whether the interaction is beneficial for the receiving end or not.  
There are two main functions: 1. `disjunctive_prime_form_text_file()` that receives the original network in DiGraph format, and generates the Boolean functions in disjunctive prime form based on the threshold model of Campbell et al.. It saves the Boolean functions in a text file. 2. `simplification_text_file()` that receives the original network in DiGraph format and generates Boolean functions of the simplified model based on the threshold model of Campbell et al.. Simplification works based on 1. sampling negative edges instead of keeping all of them, 2. having the same probabality of target node being active after simplification. It saves the Boolean functions in a text file.
3. **Name operations (name_operations.py):** It contains the functions that converts different notations and formats to each other. 
4. **Stable motif and conditionally stable motif operations (sm_csm_operations.py):** It contains the functions that build cycle graph, find stable motifs, conditionally stable motifs, and the support of conditionally stable motifs. 
5. **Relationship operarions (relationship_operations.py):** It contains the function that construct the network of functional relationships in which nodes are stable motifs or motif groups (each motif group is a conditionally stable motif with one of its supports) and edges are functional relationships logical determination and mutual exclusivity between the nodes. It also contains all other functions that find and detect logical determinations and mutual exclusivities among the nodes. 
6. **Attractor operations (attrs_operations.py):** It contains the functions that find the self consistent groups of stable motifs and motif groups that are mutually exclusive. These groups lead to distinct attractors. 

## Data files

There are model examples of plant-pollinator and plant-pollinator like networks in two directories:

1. **plant_pollinator_models** is a directory that has several plant-pollinator models stroed in text files. The name of each text file consistes of thress number separated with an undeline. These numbers show the number of plants in that network, the number of pollinators in that network, and the network ID in the ensemble respectively. Each file contains the Boolean functions of the simplified model. The species (nodes) that do not have positive regulators cannot establish and their regulatory function is equal to False, the species that only have positive regulators after simplification have a regulatory function with positive regulators connected with *OR* operator, and the species that have both positive and negative regulators after simplification have a regulatory function with positive regulators grouped together with *OR* operator in parantheses and each of the negative regulators connected by *AND NOT*. The nodes are listed in the order they were saved in the original network in DiGraph format.  
2. **other_models** is a directory that has network models with Boolean functions that follow the regularities in plant-pollinator models of Cambell et al.. They do not have self loops and negative regulation is always accomponied by positive regulation in the opposite side. The positive regulators are connected with *OR* operator and negative regulators are connected with *AND NOT*.

The Boolean regulatory functions of models used in the figures of the paper are stored in text file format in **Figures_data** directory.

## Example implementations

There are several python files with examples. These files demonstrate the workflow in simple use cases.
1. File example_1.py demostrates generation of disjunctive prime form and simplified Boolean models for plant pollinator networks in text files. There are several inputs that can be provided if necessary: perturb, wp, IC, and write_IC. These inputs provide the option of having initial conditions or/and perturbation in the Boolean model written in the text files. Functions `disjunctive_prime_form_text_file()` and `simplification_text_file()` write the Boolean functions in text files. 
2. File example_2.py demonstrates how one can use the funtion `cycle_graph_virtual_node_based()` to construct the cycle graph. The Boolean functions in booleannet format explained [above](#attr) should be read from a text file using fundtion `readlines()`. This is the only argument to the function `cycle_graph_virtual_node_based()`. It writes the cycle graph in a gml file that can be opened by YED. 
3. Files example_3.py has an example of code execution for a plant pollinator network model. First, the Boolean functions are read from a text file. The functions should be in booleannet format as noted [above](#attr). Then the network of functional relationships should be constructed using the function `construct_relationships_network()`. Then that network should be provided for the `consistent_groups_of_sms()` function as an argument. This function finds the consistent groups of stable motifs and conditionally stable motifs that are mutually exclusive. It returns the number of attractors and the list of these mutually exclusive groups. The user can also choose model_1.txt or model_2.txt from the 'other_models' directory to get the same result for a non plant-pollinator network.  


## Main analysis

### Figures in the paper
These files are used to perform the analysis reported in the paper:
* Figure_1.py illustrates the generation of the expanded network for the Boolean functions of figure 1 in the paper. 
* Figure_2.py illustrates how the conversion of threshold functions to disjunctive prime form causes edge loss for the model in figure 2. The edges in the original network can be seen in the output, and the edges in the disjunctive prime form can be seen in the text file generated. 
* Figure_4_5.py illustrates the generation of expanded network and cycle graph for the Boolean functions of the models in figure 4 and figure 5. Note that since we use cycle graph to identify conditionally stable motifs, we do not add the unconditional consistent cycles to it. 
* Figure_7.py implements the function `construct_relationships_network()` to generate the network of functional relationships for the Boolean functions in figure 7 in the paper. The edges and the nodes of this network are then printed in the output. 

### Figures in the appendices
* Figure_11.py illustrates the generation of the data in figure 11 in the appendix. It solves inequality B1 for each number of positive and negative regulators. 
* Figure_12.py illustrates the generation of expanded network and cycle graph for the Boolean functions of the models in figure 12. The code is the same as Figure_4_5.py, but it runs for a different model. 
* Figure_13_14.py illustrates the generation of the expanded network for the Boolean functions of figure 13 and figure 14 in the appendix. The code is the same as Figure_1.py, but it runs for a different model. 

### Figures in supplemental material
* Figure_S_2.py generates the network of functional relationships in figure S2 using the function `construct_relationships_network()`. 

### Performance analysis

* Performance_example.py illustrates the use of attractor identification functions of 'SM analysis 2021' and 'PyBoolNet'. This script is the example of the performance analysis we did for one of the networks. Note that 'SM analysis 2013' is written in java, so one would either need to run the java code independently or import subprocess to be able to run it from python. 

## Requirements

* Networkx 1.11
* PyBoolNet 2.3+
* NumPy 1.20+
* Pandas 1.1.5+
* Sympy 1.5.1+
* StableMotifs

## Other information

### Notations
---
The goal of this section is to provide information about different notations that are used in different steps and functions.
The function `form_network()` from the module BooleanDOI_processing reads the Boolean functions from a text file and constructs the network accordingly. The network is stored in a networkx DiGraph object. This function relabels the nodes and assigns a new name in the format of *n*+ a number to each node.

```
Example of mapping from original node names to number index notation:

pl_1*=po_1
po_1*=pl_1

original notation          number index
        pl_1                    n1
        po_1                    n2
```


* The network that is constructed from the text file of the Boolean functions is in 'number index' notation.  
* The expanded network is a network that represents both the structure and the logical functions. For each node in the interaction network, there are two virtual nodes in the expanded network each representing one of the states the node can take: The active state is shown with the name of the node, and the inactive state is shown with *~*+ the name of the node. This network is written in 'number index' notation.  
* The LDOI operations are performed on the exapnded network and hence are in 'number index' notation.  
* The stable motifs generated by the PyBoolNet function `trap_spaces()` is in the *pl-po* notation. Example: `[{pl_0:0, po_1:0},{pl_0:1, po_1:1}]`. Each dictionary represents one stable motif and consists of virtual nodes, i.e., the nodes as the keys and their state as the values in the dictionary. 
* The name of each node in the network of functional relationships is a string consisting of the virtual nodes, i.e., the nodes and their states in that specific stable motif/motif group. They are written in *pl-po* notation.  
* 'mapping' provides this information from the node names (original notation) to number index notation. 
reverse mapping provides the mapping from the number index notation to node names (original notation). 
* Consistent cycles and conditionally stable motifs are stored in this format: `[{'n1': 0, 'n2': 0}, {'~n3', '~n0'}]`.  
The first item is a dictionary of virtual nodes, i.e., the nodes and their states in that consistent cycle/consitionally stable motif, and the second item is a set of conditions with '~' showing the inactive state. It is stored in this format so that it would be easier to check if this set is a subset of LDOI of a set of node states. 

### softwares and codes used 
---
#### Ripositories
* We used the python package PyBoolNet developed in [Klarner et al. 2017](https://academic.oup.com/bioinformatics/article/33/5/770/2725550?login=true) that can be found [here](https://github.com/hklarner/PyBoolNet) to find stable motifs and minimal trap spaces. The latter was found using PyBoolNet as a reference.
* For performance comparison analyses we used the java implementation of [SM analysis 2013](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004193) that can be found [here](https://github.com/jgtz/StableMotifs), and python implementation of [SM analysis 2021](https://arxiv.org/abs/2009.05526) that can be found [here](https://github.com/jcrozum/StableMotifs). 

#### Files
* Files FVS.py and FVS_localsearch_10_python.py are taken from the [FVS python ripository](https://github.com/yanggangthu/FVS_python) by Gang Yang. 
* Files BooleanDOI_DOI.py, BooleanDOI_processing.py and qm.py are taken from the [Boolean DOI python ripository](https://github.com/yanggangthu/BooleanDOI) by Gang Yang. 
* The file cool_bool_tools.py is written by Dávid Deritei.




