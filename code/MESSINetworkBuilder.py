#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


"""
Python module for generating relevant data from a MESSI network.
TODOs:
 1. automatic checking of toricity.
 2. generation of binomials.
 3. graphs of multistationarity.
 4. nicer visualization of network?
"""


def get_relevant_matrices(debug):

    reactions = []
    if not debug:

        print("Welcome to an implementation of Algorithm 1 of the paper")
        print ("Please, input the graph G corresponding to the reaction network you want to analyze")
        print ("Input each reaction with the following format: (SOURCE, TARGET, REACTION_CONSTANT)")
        print ("One example could be: (S0 +E, ES0, k1)")
        print ("When finished, enter the command END")

    #Descomentar para Alicia

        r_input = input()
        while r_input != 'END':
            r = r_input[1:-1].split(',')
            r = [x.strip() for x in r]
            reactions.append(r)
            r_input = input()
    else:
        reactions = [['S0+E', 'ES0','k1'],
        ['ES0', 'S0+E', 'k2'],
        ['ES0', 'S1+E','k3'],
        ['S1+F', 'FS1', 'k4'],
        ['FS1', 'S1+F', 'k5'],
        ['FS1', 'S0+F', 'k6'],
        ['P0+S1','S1P0', 'k7'],
        ['S1P0', 'P0+S1', 'k8'],
        ['S1P0','P1+S1', 'k9'],
        ['P1+F','FP1', 'k10'],
        ['FP1','P1+F', 'k11'],
        ['FP1','P0+F', 'k12']]


    #This part needs some work.
    #I think we should use a multigraph for (label) visualization but a digraph for computations
    G = nx.DiGraph(directed=True)

    sources = set([reaction[0] for reaction in reactions])
    targets = set([reaction[1] for reaction in reactions])
    nodes = sources.union(targets)


    G.add_nodes_from(nodes)
    for reaction in reactions:
        G.add_edge(reaction[0], reaction[1], reaction_constant=reaction[2])

    df = pd.DataFrame(index=G.nodes(), columns=G.nodes())
    for row, data in nx.shortest_path_length(G):
        for col, dist in data.items():
            df.loc[row,col] = dist

    df = df.fillna(df.max().max())  

    layout = nx.kamada_kawai_layout(G, dist=df.to_dict())

    #pos = nx.spring_layout(G)

    nx.draw(G,layout,arrows=True, edge_color='black',width=1,linewidths=1,\
node_size=500,node_color='pink',alpha=0.9,arrowsize=12,arrowstyle='-|>',\
labels={node:node for node in G.nodes()},
)


    edge_labels = nx.get_edge_attributes(G,'reaction_constant')
    nx.draw_networkx_edge_labels(G, pos=layout, edge_labels = edge_labels)

    if not debug:
        plt.show()

    if not debug:
        answer = ""
        print("Does the network have a MESSI structure? Write YES or NO")
        while(answer!= "YES" and answer != "NO"):
            answer = input()
    else:
        answer = "YES"
        
    if answer == "YES":

        P0_intermediates = []
        P_cores = []
        
        if not debug:
            print("Great! Then if it's s-toric we will be able to check multistationarity.")
            print("Please introduce the nodes for P^0. As usual, write END when finished:")
            node = ""
            while(node != "END"):
                node = input()
                if node != "END":
                    P0_intermediates.append(node)
            
            i = 1
            #TODO: falta chequear que no sean vacios y hay al menos 1.
            node = ""
            while(node != "END_CORES"):
                print("Please introduce the nodes for P^%d. Press END to advance to the next set and END_CORES to let us know this was the last set." % i)

                node = ""
                i = i+1
                P = []
                while(node != "END" and node != "END_CORES"):
                    node = input()
                    if node != "END" and node != "END_CORES":
                        P.append(node)

                P_cores.append(P)
        else:
            P0_intermediates=['ES0','FS1','S1P0','FP1']
            P_cores=[['S0','S1'],['P0','P1'],['E'],['F']]

        #TODO: check it.
        print("I'm trusting you that this is a MESSI structure...")

        print("Let's check that it is s-toric...")

        cores = [y for x in P_cores for y in x]

        is_s_toric = True

        print("Checking C'...")
        for intermediate in P0_intermediates:
            print("checking intermediate complex %s..." % intermediate)


            simple_o_paths_from_core_source = 0
            for source in sources:

                #Deberia bastar solo los simple paths. Chequear
                simple_paths = nx.all_simple_paths(G, source, intermediate)

                L = list(simple_paths)

                for simple_path in L:
                    only_goes_through_intermediates=True
                    source = simple_path[0]

                    if source not in P0_intermediates:
                        for index, c in enumerate(simple_path):
                            if index > 0:
                                if c not in P0_intermediates:
                                    only_goes_through_intermediates=False
                                    break;
                        if only_goes_through_intermediates:
                            print("Path %s only goes through intermediates" % simple_path)

                            simple_o_paths_from_core_source = simple_o_paths_from_core_source +1


            if simple_o_paths_from_core_source!=1:
                is_s_toric=False
                break

        if is_s_toric:
            print("Condition C' holds!")
        else:
            print("Condition C' doesn't hold, so the network is not toric")

        print("")
        print("")
        print("TODO: check other conditions.")
        print("Assuming the network is s-toric, by Theorem 4.8, it's toric.")
        print("Morover, the exponents of the binomials can be calculated explicitely.")
        print("to be continued...")

        if not debug:
            var = input("Press ENTER to continue with the program.")

        #TODO: finish.

'''
Para el primer mensaje:

(S0+E, ES0,k1)
(ES0, S0+E, k2)
(ES0, S1+E,k3)
(S1+F, FS1, k4)
(FS1, S1+F, k5)
(FS1, S0+F, k6)
(P0+S1,S1P0, k7)
(S1P0, P0+S1, k8)
(S1P0,P1+S1, k9)
(P1+F,FP1, k10)
(FP1,P1+F, k11)
(FP1,P0+F, k12)
END

'''

#P0_intermediates=[ES0,FS1,S1P0,FP1]
#P_cores=[[S0,S1],[P0,P1],[E],[F]]
'''

Para el segundo mensaje:

ES0
FS1
S1P0
FP1
END
S0
S1
END
P0
P1
END
E
END
F
END_CORES
'''