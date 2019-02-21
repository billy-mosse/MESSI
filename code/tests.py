import unittest

import numpy as np
from CircuitsInformation import CircuitsInformation
import Utils

import MESSIGraphUtils
import networkx as nx
import numpy
import MESSINetworkBuilder

#TODO: write an integration test (when - some day - we write the main output in a file)


def get_24_toric_graph():
    print("(2.4) of Toric paper")


    """
    S0+E   --> y0
    S1+E   --> y1
    ES0    --> y2
    S0 + F --> y3
    S1 + F --> y4
    FS1    --> y5
    """
    reactions = [
    [0, 2,'k1'],
    [2, 0, 'k2'],
    [2, 1,'k3'],
    [4, 5, 'k4'],
    [5, 4, 'k5'],
    [5, 3, 'k6']
    ]        



    G = nx.DiGraph(directed=True)

    sources = set([reaction[0] for reaction in reactions])

    targets = set([reaction[1] for reaction in reactions])

    nodes = sources.union(targets)
    G.add_nodes_from(nodes)
    for reaction in reactions:
        G.add_edge(reaction[0], reaction[1], reaction_constant=reaction[2])
    return G


class Test1(unittest.TestCase):

    def assert_same_columns(self, np_array_1, np_array_2):
        array_1_columns = np_array_1.transpose().tolist()
        array_2_columns = np_array_2.transpose().tolist()

        for col_1 in array_1_columns:
            self.assertTrue(col_1 in array_2_columns)

        for col_2 in array_2_columns:
            self.assertTrue(col_2 in array_1_columns)


    def test_buildup_of_conformal_circuit_to_matrix_using_Lemma_A5(self):
        M1 = np.array([
            [1, 0],
            [0, 3]])


        circuits_information_M1 = CircuitsInformation(M1)

        #Supposedly, this is the output of lemma A.5
        self.assertTrue([0,-3] in circuits_information_M1.circuits)
        self.assertTrue([0,3] in circuits_information_M1.circuits)
        self.assertTrue([3,0] in circuits_information_M1.circuits)
        self.assertTrue([-3,0] in circuits_information_M1.circuits)


        M2 = np.array([
            [-1, 0, 2],
            [2, 4, 1]])


        circuits_information_M2 = CircuitsInformation(M2)

        for c in [[0, 4, 5], [0, -4, -5], [-4, 0, 8], [4, 0, -8], [-5, -8, 0], [5, 8, 0]]:
            self.assertTrue(c in circuits_information_M2.circuits)




    

    def test_buildup_of_incidence_matrix_from_network(self):
        

        """reactions = [
        ['S0+E', 'ES0','k1'],
        ['ES0', 'S0+E', 'k2'],
        ['ES0', 'S1+E','k3'],
        ['S1+F', 'FS1', 'k4'],
        ['FS1', 'S1+F', 'k5'],
        ['FS1', 'S0+F', 'k6']
        ]"""

        #The complexes are labeled y0, y1, y...

        G = get_24_toric_graph()
        messi_network = MESSINetworkBuilder.MESSINetwork(G, None)
        I = MESSIGraphUtils.build_incidence_matrix(messi_network)

        """I_solution = numpy.array([
        [-1, 1, 0, 0, 0, 0],
        [1, -1, -1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, -1, 1, 0],
        [0, 0, 0, 1, -1, -1],
        [0, 0, 0, 0, 0, 1]
        ])"""

        """

        reactions = [
    ['0', '2','k1'],
    ['2', '0', 'k2'],
    ['2', '1','k3'],
    ['4', '5', 'k4'],
    ['5', '4', 'k5'],
    ['5', '3', 'k6']
    ] """

        I_solution = numpy.array([
        [-1, 1,  0,  0,  0,  0],
        [0,  0,  1,  0,  0,  0],
        [1,  -1, -1, 0,  0,  0],
        [0,  0,   0, 0,  0,  1],
        [0,  0,   0, -1, 1,  0],
        [0,  0,   0, 1, -1, -1]])

        #The incidence matrix might have columns in another order.
        self.assert_same_columns(I, I_solution)


    def test_buildup_of_complexes_matrix_from_network(self):
        G = get_24_toric_graph()

        """
        S0+E   --> y0
        S1+E   --> y1
        ES0    --> y2
        S0 + F --> y3
        S1 + F --> y4
        FS1    --> y5
        """

        """
        S0  --> x0
        S1   --> x1
        ES0 --> x2
        FS1  --> x3
        E   --> x4
        F --> x5
        """
        complexes = [
        [0, 4], #x0+x1
        [1, 4], #x2
        [2], #x3+1
        [0, 5], #x3 + x4
        [1, 5], #x5
        [3] #x0+x4
        ]

        messi_network = MESSINetworkBuilder.MESSINetwork(G, complexes)
        complexes_matrix = MESSIGraphUtils.build_complexes_matrix(messi_network)


        """
        ordered educts:
        s0+e, es0, es0, s1+f, fs1, fs1

        that is:
        x0+x4, x2, x2, x1+x5, x3, x3

        """

        complexes_matrix_solution = numpy.array(
            [
            [1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 1, 1, 0, 0, 0],
            [0, 0, 0, 0, 1, 1],
            [1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0]
            ])

        self.assert_same_columns(complexes_matrix, complexes_matrix_solution)






if __name__ == '__main__':
    unittest.main()



"""
1) podriamos testear que los kappa sean positivos para algunos ejemplos

2) Ejemplos para testear generacion de matrices con nombre a partir de la estructura de grafo:



"""
