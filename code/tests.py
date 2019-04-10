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
    #print("(2.4) of Toric paper")


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



G = get_24_toric_graph()
complexes = [
        [0, 4], #x0+x1
        [1, 4], #x2
        [2], #x3+1
        [0, 5], #x3 + x4
        [1, 5], #x5
        [3] #x0+x4
        ]

species = ['S0', 'S1', 'ES0', 'FS1', 'E', 'F']


partitions = [['ES0', 'FS1', 'S1P0', 'FP1'],
['S0', 'S1'],
['P0', 'P1'],
['E'],
['F']]


messi_network = MESSINetworkBuilder.MESSINetwork(G, complexes, species, partitions)

class Test1(unittest.TestCase):

    def assert_same_columns(self, np_array_1, np_array_2):
        array_1_columns = np_array_1.transpose().tolist()
        array_2_columns = np_array_2.transpose().tolist()

        for col_1 in array_1_columns:
            self.assertTrue(col_1 in array_2_columns)

        for col_2 in array_2_columns:
            self.assertTrue(col_2 in array_1_columns)


    def test_buildup_of_conformal_circuit_to_matrix(self):
        matrix_1 = np.array([
            [1, 0],
            [0, 3]])


        circuits_information_matrix_1 = CircuitsInformation(matrix_1)

        #Supposedly, this is the output of lemma A.5
        self.assertTrue([0,-3] in circuits_information_matrix_1.circuits)
        self.assertTrue([0,3] in circuits_information_matrix_1.circuits)
        self.assertTrue([3,0] in circuits_information_matrix_1.circuits)
        self.assertTrue([-3,0] in circuits_information_matrix_1.circuits)


        matrix_2 = np.array([
            [-1, 0, 2],
            [2, 4, 1]])


        circuits_information_matrix_2 = CircuitsInformation(matrix_2)

        for c in [[0, 4, 5], [0, -4, -5], [-4, 0, 8], [4, 0, -8], [-5, -8, 0], [5, 8, 0]]:
            self.assertTrue(c in circuits_information_matrix_2.circuits)




    def test_buildup_of_complexes_matrix_from_network(self):
        complexes_matrix = MESSIGraphUtils.build_complexes_matrix(messi_network)
        complexes_matrix_solution = numpy.array([
            [1, 0, 0, 0, 1, 0],
            [0, 1, 0, 0, 1, 0],
            [0, 0, 1, 0, 0, 0],
            [1, 0, 0, 0, 0, 1],
            [0, 1, 0, 0, 0, 1],
            [0, 0, 0, 1, 0, 0]]
            )

        self.assert_same_columns(complexes_matrix, complexes_matrix_solution)

    

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

        
        
        incidence_matrix = MESSIGraphUtils.build_incidence_matrix(messi_network)

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

        incidence_matrix_solution = numpy.array([
        [-1, 1,  0,  0,  0,  0],
        [0,  0,  1,  0,  0,  0],
        [1,  -1, -1, 0,  0,  0],
        [0,  0,   0, 0,  0,  1],
        [0,  0,   0, -1, 1,  0],
        [0,  0,   0, 1, -1, -1]])

        #The incidence matrix might have columns in another order.
        self.assert_same_columns(incidence_matrix, incidence_matrix_solution)


    def test_buildup_of_educt_complexes_matrix_from_network(self):
        #G = get_24_toric_graph()

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
        
        educt_complexes_matrix = MESSIGraphUtils.build_educt_complexes_matrix(messi_network)


        """
        ordered educts:
        s0+e, es0, es0, s1+f, fs1, fs1

        that is:
        x0+x4, x2, x2, x1+x5, x3, x3

        """

        educt_complexes_matrix_solution = numpy.array(
            [
            [1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 1, 1, 0, 0, 0],
            [0, 0, 0, 0, 1, 1],
            [1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0]
            ])

        self.assert_same_columns(educt_complexes_matrix, educt_complexes_matrix_solution)


    def test_buildup_of_stoichiometric_matrix_from_network(self):
        incidence_matrix = MESSIGraphUtils.build_incidence_matrix(messi_network)
        complexes_matrix = MESSIGraphUtils.build_complexes_matrix(messi_network)

        stoichiometric_matrix = MESSIGraphUtils.build_stoichiometric_matrix(incidence_matrix, complexes_matrix)


        #print("stoichiometric_matrix")
        #print(stoichiometric_matrix)

        #rango 3


        

        stoichiometric_matrix_solution = numpy.array(
            [[-1,  1,  0, 0,  1,  0],
            [ 0,  0,  1, -1,  0,  1],
            [ 1, -1, -1,  0,  0,  0],
            [ 0,  0,  0,  1, -1, -1],
            [-1,  1,  1,  0,  0,  0],
            [ 0,  0,  0, -1,  1,  1]]
            )


        self.assert_same_columns(stoichiometric_matrix, stoichiometric_matrix_solution)

    """def test_buildup_of_integer_basis_matrix_of_orthogonal_complement_of_stoichiometric_matrix(self):
        ortoghonal_complement_of_stoichiometric_matrix = MESSIGraphUtils.build_integer_basis_matrix_of_orthogonal_complement_of_stoichiometric_matrix(messi_network)

        stoichiometric_matrix = MESSIGraphUtils.build_stoichiometric_matrix_from_messi_network(messi_network)

        stoch_shape = np.shape(stoichiometric_matrix)
        complement_shape = np.shape(ortoghonal_complement_of_stoichiometric_matrix)

        #same columns
        self.assertTrue(stoch_shape[1] == complement_shape[1])


        print(stoch_shape)
        print(complement_shape)

        #dimensions are OK
        #self.assertTrue(stoch_shape[0] + complement_shape[0] == stoch_shape[1])

        stoch_rank = np.linalg.matrix_rank(stoichiometric_matrix)

        complement_rank =  np.linalg.matrix_rank(ortoghonal_complement_of_stoichiometric_matrix)


        print(stoch_rank)

        print(complement_rank)

        #No da! Hay algun otro vector en el espacio estequiometrico que no viene de un ciclo.
        self.assertTrue(stoch_shape[1] == stoch_rank + complement_rank)

        for row in complement_shape.rows():
            for column in stoichiometric_matrix.transpose():
                print(numpy.dot(row, column))



        print(ortoghonal_complement_of_stoichiometric_matrix)
        print("")
        print("")
        print("")
        print(stoichiometric_matrix)

        print(ortoghonal_complement_of_stoichiometric_matrix @ stoichiometric_matrix)
        print(stoichiometric_matrix @ ortoghonal_complement_of_stoichiometric_matrix)"""


    def test_get_integer_kernel(self):
        M = numpy.array([
            [1.0, 1, 1, 1],
            [2, 2, 2, 2]])

        K = MESSIGraphUtils.get_integer_kernel(M)
        for row in K:
            for a in row:
                self.assertTrue(a.is_integer())

        product = M @ K
        self.assertTrue(product.matrix_rank() ==0 )

if __name__ == '__main__':
    unittest.main()

"""
1) podriamos testear que los kappa sean positivos para algunos ejemplos

2) Ejemplos para testear generacion de matrices con nombre a partir de la estructura de grafo:



"""
