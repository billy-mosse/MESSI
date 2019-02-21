import unittest

import numpy as np
from CircuitsInformation import CircuitsInformation
import Utils

import MessiGraphUtils
import networkx as nx
import numpy

#TODO: write an integration test (when - some day - we write the main output in a file)
class Test1(unittest.TestCase):


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
        
        print("(2.4) of Toric paper")

        """reactions = [
        ['S0+E', 'ES0','k1'],
        ['ES0', 'S0+E', 'k2'],
        ['ES0', 'S1+E','k3'],
        ['S1+F', 'FS1', 'k4'],
        ['FS1', 'S1+F', 'k5'],
        ['FS1', 'S0+F', 'k6']
        ]"""

        reactions = [
        ['0', '1','k1'],
        ['1', '0', 'k2'],
        ['1', '2','k3'],
        ['3', '4', 'k4'],
        ['4', '3', 'k5'],
        ['4', '5', 'k6']
        ]        

        G = nx.DiGraph(directed=True)

        sources = set([reaction[0] for reaction in reactions])

        targets = set([reaction[1] for reaction in reactions])

        nodes = sources.union(targets)
        G.add_nodes_from(nodes)
        for reaction in reactions:
            G.add_edge(reaction[0], reaction[1], reaction_constant=reaction[2])

        I = MessiGraphUtils.buildIncidenceMatrix(G)

        I_solution = numpy.array([
        [-1, 1, 0, 0, 0, 0],
        [1, -1, -1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, -1, 1, 0],
        [0, 0, 0, 1, -1, -1],
        [0, 0, 0, 0, 0, 1]
        ])

        #The incidence matrix might have columns in another order.
        I_columns = I.transpose().tolist()
        I_solution_columns = I_solution.transpose().tolist()

        for I_col in I_columns:
            self.assertTrue(I_col in I_solution_columns)

        for I_solution_col in I_solution_columns:
            self.assertTrue(I_solution_col in I_columns)

        #v = numpy.array([[1,2,3],[4,5,6]])
        #self.assertTrue([1,2,3] in v)

if __name__ == '__main__':
    unittest.main()



"""
1) podriamos testear que los kappa sean positivos para algunos ejemplos

2) Ejemplos para testear generacion de matrices con nombre a partir de la estructura de grafo:



"""
