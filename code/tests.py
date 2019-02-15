import unittest

import numpy as np
from CircuitsInformation import CircuitsInformation
import Utils

#TODO: write an integration test (when - some day - we write the main output in a file)
class Test1(unittest.TestCase):


	def test__buildup_of_conformal_circuit_to_matrix_using_Lemma_A5(self):
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


if __name__ == '__main__':
    unittest.main()



"""
1) podriamos testear que los kappa sean positivos para algunos ejemplos

2) Ejemplos para testear generacion de matrices con nombre a partir de la estructura de grafo:



"""
