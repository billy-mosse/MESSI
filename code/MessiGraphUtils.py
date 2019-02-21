
import numpy as np

class Complex:
	def __init__(self, species):
		"""species is a list. For example, ['S0', 'S1']
		"""
		self.species = species



def buildIncidenceMatrix(network):
    """
    the matrix whose i-th column has a 1 in
the row corresponding to the product complex of the i-th reaction and a −1 for the educt (reactant, i.e., origin)
complex.

    The network should be filled with renamed variables "0", "1", ...
    """

    #enumerate lets you loop through a structure (in this case, the edges)
    #while simultaneously keeping the index 0, 1, 2, ...
    #for index, edge in enumerate(nx.edges()):

    rows = []
    n_columns = len(network.nodes())

    #an empty row with amount of columns equal to amount of nodes of the network
    empty_row = [0]*n_columns

    for educt, product in network.edges():
        row = empty_row.copy()
        row[int(educt)] = -1
        row[int(product)] =1
        rows.append(row)


    return np.array(rows).transpose()

def buildComplexesMatrix(network):
    return None

def buildStochiometricMatrix(network):
    return None


class EductVector:
    def __init__(self, complexes):
        self.complexes = complexes

    def doSomeAlgebraStuff(self):
        print("hello")


def buildEductVector(network):
    return None


class MessiNetwork:

	#network es un networkx digraph
    def __init__(self, network):
        """
        intermediates and cores are both lists of (the class) complexes

        This function should build:
        1. Matrix N.
        2. The monomial vector phi(x)?
        3. I, D, Y.
        4. diag(k)? No.
        3. Matrix M. For that I need Mer's input.


        """
        self.incidence_matrix = buildIncidenceMatrix(network)
        self.complexes_matrix = buildComplexesMatrix(network)
        self.stochiometric_matrix = buildStochiometricMatrix(network)


        #what whould the structure of the  educt(origin vector) be?
        #I think it suffices for it to be a list of complexes with some hidden methods...
        self.educt_vector = buildEductVector(network)