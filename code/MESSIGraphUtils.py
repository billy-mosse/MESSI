
import numpy as np

class Complex:
	def __init__(self, species):
		"""species is a list. For example, ['S0', 'S1']
		"""
		self.species = species



def build_incidence_matrix(messi_network):
    """
    this function builds the incidence matrix, that is,
    the matrix whose i-th column has a 1 in
the row corresponding to the product complex of the i-th reaction and a −1 for the educt (reactant, i.e., origin)
complex, as explained in section 5.1 of the toric paper.
    """

    #The network should be filled with renamed variables "0", "1", ...
    #enumerate lets you loop through a structure (in this case, the edges)
    #while simultaneously keeping the index 0, 1, 2, ...
    #for index, edge in enumerate(nx.edges()):

    G = messi_network.nx

    rows = []
    n_columns = len(G.nodes())

    #an empty row with amount of columns equal to amount of nodes of the network
    empty_row = [0]*n_columns

    for educt, product in G.edges():
        row = empty_row.copy()
        row[educt] = -1
        row[product] =1
        rows.append(row)

    return np.array(rows).transpose()

def build_complexes_matrix(messi_network):
    """
    this function builds the educt complex matrix,
    as explained in section 5.1 of the toric paper.
    """
    rows = []
    G = messi_network.nx
    n_columns = len(G.nodes())
    empty_row = [0]*n_columns

    #an underscore is used so that we can ignore the product variable - we won't be using it.
    #for every educt of edge of index i
    for educt, _ in G.edges():
        row = empty_row.copy()

        #we write a 1 in row i. It's easier to build the matrix by rows and then transpose
        for species in messi_network.complexes[educt]:
            row[species] = 1
        rows.append(row)

    #we returned the transposed matrix,
    #that is, the matrix that in column i has a 1 in row j
    #if species j belongs to the educt of edge i.
    return np.array(rows).transpose()


def build_stochiometric_matrix(network):
    return None

'''
#Creo que esto esta totalmente de mas
#VOLAR
class EductVector:
    def __init__(self, complexes):
        self.complexes = complexes

    def doSomeAlgebraStuff(self):
        print("hello")


def buildEductVector(network):
    return None


#Idem, volar
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
        self.educt_vector = buildEductVector(network)'''