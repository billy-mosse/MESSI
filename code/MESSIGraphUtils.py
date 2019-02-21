
import numpy as np

class Complex:
    def __init__(self, species):
        """species is a list. For example, ['S0', 'S1']
        """
        self.species = species



def build_incidence_matrix(messi_network):
    """
    the matrix whose i-th column has a 1 in
the row corresponding to the product complex of the i-th reaction and a −1 for the educt (reactant, i.e., origin)
complex.

    The network should be filled with renamed variables "0", "1", ...
    """

    #enumerate lets you loop through a structure (in this case, the edges)
    #while simultaneously keeping the index 0, 1, 2, ...
    #for index, edge in enumerate(nx.edges()):

    G = messi_network.nx

    rows = []
    n_columns = len(G.nodes())

    #an empty row with amount of columns equal to amount of nodes of the network
    empty_row = [0]*n_columns

    #print(G.edges())
    for educt, product in G.edges():
        row = empty_row.copy()
        row[int(educt)] = -1
        row[int(product)] =1
        rows.append(row)


    return np.array(rows).transpose()

def build_complexes_matrix(messi_network):
    rows = []
    G = messi_network.nx
    n_columns = len(G.nodes())
    empty_row = [0]*n_columns

    for educt, _ in G.edges():
        row = empty_row.copy()

        for species in messi_network.complexes[educt]:
            row[species] = 1
        rows.append(row)

    return np.array(rows).transpose()


def build_stoichiometric_matrix(incidence_matrix, complexes_matrix):
    """
    returns the stoichiometrix matrix of the vector.
    If I is the incidence matrix, and Y is the complexes matrix,
    the stoichiometric matrix is defined as Y^t * I.
    """

    #@ is matrix multiplication
    return complexes_matrix.transpose() @ incidence_matrix

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