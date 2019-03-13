
import numpy as np
import networkx as nx

class Complex:
    def __init__(self, species):
        """species is a list. For example, ['S0', 'S1']
        """
        self.species = species


#TODO: estos van a ser metodos de la clase, me parece.
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

    G = messi_network.G

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

def build_educt_complexes_matrix(messi_network):
    rows = []
    G = messi_network.G
    n_columns = len(G.nodes())
    empty_row = [0]*n_columns

    #p#rint("Edges: %s" % G.edges())
    for educt, _ in G.edges():
        row = empty_row.copy()

        for species in messi_network.complexes[educt]:
            row[species] = 1
        rows.append(row)

    return np.array(rows).transpose()


def build_complexes_matrix(messi_network):
    #Y^t tiene, en cada columna, al complejo
    #Entonces Y tiene, en cada fila, al complejo
    #Cada fila tiene de tamaño la cantidad de especies
    n_columns = len(messi_network.species)
    empty_row = [0]*n_columns
    rows = []
    for complex in messi_network.complexes:
        row = empty_row.copy()
        for species_index in complex:
            row[species_index] = 1
        rows.append(row)
    return np.array(rows)



def build_stoichiometric_matrix(incidence_matrix, complexes_matrix):
    """
    returns the stoichiometrix matrix of the vector.
    If I is the incidence matrix, and Y is the complexes matrix,
    the stoichiometric matrix is defined as Y^t * I.
    """

    #@ is matrix multiplication

    print("complexes matrix")

    return complexes_matrix.transpose() @ incidence_matrix


def build_stoichiometric_matrix_from_messi_network(messi_network):
    incidence_matrix = build_incidence_matrix(messi_network)
    complexes_matrix = build_complexes_matrix(messi_network)

    return build_stoichiometric_matrix(incidence_matrix, complexes_matrix)

def build_integer_basis_matrix_of_orthogonal_complement_of_stoichiometric_matrix(messi_network):
    """
    returns an integer basis of the orthogonal complement of the stoichiometric matrix of a graph
    
    """
    """G = messi_network.G
    n_columns = len(G.nodes())
    empty_row = [0]*n_columns

    simple_cycles = nx.simple_cycles(G)
    rows =[]
    for cycle in simple_cycles:
        row=empty_row.copy()
        #print("cicle: %s" % cycle)
        for complex_index in cycle:
            complex = messi_network.complexes[complex_index]
            print("complex: %s" % complex)
            for species in complex:                
                row[species] +=1
        #print(row)
        rows.append(row)

    #print("rows: %s" % rows)
    return np.array(rows)"""



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