
import numpy as np
import networkx as nx
from utils import CircuitUtils, Utils

def isString(word):
    return type(word) == str


def indexOfSpecies(species_name, messi_network):

    #TODO reescribir los tests para evitar este hack
    if not isString(species_name):
        return species_name

    return messi_network.species_names.index(species_name)

def indexOfComplex(complex_name, messi_network):

    #TODO reescribir los tests para evitar este hack
    if not isString(complex_name):
        return complex_name

    return messi_network.complexes_names.index(complex_name.strip())

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
    n_columns = len(messi_network.complexes)

    #an empty row with amount of columns equal to amount of nodes of the network
    empty_row = [0]*n_columns

    #print(G.edges())

    for educt, product in G.edges():
        row = empty_row.copy()

        row[int(indexOfComplex(educt, messi_network))] = -1
        row[int(indexOfComplex(product, messi_network))] =1
        rows.append(row)

    return np.array(rows).transpose()


#Esta función también parece estar bien
def build_educt_complexes_matrix(messi_network):
    rows = []
    G = messi_network.G
    n_columns = len(messi_network.species)
    empty_row = [0]*n_columns

    #p#rint("Edges: %s" % G.edges())
    for educt, _ in G.edges():
        row = empty_row.copy()

        for species in messi_network.complexes[educt]:
            row[species] = 1
        rows.append(row)

    return np.array(rows).transpose()


def build_complexes_matrix(messi_network):
    """Y^t tiene, en cada columna, al complejo
    Entonces Y tiene, en cada fila, al complejo
    Cada fila tiene de tamaño la cantidad de especies"""
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

    #print("complexes matrix")

    #This is matrix N in toric paper
    return complexes_matrix.transpose() @ incidence_matrix


def build_stoichiometric_matrix_from_messi_network(messi_network):
    incidence_matrix = build_incidence_matrix(messi_network)

    complexes_matrix = build_complexes_matrix(messi_network)

    return build_stoichiometric_matrix(incidence_matrix, complexes_matrix)



#Creo que a esta no la tengo que llamar con la stoichiometric matrix....
def get_basis_of_columns_and_index(long_matrix):
    nrows, ncolumns = np.shape(long_matrix)
    column_basis = []
    index_column_basis = []
    for index, column in enumerate(long_matrix.transpose()):
        tmp_basis = column_basis.copy()
        tmp_basis.append(column)

        if np.linalg.matrix_rank(tmp_basis) > np.linalg.matrix_rank(column_basis):
            column_basis.append(column)
            index_column_basis.append(index)

    return np.array(column_basis).transpose(), index_column_basis



def extract_column_basis(matrix):
    column_basis = []
    for column in matrix.transpose():
        tmp_basis = column_basis.copy()
        tmp_basis.append(column)

        if np.linalg.matrix_rank(tmp_basis) > np.linalg.matrix_rank(column_basis):
            column_basis.append(column)
    return np.array(column_basis)

def extract_row_basis(matrix):
    row_basis = []
    for row in matrix:
        tmp_basis = row_basis.copy()
        tmp_basis.append(row)
        if np.linalg.matrix_rank(tmp_basis) > np.linalg.matrix_rank(row_basis):
            row_basis.append(row)
    return np.array(row_basis)

def get_positive_matrix(matrix):
    #print("Matrix: ")
    #print(matrix)
    transposed_matrix = matrix.transpose()
    positive_vector = CircuitUtils.get_positive_vector_in_matrix_rows(transposed_matrix)
    #print("positive_vector")
    #print(positive_vector)

    for i in range(0, np.shape(transposed_matrix)[0]):
        while not Utils.is_nonnegative(transposed_matrix[i]):
            transposed_matrix[i] = np.add(transposed_matrix[i], positive_vector)


    return transposed_matrix.transpose()

def build_integer_basis_of_stoichiometric_matrix(messi_network):
    stoichiometric_matrix = build_stoichiometric_matrix_from_messi_network(messi_network)

    stoichiometric_matrix_column_basis = extract_column_basis(stoichiometric_matrix)
    return stoichiometric_matrix_column_basis


#Construye la matriz M de toric, que es el nucleo de la martiz N.
#
def build_positive_integer_basis_of_kernel_of_stoichiometric_matrix(messi_network):
    #en sus columnas tiene a las reacciones
    
    #Es de SxR, asi que deberia tener R columnas, o sea, 12.
    stoichiometric_matrix = build_stoichiometric_matrix_from_messi_network(messi_network)
    #print("stoichiometric_matrix")
    #print(stoichiometric_matrix)

    #Esto, lamentablemente, esta bien segun el paper de toric...
    row_basis = extract_row_basis(stoichiometric_matrix)
    #print("row basis")
    #print(row_basis)

    return build_integer_basis_matrix_of_orthogonal_complement_of_matrix(row_basis, positive = True)


def build_positive_integer_basis_of_ortogonal_complement_of_stoichiometric_matrix(messi_network):
    integer_basis_of_orthogonal_complement_of_stoichiometric_matrix = build_integer_basis_of_orthogonal_complement_of_stoichiometric_matrix(messi_network)
    return get_positive_matrix(integer_basis_of_orthogonal_complement_of_stoichiometric_matrix.transpose())

#TODO test what happens if the basis isnt in the first columns
#I think that situation ever arises, though
#TODO also positive
def build_integer_basis_matrix_of_orthogonal_complement_of_matrix(matrix, positive = False):
    column_basis, index_column_basis = get_basis_of_columns_and_index(matrix)

    column_basis_det = np.linalg.det(column_basis)
    inverse_of_columns = np.linalg.inv(column_basis)
    
    #Tiene la identidad en las columnas index_column_basis, y M' en el resto
    IdAndMPrime = inverse_of_columns @ matrix

    d, r = np.shape(IdAndMPrime)

    index_column_basis_complement = [i for i in range(0, r) if i not in index_column_basis]
    #Devuelvo -M' en las filas buenas y la identidad en las malas

    result = [0]* r
    empty_row = [0]* (r-d)

    only_Mp = []
    for i in index_column_basis_complement:
        only_Mp.append(IdAndMPrime.transpose()[i])

    only_Mp = np.array(only_Mp).transpose()

    for index, column_index in enumerate(index_column_basis):
        new_row = np.dot(-1, only_Mp[index])
        result[column_index] = new_row
    
    for index, column_index in enumerate(index_column_basis_complement):
            new_row = empty_row.copy()
            new_row[index] = 1
            result[column_index] = new_row


    integer_matrix =  np.dot(column_basis_det, np.array(result))
    if positive:
        return get_positive_matrix(integer_matrix)

    return integer_matrix


def build_integer_basis_of_orthogonal_complement_of_stoichiometric_matrix(messi_network):
    stoichiometric_matrix = build_stoichiometric_matrix_from_messi_network(messi_network)
    
    stoichiometric_matrix_column_basis = extract_column_basis(stoichiometric_matrix)
    
    #print(stoichiometric_matrix_column_basis)

    #HACK: transponemos porque stoichiometric_matrix_column_basis esta transpuesto
    ret = build_integer_basis_matrix_of_orthogonal_complement_of_matrix(stoichiometric_matrix_column_basis, positive=False).transpose()
    #M = stoichiometric_matrix_column_basis
    #M_ort = ret.transpose()
    #print(M @ M_ort)
    return ret

def get_unique_core_reacting_through_intermediates(messi_network, intermediate_index):

    intermediate_name = messi_network.partitions[0][intermediate_index]

    #It's monomolecular
    intermediate_complex = [intermediate_name]
    intermediate_complex_index = messi_network.complexes.index(intermediate_complex)

    #TODO: probablemente el index no es la cosa correcta que tengo que agarrar
    possible_core = next(messi_network.G.predecessors(intermediate_complex_index))
    found_complex = False
    while not found_complex:

        #Iterate predecessors until find a core
        #At bifurcations, we just choose the first one
        #because we are assuming that C' holds.

        #TODO: probablemente el indice esta mal.
        #pero la funcion predecessors() efectivamente existe.

        possible_core = next(messi_network.G.predecessors(possible_core))
        if possible_core in messi_network.core_complexes():
            found_complex = True

    return possible_core



def phi(messi_network, intermediate_index):
    unique_core_reacting_through_intermediates \
    = get_unique_core_reacting_through_intermediates(messi_network, intermediate_index)

    #This should be the indices of the complex,
    #seen as a vector of species    
    return messi_network.complexes[unique_core_reacting_through_intermediates]

def get_pairs_of_binomial_exponents_of_type_1(messi_network):
    L = []

    #p es la cantidad de intermedios
    for intermediate_index, intermediate in enumerate(messi_network.intermediates()):
        #print("intermediate: %s" % intermediate)
        item = [[messi_network.species_names.index(intermediate)], phi(messi_network, intermediate_index)]
        #print("Item")
        #print(item)
        L.append(item)

    #print("L")
    #print(L)
    return L

#TODO LLAMAR A ESTA FUNCION
def check_sufficient_condition_for_s_toricity(messi_network):
    condition = True
    for connected_component in messi_network.G2.connected_components():
        for [origin, target] in connected_component.edges():

            #hay un unico simple path desde target a origin?
            if len(nx.all_simple_paths(messi_network.G2_circle, origin, target)) != 1:
                condition = False

    return condition

def get_label_from_edge(G, edge, label):
    #print("Searching for: ")
    #print(edge)
    #print("___")
    for origin, target, data in G.edges(data=True):
        #print(origin)
        #print(target)
        #print(data)
        if origin == edge[0] and target == edge[1]:
            return data[label]

def get_pairs_of_binomial_exponents_of_type_2(messi_network):
    L = []

    #TODO chequear que estoy construyendo G2 y no MG2
    for origin, target, edge_data in messi_network.G2.edges(data=True):
        h = edge_data["reaction_constant"]
        #print("h")
        #print(h)
        if origin == target:
            continue


        #print("origin: %s, target: %s" % (origin, target))
        #print("G2 circle edges:")
        #print(messi_network.G2_circle.edges())

        for unique_simple_path in nx.all_simple_paths(messi_network.G2_circle, target, origin):
            first_edge_from_simple_path = [unique_simple_path[0], unique_simple_path[1]]
            i = origin
            m = get_label_from_edge(messi_network.G2_circle, first_edge_from_simple_path, "reaction_constant")
            j = target   

            item = [[h, i], [m, j]]
            #print("Item2")
            #print(item)
            
            L.append(item)


            #Deberia haber uno solo
            break;

    #print("L2")
    #print(L)
    return L


def get_pairs_of_binomial_exponents(messi_network):
    #We assume system is s-toric

    #Phi
    L1 = get_pairs_of_binomial_exponents_of_type_1(messi_network)
   
    #Simple path
    L2 = get_pairs_of_binomial_exponents_of_type_2(messi_network)

    #print("L2")
    #print(L2)
    return L1 + L2

def vec(messi_network, binomial):
    v = [0] * len(messi_network.species)
    #print(binomial)
    #print(messi_network.G.nodes())
    #binomial puede ser, por ejemplos, [2,3]

    #binomial deberia ser un vector que indica
    #que variables estan prendidas
    #deberia ser de longitud 2?
    for var in binomial:
        if var != None:
            v[var] = 1

    #v deberia ser un vector de tamaño #species
    return v


def get_binomial_basis(messi_network):
    #TODO separar en 2, para testear mas facil.
    L = []
    pairs_of_binomial_exponents = get_pairs_of_binomial_exponents(messi_network)

    #En el test esto se va a romper.
            
    #
    for pair in pairs_of_binomial_exponents:
        #print("pair:")
        #print(pair)
        vec1 = vec(messi_network, pair[0])
        vec2 = vec(messi_network, pair[1])
        res = list(map(lambda x, y: x - y, vec1, vec2))

        L.append(res)

    return L

def build_binomial_matrix(messi_network):
    """
    builds the binomial matrix B, as described in the MESSI paper
    """
    binomial_basis = get_binomial_basis(messi_network)
    #print("Binomial basis: ")
    #print(binomial_basis)

    #Las filas de B tienen los binomios
    column_basis =  extract_column_basis(np.array(binomial_basis).transpose())

    #print("Column basis: ")
    #print(column_basis)
    return column_basis.transpose()
    #return get_positive_matrix(column_basis.transpose())

#Aca se esta armando el ortogonal de la matriz B
#esta funcion esta mal
def build_integer_basis_of_orthogonal_complement_of_binomial_matrix(messi_network):
    """
    builds orthogonal complement of the binomial matrix - needed for Sigma_perp
    """
    binomial_matrix = build_binomial_matrix(messi_network)
    #print("binomial matrix")
    #print(binomial_matrix)
    ret =  build_integer_basis_matrix_of_orthogonal_complement_of_matrix(binomial_matrix.transpose(), positive=False).transpose()
    #M = binomial_matrix.transpose()
    #M_ort = ret.transpose()
    #print(M @ M_ort)
    return ret