from utils import Utils
import numpy as np
import copy
import itertools


class CircuitsInformation:

    def __init__(self, M):
        #TODO assert rank M is max

        self.matrix = M
        self.circuits = self.get_circuits(M)

    def get_circuits(self, M):
        """
        Gets all the circuits of a matrix, as explained in Lemma 50 of the paper

        Parameters
        ----------
        M: numpy.array
            A matrix

        returns: a list of circuits of the matrix

        """

        d = np.shape(self.matrix)[0]
        s= np.shape(self.matrix)[1]

        #range(1,s+1) returns [1,..., s]
        SList = range(1,s+1)
        circuits = []

        #Obtengo todos los subconjuntos de tamanio d-1, donde d es la cantidad de filas de la matriz
        L = Utils.get_r_subsets(SList,d-1)
        for J in L:
            circuit = []
            for l in SList:
                if l in J:
                    circuit.append(0)
                else:
                    Jl = copy.copy(J)
                    Jl.append(l)

                    Jl.sort()

                    fl =(-1)**self.mu(l,J) * Utils.minor(self.matrix, Jl)

                    '''if 4 in J and l==2:
                        print Jl
                        print mu(l,J)
                        print Utils.minor(self.matrix, sort(Jl))
                    '''

                    circuit.append(int(round(fl,0)))
            #print circuit

            circuits.append(circuit)

            #var = input("Press ENTER to continue")

        #Agrego a los opuestos
        all_circuits = []
        for circuit in circuits:
            if circuit not in all_circuits:
                all_circuits.append(circuit)
            opposite_circuit = [(-1)*x for x in circuit]
            if opposite_circuit not in all_circuits:
                all_circuits.append(opposite_circuit)

        #Elimino a repetidos de una manera probablemente muy costosa
        #print all_circuits
        #var = input("Press ENTER to continue")
        return all_circuits


    def get_conformal_circuits(self, orthant):
        """
        Gets all circuits that are conformal
        to the orthant received as a parameter
        """

        conformal_circuits = []
        for circuit in self.circuits:
            if Utils.is_conformal(orthant, circuit):
                conformal_circuits.append(circuit)

        return conformal_circuits


    def mu(self, l, J):
        """
        The mu function, as described in the MESSI paper
        """

        ret = sum(1 for j in J if j < l) % 2
        return ret


def get_only_one_equal_sign_vector(s, circuits_information_M1, circuits_information_M2):
    """
    Gets 2 vectors with equal signs from the circuits information of M1 and M2

    Brute-force searches through all orthants for 2 vectors we know exist with equal signs and returns the first one

    TODO: optimize! This is an exponetial algorithm.
    """

    equal_sign_vectors = []

    #TODO: create them dynamically
    orthants = list(itertools.product([-1,0,1],repeat=s))

    #Steps 2,3,4
    for orthant in orthants:
        #print(orthant)
        conformal_circuits_M1 = circuits_information_M1.get_conformal_circuits(orthant)

        is_zero_or_invalid = True
        #union of the circuits conformal to the orthant
        U_M1 = Utils.union(conformal_circuits_M1)


        is_zero_or_invalid = True
        if U_M1 != None:
            for elt in U_M1:
                if elt != 0:
                    is_zero_or_invalid = False

        #print(U_M1)


        if U_M1 != None and Utils.has_equal_sign(orthant, U_M1) and not is_zero_or_invalid:

            #Si el ortante tiene soporte igual a la union de los circuitos de M1, sigo con el paso 3 del algoritmo.
            conformal_circuits_M2 = circuits_information_M2.get_conformal_circuits(orthant)


            U_M2 = Utils.union(conformal_circuits_M2)


            if U_M2 != None and Utils.has_equal_sign(orthant, U_M2):
                equal_sign_vectors.append([U_M1, U_M2])

                print("Two vectors with the same sign, corresponding to the orthant %s, are %s, from M1, and %s, from M2." % (orthant, U_M1, U_M2))
                print("As requested, now stopping.")

                return equal_sign_vectors

        else:
            continue
            #Isn't useful
    return equal_sign_vectors

#Esta función parece estar bien, ver tests.ipynb
def get_all_equal_sign_orthants(s, circuits_information_M1, circuits_information_M2):
    """
    Gets 2 vectors with equal signs from the circuits information of M1 and M2

    Brute-force searches through all orthants for 2 vectors we know exist with equal signs

    TODO: optimize! This is an exponetial algorithm.
    """

    equal_sign_vectors = []

    #TODO: create them dynamically
    orthants = list(itertools.product([-1,0,1],repeat=s))

    solution_orthants = []
    #Steps 2,3,4
    for orthant in orthants:
        conformal_circuits_M1 = circuits_information_M1.get_conformal_circuits(orthant)

        #union of the circuits conformal to the orthant
        U_M1 = Utils.union(conformal_circuits_M1)

        is_zero_or_invalid = True
        if U_M1 != None:
            for elt in U_M1:
                if elt != 0:
                    is_zero_or_invalid = False


        #print(U_M1)

        if U_M1 != None and not is_zero_or_invalid and Utils.has_equal_sign(orthant, U_M1):
            #Si el ortante tiene soporte igual a la union de los circuitos de M1, sigo con el paso 3 del algoritmo.
            conformal_circuits_M2 = circuits_information_M2.get_conformal_circuits(orthant)
            U_M2 = Utils.union(conformal_circuits_M2)
            if U_M2 != None and Utils.has_equal_sign(orthant, U_M2):
                #equal_sign_vectors.append([U_M1, U_M2])
                solution_orthants.append(orthant)
                equal_sign_vectors.append([U_M1, U_M2])

                #Stop searching in this orthant
                continue

                print("Two vectors with the same sign, corresponding to the orthant %s, are %s, from M1, and %s, from M2." % (orthant, U_M1, U_M2))
        else:
            continue
            #Isn't useful
    return equal_sign_vectors

def get_count_equal_sign_orthants(s, circuits_information_M1, circuits_information_M2):
    """
    Gets 2 vectors with equal signs from the circuits information of M1 and M2

    Brute-force searches through all orthants for 2 vectors we know exist with equal signs

    TODO: optimize! This is an exponetial algorithm.
    """

    equal_sign_vectors = []

    #TODO: create them dynamically
    orthants = list(itertools.product([-1,0,1],repeat=s))

    solution_orthants = 0
    #Steps 2,3,4
    for orthant in orthants:
        conformal_circuits_M1 = circuits_information_M1.get_conformal_circuits(orthant)

        #union of the circuits conformal to the orthant
        U_M1 = Utils.union(conformal_circuits_M1)

        is_zero_or_invalid = True
        if U_M1 != None:
            for elt in U_M1:
                if elt != 0:
                    is_zero_or_invalid = False


        #print(U_M1)

        if U_M1 != None and not is_zero_or_invalid and Utils.has_equal_sign(orthant, U_M1):
            #Si el ortante tiene soporte igual a la union de los circuitos de M1, sigo con el paso 3 del algoritmo.
            conformal_circuits_M2 = circuits_information_M2.get_conformal_circuits(orthant)
            U_M2 = Utils.union(conformal_circuits_M2)
            if U_M2 != None and Utils.has_equal_sign(orthant, U_M2):
                #equal_sign_vectors.append([U_M1, U_M2])
                solution_orthants += 1
                #equal_sign_vectors.append([U_M1, U_M2])

                #Stop searching in this orthant
                continue

                print("Two vectors with the same sign, corresponding to the orthant %s, are %s, from M1, and %s, from M2." % (orthant, U_M1, U_M2))
        else:
            continue
            #Isn't useful
    return solution_orthants



def get_equal_sign_orthants(s, circuits_information_M1, circuits_information_M2):
    """
    Gets 2 vectors with equal signs from the circuits information of M1 and M2

    Brute-force searches through all orthants for 2 vectors we know exist with equal signs

    TODO: optimize! This is an exponetial algorithm.
    """

    equal_sign_vectors = []

    #TODO: create them dynamically
    orthants = list(itertools.product([-1,0,1],repeat=s))

    solution_orthants = []
    #Steps 2,3,4
    for orthant in orthants:
        print(orthant)
        conformal_circuits_M1 = circuits_information_M1.get_conformal_circuits(orthant)

        #union of the circuits conformal to the orthant
        U_M1 = Utils.union(conformal_circuits_M1)

        is_zero_or_invalid = True
        if U_M1 != None:
            for elt in U_M1:
                if elt != 0:
                    is_zero_or_invalid = False


        #print(U_M1)

        if U_M1 != None and not is_zero_or_invalid and Utils.has_equal_sign(orthant, U_M1):
            #Si el ortante tiene soporte igual a la union de los circuitos de M1, sigo con el paso 3 del algoritmo.
            conformal_circuits_M2 = circuits_information_M2.get_conformal_circuits(orthant)
            U_M2 = Utils.union(conformal_circuits_M2)
            if U_M2 != None and Utils.has_equal_sign(orthant, U_M2):
                #equal_sign_vectors.append([U_M1, U_M2])
                solution_orthants.append(orthant)
                #equal_sign_vectors.append([U_M1, U_M2])

                #Stop searching in this orthant
                continue

                print("Two vectors with the same sign, corresponding to the orthant %s, are %s, from M1, and %s, from M2." % (orthant, U_M1, U_M2))
        else:
            continue
            #Isn't useful
    return solution_orthants


#Esta función parece estar bien, ver tests.ipynb
def get_unique_equal_sign_vector(s, circuits_information_M1, circuits_information_M2):
    """
    Gets 2 vectors with equal signs from the circuits information of M1 and M2

    Brute-force searches through all orthants for 2 vectors we know exist with equal signs

    TODO: optimize! This is an exponetial algorithm.
    """

    equal_sign_vectors = []

    #TODO: create them dynamically
    orthants = list(itertools.product([-1,0,1],repeat=s))

    solution_orthants = []
    #Steps 2,3,4
    for orthant in orthants:
        conformal_circuits_M1 = circuits_information_M1.get_conformal_circuits(orthant)

        #union of the circuits conformal to the orthant
        U_M1 = Utils.union(conformal_circuits_M1)

        is_zero_or_invalid = True
        if U_M1 != None:
            for elt in U_M1:
                if elt != 0:
                    is_zero_or_invalid = False


        #print(U_M1)

        if U_M1 != None and not is_zero_or_invalid and Utils.has_equal_sign(orthant, U_M1):
            #Si el ortante tiene soporte igual a la union de los circuitos de M1, sigo con el paso 3 del algoritmo.
            conformal_circuits_M2 = circuits_information_M2.get_conformal_circuits(orthant)
            U_M2 = Utils.union(conformal_circuits_M2)
            if U_M2 != None and Utils.has_equal_sign(orthant, U_M2):
                #equal_sign_vectors.append([U_M1, U_M2])
                solution_orthants.append(orthant)
                if len(solution_orthants) > 1:
                    return None

                #Stop searching in this orthant
                continue

                print("Two vectors with the same sign, corresponding to the orthant %s, are %s, from M1, and %s, from M2." % (orthant, U_M1, U_M2))
        else:
            continue
            #Isn't useful
    if len(solution_orthants) == 1:
        return solution_orthants[0]


#Esta función parece estar bien, ver tests.ipynb
def get_equal_sign_vectors(s, circuits_information_M1, circuits_information_M2):
    """
    Gets 2 vectors with equal signs from the circuits information of M1 and M2

    Brute-force searches through all orthants for 2 vectors we know exist with equal signs

    TODO: optimize! This is an exponetial algorithm.
    """

    equal_sign_vectors = []

    #TODO: create them dynamically
    orthants = list(itertools.product([-1,0,1],repeat=s))

    #Steps 2,3,4
    for orthant in orthants:
        conformal_circuits_M1 = circuits_information_M1.get_conformal_circuits(orthant)

        #union of the circuits conformal to the orthant
        U_M1 = Utils.union(conformal_circuits_M1)

        is_zero_or_invalid = True
        if U_M1 != None:
            for elt in U_M1:
                if elt != 0:
                    is_zero_or_invalid = False

        #print(U_M1)

        if U_M1 != None and not is_zero_or_invalid and Utils.has_equal_sign(orthant, U_M1):
            #Si el ortante tiene soporte igual a la union de los circuitos de M1, sigo con el paso 3 del algoritmo.
            conformal_circuits_M2 = circuits_information_M2.get_conformal_circuits(orthant)
            U_M2 = Utils.union(conformal_circuits_M2)
            if U_M2 != None and Utils.has_equal_sign(orthant, U_M2):
                equal_sign_vectors.append([U_M1, U_M2])

                print("Two vectors with the same sign, corresponding to the orthant %s, are %s, from M1, and %s, from M2." % (orthant, U_M1, U_M2))
        else:
            continue
            #Isn't useful
    return equal_sign_vectors


def get_positive_vector_in_matrix_rows(matrix):
    #print("supuetamente las filas con binomios")
    #print(matrix)
    circuit_information = CircuitsInformation(matrix)
    nonnegative_circuit = [0] * np.shape(matrix)[1]
    for circuit in circuit_information.circuits:
        #print(circuit)
        if Utils.is_nonnegative(circuit):
            nonnegative_circuit = np.add(nonnegative_circuit, circuit)
    return nonnegative_circuit

"""
def get_positive_vector_in_kernel(matrix):
    #Los tiene en las columnas
    kernel = MatrixUtils.build_integer_basis_matrix_of_orthogonal_complement_of_matrix(matrix)
    d = shape(kernel)[0]
    n = shape(kernel)[1]
    J = Utils.nchoosek(range(1, n+1), d-1)
    W = []
    #We look for circuits of M:
    for i in range(1, len(J)):
        if numpy.matrix_rank(kernel[:][J[i]])==d-1:
            w = []
            for l in range(1, n+1):
                w=[w,(-1)^(sum(j for j in J if j < l))*numpy.det(M[:][sort([J(i,:),l])])];
            W = [W, w] #?



for i=1:size(J,1)
  if rank(M(:,J(i,:)))==d-1
    w=[];
    for l=1:n
      w=[w,(-1)^(sum(J(i,:)<l))*det(M(:,sort([J(i,:),l])))];
    endfor %l
    W=[W;w];;%we keep the circiuts of M in each row of W
  endif
endfor
%we make the first nonzero coordinate of each circuit positive:
for i=1:size(W,1)
  ww=find(W(i,:));
  if length(ww)!=0
    if W(i,ww(1))<0
      W(i,:)=-W(i,:);
    endif
  endif
endfor
%we keep in W1 different circuits of M with first nonzero coordinate positive:
W1=W(1,:);
for i=2:size(W,1)
  l=size(W1,1);
  m1=[];
  for k=1:l
    m1(k)=prod(abs(W(i,:)-W1(k,:))<1e-14);
  endfor
  if sum(m1)==0
    W1=[W1;W(i,:)];
  endif
endfor
%we keep in W2 different circuits of M with all of its coordinates positive:
W2=[];
for i=1:size(W1,1)
  if sum(W1(i,:)<-1e-14)==0
    W2 = [W2;W1(i,:)];
  endif
endfor
v=sum(W2);
Nv=N*v';"""
