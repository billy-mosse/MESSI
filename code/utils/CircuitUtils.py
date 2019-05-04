from utils import Utils
import numpy as np
import copy


class CircuitsInformation:

    def __init__(self, M):
        #TODO assert rank M is max
        
        self.matrix = M
        self.circuits = self.get_circuits(M)

    def get_circuits(self, M):
        """Gets all the circuits of matrix M
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
        """Gets all circuits that are conformal
        to the orthant received as a parameter"""

        conformal_circuits = []
        for circuit in self.circuits:
            if Utils.is_conformal(orthant, circuit):
                conformal_circuits.append(circuit)
        
        return conformal_circuits


    def mu(self, l, J):
        ret = sum(1 for j in J if j < l) % 2
        return ret

def get_equal_sign_vectors(s, circuits_information_Bperp, circuits_information_Mt):
    equal_sign_vectors = []

    #TODO: create them dynamically
    orthants = list(itertools.product([-1,0,1],repeat=s))

    #Steps 2,3,4
    for orthant in orthants:
        conformal_circuits_Bperp = circuits_information_Bperp.get_conformal_circuits(orthant)

        #union of the circuits conformal to the orthant
        U_Bperp = Utils.union(conformal_circuits_Bperp)
        
        if U_Bperp != None and Utils.has_equal_sign(orthant, U_Bperp):
            #Si el ortante tiene soporte igual a la union de los circuitos de Bperp, sigo con el paso 3 del algoritmo.
            conformal_circuits_Mt = circuits_information_Mt.get_conformal_circuits(orthant)
            U_Mt = Utils.union(conformal_circuits_Mt)
            if U_Mt != None and Utils.has_equal_sign(orthant, U_Mt):
                equal_sign_vectors.append([U_Bperp, U_Mt])
                print("Two vectors with the same sign, corresponding to the orthant %s, are %s, from T^perp, and %s, from S." % (orthant, U_Bperp, U_Mt))
        else:
            continue
            #Isn't useful
    return equal_sign_vectors