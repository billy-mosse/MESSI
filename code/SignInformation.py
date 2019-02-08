import Utils

import numpy as np
import copy

class SignInformation:
    def getCircuits(self, M):

        d = np.shape(self.matrix)[0]
        s= np.shape(self.matrix)[1]

        #Para que vaya de 1 a s.
        SList = range(1,s+1)
        circuits = []

        print(M)
        #Obtengo todos los subconjuntos de tamanio d-1, donde d es la cantidad de filas de la matriz
        L = Utils.getRsubsets(SList,d-1)
        for J in L:
            circuit = []
            for l in SList:
                if l in J:
                    circuit.append(0)
                else:
                    Jl = copy.copy(J)
                    Jl.append(l)

                    Jl.sort()

                    fl =(-1)**mu(l,J) * Utils.minor(self.matrix, Jl) 

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


    def __init__(self, M):
        self.matrix = M
        self.circuits = self.getCircuits(M)


    def get_conformal_circuits(self, orthant):
        conformal_circuits = []
        for circuit in self.circuits:
            if Utils.isConformal(orthant, circuit):
                conformal_circuits.append(circuit)

        return conformal_circuits


def mu(l, J):
    rta = sum(1 for j in J if j < l) % 2
    return rta
