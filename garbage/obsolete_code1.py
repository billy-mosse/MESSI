import sys

import numpy as np

from numpy.linalg import det
from numpy import shape
import copy

#Pregunta: tengo que guardar la info de los circuitos, no? Para calcular los witnesses despues

def combinationUtil(arr, n, r,  
                    index, data, i, L): 
    # Current combination is  
    # ready to be printed, 
    # print it 
    if(index == r):
        L.append(list(data))
        return
  
    # When no more elements  
    # are there to put in data[] 
    if(i >= n):
        return
  
    # current is included,  
    # put next at next 
    # location  
    data[index] = arr[i] 
    combinationUtil(arr, n, r,  
                    index + 1, data, i + 1, L) 
      
    # current is excluded,  
    # replace it with 
    # next (Note that i+1  
    # is passed, but index  
    # is not changed) 
    combinationUtil(arr, n, r, index,  
                    data, i + 1, L) 
  
  
# The main function that 
# prints all combinations 
# of size r in arr[] of 
# size n. This function  
# mainly uses combinationUtil() 
def getRsubsets(arr, r):
    n = len(arr) 
    L = []
    # A temporary array to 
    # store all combination 
    # one by one 
    data = list(range(r)) 
      
    # Print all combination  
    # using temporary  
    # array 'data[]' 
    combinationUtil(arr, n, r,  
                    0, data, 0, L) 

    #print L
    return L

#Deberia hacer la union "circuital" de vectores
def union(circ1,circ2):
    print "hola"

#Devuelve una tabla?
def Sigma_perp(Mt,Bperp):
    print "hola"

#Calcula a mano la clausura de circuitos
def clau_circ(circs):
    print "hola"

def chequear_interseccion(circs1, circs2):
    print "hola"

#returns indexes for matrix
def i(cols):
    return [col-1 for col in cols]

#La cantidad de columnas es tal que cuando sliceamos, obtenemos una matriz cuadrada 
#OJO que los indices empiezan en cero
def minor(A, cols):
    A_sliced = A[:, i(cols)]

    #Deberia ser cuadrada
    assert np.shape(A_sliced)[0] == np.shape(A_sliced)[1]
    return det(A_sliced)
  

s = 4
d = 2

#TODO: podrian ser SETS mejor.
def get_full_Sigma_subperp(Bperp, Mt):

    signs = []
    witnesses = []
    print "entre a la funcion"
    
    #Bperp tiene d filas y s columnas.
    #Mt tiene s-d filas y s columnas.
    SList = range(1,s+1)

    #Obtengo todos los subconjuntos de tamanio d
    L = getRsubsets(SList,d)

    S = set(SList)
    Sigma_subperp = [["J^c", "menor M^t", "J", "menor B^\\perp", "Result"]]
    for JList in L:
        J = set(JList)
        #Jc = set(S)[x for x in S if x not in J]
        Jc = S - J
        multiplier = (-1)**sum(J)
        mB = minor(Bperp,J)
        mM = minor(Mt, Jc)
        result = mB * mM * multiplier
        if(result > 0): #and '+' not in signs:
            signs.append('+')
            witnesses.append(['+', J, Jc])
        elif(result < 0):# and '-' not in signs:
            signs.append('-')
            witnesses.append(['-', J, Jc])
        Sigma_subperp.append([Jc, mM, J, mB, result])

    
    if len(signs) > 1:
        print "Esta mezclada"
        #print witnesses
    return Sigma_subperp    
    #Faltaria armar la tabla, que es bastante grande, no? No se si la voy a armar
    #Pero si voy  devolver

    #for i in range(1, 2**d):

def print_table(table):
    longest_cols = [
        (max([len(str(row[i])) for row in table]) + 3)
        for i in range(len(table[0]))
    ]
    row_format = "".join(["{:>" + str(longest_col) + "}" for longest_col in longest_cols])
    for row in table:
        print(row_format.format(*row))

def main(args):
    print args[1]
    Bperp = np.array([[1, 0, -3, 2],
                     [0, 1, -2, 1]])
    Mt = np.array([[1, 1, 1, 1], [5, 2, 0, 5]])
    assert d<=s
    FULL_Sigma_subperp = get_full_Sigma_subperp(Bperp, Mt)
    
    print_table(FULL_Sigma_subperp)


if __name__=="__main__":
    main(sys.argv)