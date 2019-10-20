#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
#from numpy.linalg import det

from scipy import linalg, matrix
import scipy


debug=False

#This function is not mine! I stole it.
#Billy.
def combination_util(arr, n, r,  
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
    combination_util(arr, n, r,  
                    index + 1, data, i + 1, L) 
      
    # current is excluded,  
    # replace it with 
    # next (Note that i+1  
    # is passed, but index  
    # is not changed) 
    combination_util(arr, n, r, index,  
                    data, i + 1, L) 
  
#This function is also not mine!
#Billy.
def get_r_subsets(arr, r):
    """The main function that 
    prints all combinations 
    of  size r in arr[] of 
    size n. This function  
    mainly uses combination_util() 
    """
    n = len(arr) 
    L = []
    # A temporary array to 
    # store all combination 
    # one by one 
    data = list(range(r)) 
      
    # Print all combination  
    # using temporary  
    # array 'data[]' 
    combination_util(arr, n, r,  
                    0, data, 0, L) 

    #print L
    return L


#returns indexes for matrix
def i(cols):
    return [col-1 for col in cols]


#La cantidad de columnas es tal que cuando sliceamos, obtenemos una matriz cuadrada 
#OJO que los indices empiezan en cero
def minor(A, cols):

    A_sliced = A[:, i(cols)]
    #Deberia ser cuadrada
    assert np.shape(A_sliced)[0] == np.shape(A_sliced)[1]
    return np.linalg.det(A_sliced)

#Asume que son del mismo tamaño
def is_conformal(orthant, circuit):
    for i in range(0, len(orthant)):
        
        #No son conformes si tienen signos opuestos, distintos de cero, en alguna coordenada.
        if orthant[i] * circuit[i] < 0:
            return False
    return True

#def sign(x):
#      return 1-(x<=0)

def has_equal_sign(orthant, circuit):
    """
    returns true iff the circuit parameter has equal sign as the orthant
    """
    for i in range(0, len(orthant)):
        
        #They DON'T have equal sign if they have opposite nonzero signs, or if XOR is false
        #

        if (orthant[i] * circuit[i] < 0) or (bool(orthant[i]) != bool(circuit[i])):
            return False

        #if (orthant[i] * circuit[i] < 0) or (orthant[i] == 0 and circuit[i] != 0) or (orthant[i] != 0 and circuit[i] == 0):
        #    return False
    return True

#Esta funcion hace la union de unos conformal circuits.
#La union esta definida coordenada a coordenada:
#si algun circuito tiene la coordenada i distinta de cero,
#la union tiene esa coordenada igual a la del circuito.

def union(conformal_circuits):
    """returns the sum of the conformal circuits list.
    """
    if len(conformal_circuits) == 0:
        return None
    
    #No se si se deberia llamar s o distinto...
    s = len(conformal_circuits[0])
    U = [0 for i in range(s)]
    for j in range(0,s):
        for circuit in conformal_circuits:
            U[j]=U[j]+circuit[j]
    return U

'''
useless function
def union(conformalCircuit1,conformalCircuit2):
    """Returns the union of 2 cicuits.
    The circuits are assumed to be conformal.
    """
    ret = []
    for i in range(0, len(circ1)):
        c1i = conformalCircuit1[i]
        c2i = conformalCircuit2[i]
        if c1i==1 or c2i==1:
            ret.append(1)
        elif c1i==-1 or c2i==-1:
            ret.append(-1)
        else:
            ret.append(0)
    return ret'''

def orthant_has_conformal_circuit(orthant, sign_information_M):
    """unused function that checks if the orthant parameter is conformal
    to any of M's circuits"""

    for circuit in sign_information_M.circuits:
        if is_conformal(orthant, circuit):
            if debug:
                print("The orthant %s is conformal to the circuit %s" % (orthant, circuit))
            return True
    return False


#La Phi parece estar bien
def Phi(x, educt_complexes_matrix):

    #print("educt complexes matrix", educt_complexes_matrix)
    #print(educt_complexes_matrix)
    #Educt complexes matrix itene a y_k en cada columna.
    #phi(x1) deberia tener en la primera coordenada al producto de los momomios correspondientes
    ret = []
    for column in educt_complexes_matrix.transpose():
        #print("column")
        #print(column)
        #print("x")
        #print(x)
        #Es 1 o 0?
        num = 1
        for index, x_i in enumerate(x):
            if column[index]:
                num = num*pow(x_i, column[index])
        #print("num")
        #print(num)
        ret.append(num)
    return ret


#La educt complexes matrix tambien parece estar bien
#Phi tambien parece estar bien
def get_kappa2(x1, x2, positive_Mperp, educt_complexes_matrix, messi_network, toric_N):
    #k = diag(Phi(x))\m M \lambda
    #print("x1, x2:")
    #print(x1)
    #print(x2)

    #positive M perp tiene un error...
    #positive_Mperp = positive_Mperp.transpose()

    p1 = Phi(x1, educt_complexes_matrix)

    #No hace falta calcular los dos
    p2 = Phi(x2, educt_complexes_matrix)

    if False:
        print("p1: ")
        print(p1)

        print("p2: ")
        print(p2)


    temp_matrix = np.linalg.inv(np.diag(p1))
    size_l = np.shape(positive_Mperp)[1]
    l = np.array([[1] * size_l]).transpose()


    """print("Experimento nuevo!")
    print(positive_Mperp @ l)

    print(temp_matrix)
    print(positive_Mperp)
    print(l)

    print("temp matrix")
    print(np.shape(temp_matrix))

    print("p1")
    print(np.shape(p1))

    print("positive m perp")
    print(np.shape(positive_Mperp))

    print("l")
    print(np.shape(l))"""

    if True:
        print("temp matrix", temp_matrix)
        print("positive M perp", positive_Mperp)
        print("l", l)

    k = temp_matrix @ positive_Mperp @ l

    list_k = list(k)
    #for k_i in list_k:
    #    if k_i < 0:
    #        #Alguno da cero
    #        return False

    print("k")
    print(k)
    #input("")

    n = np.shape(p2)[0]
    k_vec = []
    for vec in list_k:
        for j in vec:
            k_vec.append(j)

    matrix_k = np.diag(k_vec)

    #print(matrix_k)

    #La cuenta está bien, porque esto da cero
    res1 = toric_N @ matrix_k @ p1

    #Pero esto no da cero...
    res2 = toric_N @ matrix_k @ p2


    if np.linalg.norm(res1) < 0.0001 and\
    np.linalg.norm(res2) < 0.0001:
        print("These are indeed multistationarity witnesses.")
        return True
    else:
        print("Error")
        return False

#deprecated
def get_kappa(x1, x2):
    '''M = matrix([
    [-x1[0] * x1[8], x1[4], 0, 0, 0, 0, 0, 0, 0, -x1[3] * x1[0], -x1[7], 0],
    [0, 0, x1[4], -x1[1] * x1[0],  x1[5], 0, -x1[2]*x1[1], x1[6], x1[6], 0, 0, 0],
    [0, 0, 0, 0, 0, 0, -x1[2] * x1[1], x1[6], 0, 0, 0, -x1[7]],
    [0, 0, 0, 0, 0, 0, 0, 0, x1[6], x1[3]*x1[0], x1[7], 0],
    [x1[0]*x1[8], -x1[4], -x1[4], 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, x1[1]*x1[0], -x1[5], -x1[5], 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, x1[2]*x1[1], -x1[6], -x1[6], 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, x1[3]*x1[0], -x1[7], x1[7]],
    [-x1[0]*x1[8], x1[4], x1[4], 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -x1[1]*x1[0], x1[5], x1[5], 0, 0, 0, 0, 0, x1[7]]
    ])'''


    #TODO: this shouldn't be hardcoded.
    M = matrix([
    [-x1[0]*x1[8], x1[4], 0, 0, 0, x1[5], 0, 0, 0, 0, 0, 0], #1
    [0, 0, x1[4], -x1[1]*x1[9], x1[5], 0, -x1[2]*x1[1], x1[6], x1[6], 0, 0, 0], #2
    [0, 0, 0, 0, 0, 0, -x1[2]*x1[1], x1[6], 0, 0, 0, x1[7]], #3
    [0, 0, 0, 0, 0, 0, 0, 0, x1[6], -x1[3]*x1[9], x1[7], 0], #4
    [x1[0]*x1[8], -x1[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #5
    [0, 0, 0, x1[1]*x1[9], -x1[5], -x1[5], 0, 0, 0, 0, 0, 0],#6
    [0, 0, 0, 0, 0, 0, x1[2]*x1[1], -x1[6], -x1[6], 0, 0, 0],#7
    [0, 0, 0, 0, 0, 0, 0, 0, 0, x1[3]*x1[9], -x1[7], -x1[7]], #8
    [-x1[0]*x1[8], x1[4], x1[4], 0, 0, 0, 0, 0, 0, 0, 0, 0], #9
    [0, 0, 0, -x1[1]*x1[9], x1[5], x1[5], 0, 0, 0, -x1[3]*x1[9], x1[7], x1[7]]])#10
    

    #HAY QUE ASEGURAR QUE SEAN POSITIVOS.
    #https://arxiv.org/abs/1102.1590
    print("The kappas which are solution of f(x1, k) = 0 are (in columns):")
    ARR = scipy.linalg.null_space(M)
    #print(ARR)

    #print("Why should be the vectors in the nullspace positive?")

    print("Press ENTER to continue")

    k = []
    for i in range(0, np.shape(ARR)[0]):
        k.append(ARR[i][0])

    #TODO: this shouldn't be hardcoded
    A = matrix([
    [-x2[0]*x2[8], x2[4], 0, 0, 0, x2[5], 0, 0, 0, 0, 0, 0], #1
    [0, 0, x2[4], -x2[1]*x2[9], x2[5], 0, -x2[2]*x2[1], x2[6], x2[6], 0, 0, 0], #2
    [0, 0, 0, 0, 0, 0, -x2[2]*x2[1], x2[6], 0, 0, 0, x2[7]], #3
    [0, 0, 0, 0, 0, 0, 0, 0, x2[6], -x2[3]*x2[9], x2[7], 0], #4
    [x2[0]*x2[8], -x2[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #5
    [0, 0, 0, x2[1]*x2[9], -x2[5], -x2[5], 0, 0, 0, 0, 0, 0],#6
    [0, 0, 0, 0, 0, 0, x2[2]*x2[1], -x2[6], -x2[6], 0, 0, 0],#7
    [0, 0, 0, 0, 0, 0, 0, 0, 0, x2[3]*x2[9], -x2[7], -x2[7]], #8
    [-x2[0]*x2[8], x2[4], x2[4], 0, 0, 0, 0, 0, 0, 0, 0, 0], #9
    [0, 0, 0, -x2[1]*x2[9], x2[5], x2[5], 0, 0, 0, -x2[3]*x2[9], x2[7], x2[7]]])#10

    nk = np.array(k)
    print("")
    #input("The moment of truth has arrived.")
    input("We must check that both x1 and x2 are (positive) steady states.")
    input("Let's do it with the first choice of kappa, i.e., the first column.")
    #var = input("We just gotta check that f(x2,k) is also zero.")
    #var = input("Ok, let's do it!")
    #var = input("The following matrix A was calculated by hand.")
    #var = input("The condition f(x2,k)=0 is equivalent to A k^t=0")
    #print(A)
    
    #var = input("kappa is:")
    #print(nk)

    #print(input("Product:"))
    z = A.dot(nk)
    #print(z)

    if np.linalg.norm(A.dot(nk) < 0.0001):
        print("f(x2, k) = 0. It's a steady state!")
        return nk
    else:
        print("f(x2, k) != 0. It's not a steady state.")
        return None




#Theorem 5.8
#Esto deberia andar porque ya andaba para el final de Alicia
def get_multistationarity_witnesses(w, v, s, d):
    print("w", w)
    print("v", v)
    x1 = []
    x2 = []
    for i in range (0, s):
        if v[i] != 0:
            x1.append(
                (w[i] * 1.0) / (np.exp(v[i])-1.0))
        else:
            x1.append(1)

        x2.append(np.exp(v[i]) * x1[i])
    
    """print("Hola")
    print("v", v)
    print("Hola2")
    print("w", w)

    print("x1",x1)
    
    print("x2", x2)
    exit(0)"""

    return x1,x2



def nchoosek(n, k):
    if k == 0:
        r = 1
    else:
        r = n/k * nchoosek(n-1, k-1)
    return round(r)



def is_nonnegative(vec):
    for v_i in vec:
        if v_i < 0:
            return False
    return True