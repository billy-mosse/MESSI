#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
#from numpy.linalg import det

from scipy import linalg, matrix
import scipy


debug=False

#This function is not mine! I stole it.
#Billy.
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
  
#This function is also not mine!
#Billy.
def getRsubsets(arr, r):
    """The main function that 
    prints all combinations 
    of  size r in arr[] of 
    size n. This function  
    mainly uses combinationUtil() 
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
    combinationUtil(arr, n, r,  
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

#Asume que son del mismo tamanio
def isConformal(orthant, circuit):
    for i in range(0, len(orthant)):
        
        #No son conformes si tienen signos opuestos, distintos de cero, en alguna coordenada.
        if orthant[i] * circuit[i] < 0:
            return False
    return True

#def sign(x):
#      return 1-(x<=0)

def hasEqualSign(orthant, circuit):
    """
    returns true iff the circuit parameter has equal sign as the orthant
    """
    for i in range(0, len(orthant)):
        
        #They DON'T have equal sign if they have opposite nonzero signs, or if XOR is false
        #

        if (orthant[i] * circuit[i] < 0) or (bool(orthant[i]) != bool(citcuit[i])):
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

def orthantHasConformalCircuit(orthant, sign_information_M):
    """unused function that checks if the orthant parameter is conformal
    to any of M's circuits"""

    for circuit in sign_information_M.circuits:
        if isConformal(orthant, circuit):
            if debug:
                print("The orthant %s is conformal to the circuit %s" % (orthant, circuit))
            return True
    return False


def getKappa(x1, x2):
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
    

    print("The kappas which are solution of f(x1, k) = 0 are (in columns):")
    ARR = scipy.linalg.null_space(M)
    print(ARR)

    input()

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
    var = input("We just gotta check that f(x2,k) is also zero.")
    #var = input("Ok, let's do it!")
    var = input("The following matrix A was calculated by hand.")
    var = input("The condition f(x2,k)=0 is equivalent to A k^t=0")
    print(A)
    var = input("On the other hand, the first kappa is:")
    print(nk)
    print(input("And the product, ladies and gentelmen, is:"))
    z = A.dot(nk)
    print(z)
    if np.linalg.norm(A.dot(nk) < 0.0001):
        print("Mazel tov!")
    else:
        print("Keep trying...")

