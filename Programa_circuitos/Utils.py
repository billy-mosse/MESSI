#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
#from numpy.linalg import det

debug=False
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


def hasEqualSign(orthant, circuit):
    for i in range(0, len(orthant)):
        
        #No tienen el mismo signo si tienen signos opuestos distintos de cero,
        #O alguno es cero y el otro no 
        if (orthant[i] * circuit[i] < 0) or (orthant[i] == 0 and circuit[i] != 0) or (orthant[i] != 0 and circuit[i] == 0):
            return False
    return True

#Esta funcion hace la union de unos conformal circuits.
#La union esta definida coordenada a coordenada:
#si algun circuito tiene la coordenada i distinta de cero,
#la union tiene esa coordenada igual a la del circuito.

#La union podria ser la suma, tranquilamente, pero prefiero hacer esto.
def union(conformal_circuits):
    if len(conformal_circuits) == 0:
        return None
    
    #No se si se deberia llamar s o distinto...
    s = len(conformal_circuits[0])
    U = [0 for i in range(s)]
    for j in range(0,s):
        for circuit in conformal_circuits:
            if circuit[j] != 0:
                U[j]=circuit[j]
                break
    return U

#Esta funcion chequea si el ortante orthant
#tiene algun circuito conforme a los circuitos de M
#
#Al final no la necesito.
def orthantHasConformalCircuit(orthant, sign_information_M):
    for circuit in sign_information_M.circuits:
        if isConformal(orthant, circuit):
            if debug:
                print "El ortante %s es conforme al circuito %s" % (orthant, circuit)
            return True
    return False