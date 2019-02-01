#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Ver https://wiki.sagemath.org/Python3-compatible%20code para sage

from __future__ import division, absolute_import, print_function

import sys

import numpy as np

from numpy.linalg import det
from numpy import shape
import copy

import Utils
from SignInformation import SignInformation
import itertools
import MESSIUtils

#NICETOHAVE: GUI para ingresar un grafo y que se fije si admite estructura MESSI y haga todas las cuentitas.

#Objetivo para hoy: hacer que calcule los circuitos.
#Objetivo secundario: hacer la otra parte.

#Pregunta: tengo que guardar la info de los circuitos, no? Para calcular los witnesses despues

#Deberia hacer la union "circuital" de vectores

#./sage -python circuits.py
def union(circ1,circ2):
    d = []
    for i in range(0, len(circ1)):
        c1 = circ1[i]
        c2 = circ2[i]
        if c1==1 or c2==1:
            d.append(1)
        elif c1==-1 or c2==-1:
            d.append(-1)
        else:
            d.append(0)
    return d

#TODO: podrian ser SETS mejor.
def getFullSigmaSubperp(Bperp, Mt, s, d):

    signs = []
    witnesses = []
    
    #Bperp tiene d filas y s columnas.
    #Mt tiene s-d filas y s columnas.
    SList = range(1,s+1)

    #Obtengo todos los subconjuntos de tamanio d
    L = Utils.getRsubsets(SList,d)

    S = set(SList)
    Sigma_subperp = [["J^c", "menor M^t", "J", "menor B^\\perp", "Result"]]
    for JList in L:
        J = set(JList)
        #Jc = set(S)[x for x in S if x not in J]
        Jc = S - J
        multiplier = (-1)**sum(J)

        #No hay nada mas horrible que transformar una lista a un set y luego de nuevo a una lista :-)
        
        JL = list(J)
        JL.sort()
        JcL = list(Jc)
        JcL.sort()


        mB = Utils.minor(Bperp,JL)
        mM = Utils.minor(Mt, JcL)
        result = int(round(mB * mM * multiplier,0)) #son todos enteros, ver paper.
        if(result > 0): #and '+' not in signs:
            signs.append('+')
            witnesses.append(['+', J, Jc])
        elif(result < 0):# and '-' not in signs:
            signs.append('-')
            witnesses.append(['-', J, Jc])
        Sigma_subperp.append([Jc, mM, J, mB, result])

    
    if len(signs) > 1:
        print "Sigma_perp is mixed! Great, then we can find v in T^perp and w in S with same sign."
        #print witnesses
    else:
        print "Sigma_perp is NOT mixed"
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

def main(debug):

    relevant_matrices = MESSIUtils.getRelevantMatrices(debug)


    Bperp = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, -1, 0],
    [0, 1, 0, 0, 2, 2, 1, 1, 2, 1],
    [-0, 0, 1, 0, 1, 1, 1, 1, 1, 1],
    [-0, 0, 0, 1, -1, -1, 0, -0, -1, -1]])

    Mt = np.array([
    [1, 0, 0, 0, 0, 0, -1, 1, 0, -1],
    [0, 1, 0, 0, 0, 0, -1, 1, 0, -1],
    [0, 0, 1, 0, 0, 0, 0, -1, 0, 1],
    [0, 0, 0, 1, 0, 0, 0, -1, 0, 1],
    [0, 0, 0, 0, 1, 0, -1, 1, -1, -1],
    [0, 0, 0, 0, 0, 1, -1, 1, -0, -2]])


    s = np.shape(Bperp)[1]
    assert s == np.shape(Mt)[1]
    d = np.shape(Bperp)[0]
    assert np.shape(Mt)[0] == s-d

    #Bperp = np.array([[1, 0, -3, 2],
                     #[0, 1, -2, 1]])
    #Mt = np.array([[1, 1, 1, 1], [5, 2, 0, 5]])


    print ""
    print ""
    print "For the second part of the program, we will be working with toy matrices"
    print "B^perp is:"
    print Bperp
    print "M^t is:"
    print Mt

    assert d<=s

    print "First we will compute Sigma_perp to check if it's mixed."
    FULL_Sigma_subperp = getFullSigmaSubperp(Bperp, Mt, s, d)

    sign_information_Bperp = SignInformation(Bperp)
    sign_information_Mt = SignInformation(Mt)
    #print(sign_information_Bperp.circuits)

    if not debug or True:
        var = raw_input("Press ENTER to continue.")

    orthants = list(itertools.product([-1,0,1],repeat=s))



    for orthant in orthants:
        conformal_circuits_Bperp = sign_information_Bperp.get_conformal_circuits(orthant)
        #print "Los siguientes circuitos son conformes a..."
        #print orthant
        #print ".."
        #print conformal_circuits

        U_Bperp = Utils.union(conformal_circuits_Bperp)

        '''print "orthant:"
        print orthant
        print "U:"
        print U
        print "______"'''

        if U_Bperp != None and Utils.hasEqualSign(orthant, U_Bperp):
            #Si el ortante tiene soporte igual a la union de los circuitos de Bperp, sigo con el paso 3 del algoritmo.
            conformal_circuits_Mt = sign_information_Mt.get_conformal_circuits(orthant)
            U_Mt = Utils.union(conformal_circuits_Mt)
            if U_Mt != None and Utils.hasEqualSign(orthant, U_Mt):
                print "Two vectors with the same sign, corresponding to the orthant %s, are %s, from T^perp, and %s, from S." % (orthant, U_Bperp, U_Mt)
        else:
            continue
            #No hago nada.



    sign_information_Mt = SignInformation(Mt)
    #print(sign_information_Mt.circuits)

    #print_table(FULL_Sigma_subperp)


if __name__=="__main__":
    debug  = False
    if len(sys.argv) == 2:
        if sys.argv[1] == "debug":
            debug = True
    main(debug)
