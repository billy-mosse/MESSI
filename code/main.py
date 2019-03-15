#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Check https://wiki.sagemath.org/Python3-compatible%20code
#if someday we want to use sage again

#This shouldn't be necessary, as we are using python 3 now
from __future__ import division, absolute_import, print_function

import sys
import os
import numpy as np

from numpy.linalg import det
from numpy import shape
import copy

#our library
import Utils

#our library
from CircuitsInformation import CircuitsInformation
import itertools

#our libary
import MESSINetworkBuilder

#For nice-printing of tables
import texttable as tt

############################################################
#Old obsolete comments:
#NICETOHAVE: GUI para ingresar un grafo y que se fije si admite estructura MESSI y haga todas las cuentitas.
#Objetivo para hoy: hacer que calcule los circuitos.
#Objetivo secundario: hacer la otra parte.
#Pregunta: tengo que guardar la info de los circuitos, no? Para calcular los witnesses despues
#Deberia hacer la union "circuital" de vectores

#For running sage:
#./sage -python circuits.py
############################################################




#TODO: maybe we could use sets.
#TODO: names of parameters could be nicer.
def check_if_sigma_subperp_is_mixed(Bperp, Mt, s, d):
    """
    This function checks if SigmaSubperp is mixed.
    TODO: the name of the function is wrong!
    """
    signs = []
    witnesses = []
    
    #Bperp has d rows and s columns.
    #Mt has s-d rows and s columns.
    SList = range(1,s+1)

    #First we get all subsets of size d.
    #TODO: an obvious optimization:
    #Don't store all subsets, but generate them dynamically.
    #Stop generathing them if we get mixed signs for a subset.
    L = Utils.get_r_subsets(SList,d)

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
        Sigma_subperp.append([JcL, mM, JL, mB, result])

    #print_Sigma_subperp_table(Sigma_subperp)
    
    if len(signs) > 1:
        print("Sigma_perp is mixed! Great, then we can find v in T^perp and w in S with same sign.")
        #print witnesses
    else:
        print("Sigma_perp is NOT mixed")
        exit(0)

    #Faltaria armar la tabla, que es bastante grande, no? No se si la voy a armar
    #Pero si voy  devolver

    #for i in range(1, 2**d):



def print_Sigma_subperp_table(table):
    """
    Prints the Sigma_subperp table.
    """
    headings = table.pop(0)

    table2 = []
    map(table2, zip(*table))
    
    table2 = list(map(list, np.transpose(table)))
    
    tab = tt.Texttable()
    tab.header(headings)
    Jc = table2[0]
    menorMt = table2[1]
    J = table2[2]
    menorBPerp = table2[3]
    result = table2[4]

    for row in zip(Jc, menorMt, J, menorBPerp, result):
        tab.add_row(row)
    s = tab.draw()
    print (s)

#Theorem 5.8
def get_multistationarity_witnesses(v, w, s, d):
    x1 = []
    x2 = []
    for i in range (0, s):
        if v[i] != 0:
            x1.append(w[i] * 1.0 / (np.exp(v[i])-1))
        else:
            x1.append(1)

        x2.append(np.exp(v[i]) * x1[i])
    return x1,x2


def get_hardcoded_matrices():
    Bperp = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, -1, 0],
    [0, 1, 0, 0, 2, 2, 1, 1, 2, 1],
    [-0, 0, 1, 0, 1, 1, 1, 1, 1, 1],
    [-0, 0, 0, 1, -1, -1, 0, -0, -1, -1]])

    #Hardcoded from example 2.7
    Mt = np.array([
    [1, 0, 0, 0, 0, 0, -1, 1, 0, -1],
    [0, 1, 0, 0, 0, 0, -1, 1, 0, -1],
    [0, 0, 1, 0, 0, 0, 0, -1, 0, 1],
    [0, 0, 0, 1, 0, 0, 0, -1, 0, 1],
    [0, 0, 0, 0, 1, 0, -1, 1, -1, -1],
    [0, 0, 0, 0, 0, 1, -1, 1, -0, -2]])

    return Bperp, Mt



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


def main(debug):
    """
    gets multistationarity witnesses x^1, x^2, \\kappa or exits.
    the debug flag is used for fast computation
    """
    MESSINetworkBuilder.get_relevant_matrices(False)

    #for easier reading
    Bperp, Mt = get_hardcoded_matrices()

    #TODO: finish this MESSINetworkBuilder function instead of using hardcoded Bperp and Mt.
    #The function should return BPerp and Mt from user input of the MESSI system.
    #MESSINetworkBuilder.get_relevant_matrices(debug)


    #From now on, almost everything is automated

    #Columns of B^\perp
    s = np.shape(Bperp)[1]
    assert s == np.shape(Mt)[1]

    d = np.shape(Bperp)[0]    
    assert np.shape(Mt)[0] == s-d
    assert d<=s

    print("We now compute Sigma_perp and check if it is mixed. Press ENTER to continue.")
    input()
    
    #exits if false
    check_if_sigma_subperp_is_mixed(Bperp, Mt, s, d)

    #circuits of matrices Bperp, Mt
    circuits_information_Bperp = CircuitsInformation(Bperp)
    circuits_information_Mt = CircuitsInformation(Mt)

    if not debug:
        input("Press ENTER to continue.")

    #Steps 2-5
    equal_sign_vectors = get_equal_sign_vectors(s, circuits_information_Bperp, circuits_information_Mt)    

    input("We now compute multistationarity witnesses.")
    if len(equal_sign_vectors) == 0:
        print("No solutions were found. Was the system not s-toric?")
        print("TODO: this should be ckecked automatically")
    else:
        #input("Step 5) Conformal vectors v and w")
        print("For each solution orthanth, we will now produce the steady states x1 and x2, and the reaction constants kappa")

        #We iterate the list equal_sign_vectors, while simultaneously generating counter "index"
        for index, L in enumerate(equal_sign_vectors):

            #Here, %d is replaced by index+1.
            print("Producing witness N° %d" % (index+1))

            first_solution = L
            v = first_solution[0]
            w = first_solution[1]

            #input("6) x^1, x^2, \\kappa")
            x1, x2 = get_multistationarity_witnesses(v, w, s, d)
            print("x1 is %s" % x1)
            print("x2 is %s" % x2)

            #Step 6
            #TODO: there are some hardcoded stuff inside this function
            Utils.get_kappa(x1,x2)
            print("_______________________________________")
            input("Press ENTER to continue.")



if __name__=="__main__":
    debug  = False
    if len(sys.argv) == 2:
        if sys.argv[1] == "debug":
            debug = True
    main(debug)
