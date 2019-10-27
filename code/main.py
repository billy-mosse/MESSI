#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Check https://wiki.sagemath.org/Python3-compatible%20code
#if someday we want to use sage again

#Tutorial for documentation using sphinx: https://pythonhosted.org/an_example_pypi_project/sphinx.html

#This shouldn't be necessary, as we are using python 3 now
from __future__ import division, absolute_import, print_function

#For documentation
from pygments import highlight
from pygments.lexers import PythonLexer
from pygments.formatters import HtmlFormatter
#Tutorial! https://codeandchaos.wordpress.com/2012/07/30/sphinx-autodoc-tutorial-for-dummies/

import sys
import os
import numpy as np

from numpy.linalg import det
from numpy import shape
import copy

#For nice-printing of tables
import texttable as tt


#Custom libraries
from utils import SigmaUtils, CircuitUtils, HardcodedUtils, MatrixUtils, Utils

from messi import MESSINetworkBuilder

def main(debug=False):
    """
    gets multistationarity witnesses x^1, x^2, \\kappa or exits.
    the debug flag is used for fast computation
    """
    messi_network = MESSINetworkBuilder.get_network(debug)


    #print(messi_network.species_names)
    #exit(0)

    if debug:
        print("Complexes: ")
        print(messi_network.complexes)

        print("Complexes names:")
        print(messi_network.complexes_names)

        print("Species:")
        print(messi_network.species)

        print("Species names:")
        print(messi_network.species_names)
    #for easier reading
    #Bperp, Mt = HardcodedUtils.get_hardcoded_matrices()

    Mperp = MatrixUtils.build_integer_basis_of_orthogonal_complement_of_stoichiometric_matrix(messi_network)

    Bt = MatrixUtils.build_binomial_matrix(messi_network).transpose()
    #print("species names")
    #print(messi_network.species_names)
    #print("______________________")

    toric_N = MatrixUtils.build_stoichiometric_matrix_from_messi_network(messi_network)

    M = MatrixUtils.build_integer_basis_of_stoichiometric_matrix(messi_network)
    Bperp = MatrixUtils.build_integer_basis_of_orthogonal_complement_of_binomial_matrix(messi_network)
    #print("M - sus filas son una base de S")
    #print(M)

    #print("Stochiometric matrix - las filas generan S")
    stochiometric_matrix2 = MatrixUtils.build_stoichiometric_matrix_from_messi_network(messi_network)

    #print(stochiometric_matrix2.transpose())


    if debug and False:
        print("Bperp")
        print(Bperp)
    #TODO: finish this MESSINetworkBuilder function instead of using hardcoded Bperp and Mt.
    #The function should return BPerp and Mt from user input of the MESSI system.
    #MESSINetworkBuilder.get_relevant_matrices(debug)


    #From now on, almost everything is automated

    #Columns of B^\perp
    #Poner asserts correctos...

    #s = #columnas Bperp
    #s es la cantidad de columnas de Bt
    #d es la cantidad de filas de Mperp

    d = np.shape(Mperp)[0]
    s = np.shape(Mperp)[1]


    #toric M
    positive_Mperp = MatrixUtils.build_positive_integer_basis_of_kernel_of_stoichiometric_matrix(messi_network)
    educt_complexes_matrix = MatrixUtils.build_educt_complexes_matrix(messi_network)
    if debug and False:
        print("positive M perp")
        print(positive_Mperp)
        print("M")
        print(M)
        print("educt complexes matrix")
        print(educt_complexes_matrix)



        print("Mperp")
        print(Mperp)


        print("Bt transpose")
        print(Bt.transpose())
    #input("")
    #print("d")
    #print(d)
    #print("s")
    #print(s)




    """s = np.shape(Bperp)[1]
    assert s == np.shape(Mt)[1]

    d = np.shape(Bperp)[0]
    assert np.shape(Mt)[0] == s-d
    assert d<=s"""

    #print("We now compute Sigma_perp and check if it is mixed. Press ENTER to continue.")
    #input()

    #exits if false
    SigmaUtils.check_if_sigma_subperp_is_mixed(Mperp, Bt, s, d)

    circuits_information_M = CircuitUtils.CircuitsInformation(M)
    #print(circuits_information_M.circuits)

    circuits_information_Bperp = CircuitUtils.CircuitsInformation(Bperp)
    #print(circuits_information_Bperp.circuits)

    if not debug:
        input("Press ENTER to continue.")

    #Steps 2-5

    print("Write YES if you want to fastly have only 1 multistationarity witness (instead of the many)")
    only_one_string = input()
    print("Getting equal sign vectors...")
    only_one = False
    equal_sign_vectors = []
    if only_one_string == 'YES' or 'yes':
        equal_sign_vectors = CircuitUtils.get_only_one_equal_sign_vector(s, circuits_information_M, circuits_information_Bperp)
    else:
        equal_sign_vectors = CircuitUtils.get_equal_sign_vectors(s, circuits_information_M, circuits_information_Bperp)

    #print(M)
    #print(Bperp)
    #exit(0)

    input("We now compute multistationarity witnesses...")

    if len(equal_sign_vectors) == 0:
        print("No solutions were found. Was the system not s-toric?")
        print("TODO: this should be ckecked automatically")
    else:
        #input("Step 5) Conformal vectors v and w")
        if len(equal_sign_vectors) > 1:
            print("For each solution orthanth, we will now produce the steady states x1 and x2, and the reaction constants kappa...")
        else:
            print("We will not produce the steady states x1 and x2, and the reaction constants kappa, for the pair of equal sign vectors found...")

        #We iterate the list equal_sign_vectors, while simultaneously generating counter "index"
        for index, L in enumerate(equal_sign_vectors):

            #Here, %d is replaced by index+1.

            first_solution = L
            w = first_solution[0] #viene de Stoc
            v = first_solution[1] #Viene de binomios

            #input("6) x^1, x^2, \\kappa")


            x1, x2 = Utils.get_multistationarity_witnesses(w, v, s, d)

            #Step 6
            #TODO: there are some hardcoded stuff inside this function
            ret = Utils.get_kappa2(x1,x2, positive_Mperp, educt_complexes_matrix, messi_network, toric_N)

            if ret:
                print("x1 is %s" % x1)
                print("x2 is %s" % x2)
                print("_______________________________________")
            input("Press ENTER to continue.")



if __name__=="__main__":
    debug  = False
    if len(sys.argv) == 2:
        if sys.argv[1] == "debug":
            debug = True
    main(debug)
