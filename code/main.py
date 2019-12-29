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


    show_matrices = input('Do you want the program to show debug info? (YES/NO. Default: NO)')

    if 'Y' in show_matrices:
        show_matrices = True
    else:
        show_matrices = False
    #print(messi_network.species_names)
    #exit(0)


    #for easier reading
    #Bperp, Mt = HardcodedUtils.get_hardcoded_matrices()

    Mperp = MatrixUtils.build_integer_basis_of_orthogonal_complement_of_stoichiometric_matrix(messi_network)

    Bt = MatrixUtils.build_binomial_matrix(messi_network).transpose()


    #print("Mperp", Mperp)
    #print("Bt", Bt)

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


    if show_matrices:
        print('.' * 30)
        print('.' * 30)
        print("Complexes: ")
        print(messi_network.complexes)

        print("Complexes names:")
        print(messi_network.complexes_names)

        print("Species:")
        print(messi_network.species)

        print("Species names:")
        print(messi_network.species_names)

        print("Mperp")
        print(Mperp)

        print("positive M perp")
        print(positive_Mperp)

        print("M")
        print(M)

        print("educt complexes matrix")
        print(educt_complexes_matrix)

        print("Bt transpose")
        print(Bt.transpose())

        print("d", d)

        print("s", s)

        print("stoichiometric matrix", stochiometric_matrix2)

        print("Bperp", Bperp)

        print("toric N", toric_N)

        print('.' * 30)
        print('.' * 30)

    #input("")
    #print("d")
    #print(d)
    #print("s")
    #print(s)


    Mt = M.transpose()

    """s = np.shape(Bperp)[1]
    assert s == np.shape(Mt)[1]

    d = np.shape(Bperp)[0]
    assert np.shape(Mt)[0] == s-d
    assert d<=s"""

    #print("We now compute Sigma_perp and check if it is mixed. Press ENTER to continue.")
    #input()

    #exits if false
    print("Checking if the system is multisatationary by checking Sigma's signs...")
    is_mixed = SigmaUtils.check_if_sigma_subperp_is_mixed(M, Bperp, s, d)
    if not is_mixed:
        return False

    circuits_information_M = CircuitUtils.CircuitsInformation(M)
    #print(circuits_information_M.circuits)

    circuits_information_Bperp = CircuitUtils.CircuitsInformation(Bperp)
    #print(circuits_information_Bperp.circuits)

    #Steps 2-5

    print("Do you prefer to get only 1 multistationarity witness instead of all of them? It's faster this way. (YES/NO. Default: YES)")
    only_one_string = input()
    print("Getting equal sign vectors...")
    #only_one = False
    equal_sign_vectors = []
    if 'y' in only_one_string.lower() or len(only_one_string) == 0:
        equal_sign_vectors = CircuitUtils.get_only_one_equal_sign_vector(s, circuits_information_M, circuits_information_Bperp)
    else:
        equal_sign_vectors = CircuitUtils.get_equal_sign_vectors(s, circuits_information_M, circuits_information_Bperp)

    #print(M)
    #print(Bperp)
    #exit(0)

    print("Mperp", Mperp)

    print("We now compute multistationarity witnesses...")
    found_kappa = False
    if len(equal_sign_vectors) == 0:
        print("No solutions were found. Was the system not s-toric?")
        print("TODO: this should be ckecked automatically")
    else:
        #input("Step 5) Conformal vectors v and w")
        if len(equal_sign_vectors) > 1:
            print("For each solution orthant, we will now produce the steady states x1 and x2, and the reaction constants kappa...")
        else:
            print("We will now produce the steady states x1 and x2, and the reaction constants kappa, for the pair of equal sign vectors found...")

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
            kappa = Utils.get_kappa2(x1,x2, positive_Mperp, educt_complexes_matrix, messi_network, toric_N)

            if isinstance(kappa, list):
                found_kappa = True
                print('Concentrations for x1:')
                CX1 = ''
                for index, val in enumerate(x1):
                    CX1 +='%s: %s | ' % (messi_network.species_names[index], str(val))
                CX1 = CX1[:-3]
                print(CX1)

                print('')
                print('Concentrations for x2:')
                CX2 = ''
                for index, val in enumerate(x2):
                    CX2 +='%s: %s | ' % (messi_network.species_names[index], str(val))
                CX2 = CX2[:-3]
                print(CX2)


                print('')
                print('Reaction constants')
                Ck = ''
                for index, val in enumerate(kappa):
                    print(index, val)
                    Ck +='%s: %s | ' % (messi_network.constants_names[index], str(val))
                Ck = Ck[:-3]
                print(Ck)
                print("_______________________________________")
            input("Press ENTER to continue.")
        return found_kappa



if __name__=="__main__":
    debug  = False
    if len(sys.argv) == 2:
        if sys.argv[1] == "debug":
            debug = True
    main(debug)
