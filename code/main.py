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
import itertools

#For nice-printing of tables
import texttable as tt


#Custom libraries
import SigmaUtils
import CircuitUtils
from messi import MESSINetworkBuilder
import Utils

############################################################
#Old obsolete comments:
#NICETOHAVE: GUI para ingresar un grafo y que se fije si admite estructura MESSI y haga todas las cuentitas.
#Objetivo para hoy: hacer que calcule los circuitos.
#Objetivo secundario: hacer la otra parte.
#Pregunta: tengo que guardar la info de los circuitos, no? Para calcular los witnesses despues
#Deberia hacer la union "circuital" de vectores
############################################################

Main function. It gets multistationarity witnesses::
def main(debug):
    """
    gets multistationarity witnesses x^1, x^2, \\kappa or exits.
    the debug flag is used for fast computation
    """
    messi_network = MESSINetworkBuilder.get_network(False)

    #for easier reading
    Bperp, Mt = HardcodedUtils.get_hardcoded_matrices()

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
    SigmaUtils.check_if_sigma_subperp_is_mixed(Bperp, Mt, s, d)

    #circuits of matrices Bperp, Mt
    circuits_information_Bperp = CircuitUtils.CircuitsInformation(Bperp)
    circuits_information_Mt = CircuitUtils.CircuitsInformation(Mt)

    if not debug:
        input("Press ENTER to continue.")

    #Steps 2-5
    equal_sign_vectors = CircuitUtils.get_equal_sign_vectors(s, circuits_information_Bperp, circuits_information_Mt)    

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
            x1https://pythonhosted.org/an_example_pypi_project/sphinx.html, x2 = Utils.get_multistationarity_witnesses(v, w, s, d)
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
