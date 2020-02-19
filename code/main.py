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
import os, os.path




#Custom libraries
from utils import SigmaUtils, CircuitUtils, HardcodedUtils, MatrixUtils, Utils

from messi import MESSINetworkBuilder

def main(debug=False):
    """
    gets multistationarity witnesses x^1, x^2, \\kappa or exits.
    the debug flag is used for fast computation
    """
    messi_network = MESSINetworkBuilder.get_network(debug)


    show_matrices = input('Do you want the program to show debug info, i.e., additional information computed by the program? (YES/NO. Default: NO)\n')

    if 'y' in show_matrices.lower():
        show_matrices = True
    else:
        show_matrices = False
    #print(messi_network.species_names)
    #exit(0)


    #for easier reading
    #Bperp, Mt = HardcodedUtils.get_hardcoded_matrices()

    Mperp = MatrixUtils.build_integer_basis_of_orthogonal_complement_of_stoichiometric_matrix(messi_network)

    #rint('Mperp', Mperp)

    Bt = MatrixUtils.build_binomial_matrix(messi_network).transpose()

    #print('Bt', Bt)


    #print("Mperp", Mperp)
    #print("Bt", Bt)

    #print("species names")
    #print(messi_network.species_names)
    #print("______________________")

    toric_N = MatrixUtils.build_stoichiometric_matrix_from_messi_network(messi_network)

    #print('toric_N', toric_N)

    M = MatrixUtils.build_integer_basis_of_stoichiometric_matrix(messi_network)

    #print('M', M)


    Bperp = MatrixUtils.build_integer_basis_of_orthogonal_complement_of_binomial_matrix(messi_network)
    #print('Bperp', Bperp)

    #print("M - sus filas son una base de S")
    #print(M)

    #print("Stochiometric matrix - las filas generan S")
    stochiometric_matrix2 = MatrixUtils.build_stoichiometric_matrix_from_messi_network(messi_network)

    #print('stoichiometrix matrix', stochiometric_matrix2)

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


    #alternative_Mperp = MatrixUtils.build_alternative_positive_integer_basis_of_ortogonal_complement_of_stoichiometric_matrix(messi_network)
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

        print('Partitions')
        print(messi_network.partitions_names)

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

        print('alternative Mperp')

        #print(alternative_Mperp)

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

    circuits_information_Bperp = CircuitUtils.CircuitsInformation(Bperp)

    #Steps 2-5

    print("Do you prefer to get only 1 multistationarity witness instead of all of them? It's faster this way. (YES/NO. Default: YES)")
    only_one_string = input()


    #only_one = False
    equal_sign_vectors = []
    only_one = False
    if 'y' in only_one_string.lower() or len(only_one_string) == 0:
        only_one = True
        print("Getting only one pair of equal sign vectors...")
        equal_sign_vectors = CircuitUtils.get_only_one_equal_sign_vector(s, circuits_information_Bperp, circuits_information_M)
    else:
        print("Getting all equal sign vectors...")
        equal_sign_vectors = CircuitUtils.get_equal_sign_vectors(s, circuits_information_Bperp, circuits_information_M)

    #print(M)
    #print(Bperp)
    #exit(0)


    print("We now compute multistationarity witnesses...")
    found_kappa = False
    if len(equal_sign_vectors) == 0:
        print("No solutions were found. Was the system not s-toric?")
        print("TODO: this should be ckecked automatically")
    else:

        max_number = 0
        files =  [name for name in os.listdir('./outputs') if os.path.isfile('./outputs/' + name)]
        for filename in files:
            try:
                dot_index = filename.index('.')
                number = int(filename[7:dot_index])
                max_number = max(max_number, number)
            except ValueError:
                continue
        output_filename = 'outputs/output_%d.log' % (max_number + 1)

        print('Saving the output in ' + output_filename + '\n')
        with open(output_filename, 'w') as f:
            #input("Step 5) Conformal vectors v and w")
            if not only_one:
                print("For each solution orthant, we will now produce the steady states x1 and x2, and the reaction constants kappa...")
            else:
                print("We will now produce the steady states x1 and x2, and the reaction constants kappa, for the pair of equal sign vectors found...")
                f.write('(Only 1 pair of steady states was searched for)\n')



            f.write('Species'  + '\n')
            f.write(str(messi_network.species_names) + '\n\n' )
            f.write('Complexes'  + '\n')
            f.write(str(messi_network.complexes_names) + '\n\n' )
            f.write('Partitions'  + '\n')
            f.write(str(messi_network.partitions_names)  + '\n\n')
            f.write('G'  + '\n')
            f.write(str(messi_network.G.edges(data=True))  + '\n\n')
            #f.write('G1', messi_network.G1.edges(data=True))
            #f.write('G1', messi_network.G2.edges(data=True))
            #f.write('G2_circle', messi_network.G2_circle.edges(data=True))


            #We iterate the list equal_sign_vectors, while simultaneously generating counter "index"
            for index, L in enumerate(equal_sign_vectors):

                #Here, %d is replaced by index+1.

                first_solution = L
                v = first_solution[0] #viene de Stoc
                w = first_solution[1] #Viene de binomios

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
                    print('Reaction constants:')
                    Ck = ''
                    for index, val in enumerate(kappa):
                        Ck +='%s: %s | ' % (messi_network.constants_names[index], str(val))
                    Ck = Ck[:-3]
                    print(Ck)


                    print('')
                    print('Total Amounts:')

                    total_amounts_vec = []
                    TA_toprint = ''
                    for relation in messi_network.linear_relations:
                        species_text = ''
                        total_value = 0
                        for species in relation:
                            total_value += x1[species]
                            species_text += messi_network.species_names[species] + ' + '
                        total_amounts_vec.append(total_value)
                        species_text = species_text[:-3]
                        TA_toprint += ('%s = ' % species_text +  "{:.6f}".format(total_value) + '\n')
                        print('%s = ' % species_text, "{:.6f}".format(total_value))

                    all_OK = True
                    for index, relation in enumerate(messi_network.linear_relations):
                        species_text = ''
                        total_value = 0
                        for species in relation:
                            total_value += x2[species]
                            species_text += messi_network.species_names[species] + ' + '
                        species_text = species_text[:-3]

                        if abs(total_value - total_amounts_vec[index]) > 1e-3:
                            print('WARNING', species_text, 'got a different result with x2.')
                            print('Value:', "{:.6f}".format(total_value))
                            all_OK = False
                        #print('%s = ' % species_text, "{:.6f}".format(total_value))
                    if all_OK:
                        print('Total Amounts were calculated with both x1 and x2 and the result did not change.')




                    f.write('Concentrations for x1:'  + '\n')
                    f.write(CX1  + '\n\n')


                    f.write('Concentrations for x2:' + '\n')
                    f.write(CX2 + '\n\n')


                    f.write('Reaction constants:' + '\n')
                    f.write(Ck  + '\n\n')


                    f.write('Total Amounts:'  + '\n')
                    f.write(TA_toprint  + '\n')
                    f.write("_______________________________________\n\n\n")

                    print("_______________________________________")
                input("Press ENTER to continue.\n")



        return found_kappa



if __name__=="__main__":
    debug  = False
    if len(sys.argv) == 2:
        if sys.argv[1] == "debug":
            debug = True
    main(debug)
