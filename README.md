# A MESSI system analyzer

[![Build Status](https://travis-ci.com/billy-mosse/MESSI.png)](https://travis-ci.com/billy-mosse/MESSI)

A simple program that provides rate constants which witness multistationarity of MESSI systems that satisfy conditions stated in Theorem 5.4 of "[The Structure of MESSI Biological Systems](https://arxiv.org/abs/1612.08763)", by Alicia Dickenstein and Mercedes Pérez Millán. (Under construction.)

## What is it useful for?

TODO
<!-- This program does amazing stuff -->

### Installation

TODO

### Example usage (it's really easy!):

#### How to analyze a custom network

TODO

##### Is my system MESSI?

TODO

##### How to test the project

Assuming you have python 3 installed, just run `python tests.py`.

Current (passing) tests:

- test_buildup_of_conformal_circuit_to_matrix_using_Lemma_A5:

Asserts that we can correctly build a circuit conformal to a given matrix.

- test_buildup_of_incidence_matrix_from_network

Asserts that we can correctly build the incidence matrix from a given (MESSINetwork) network.

- test_buildup_of_complexes_matrix_from_network

Asserts that we can correctly build the educt complexes matrix from a given (MESSINetwork) network.

- test_buildup_of_stoichiometric_matrix_from_network

Asserts that we can correctly build the stoichiometric matrix from a given (MESSINetwork) network.

(See [this](https://travis-ci.com/billy-mosse/MESSI).)

#### Examples (built-in systems)

TODO

<!-- ##### How to process a big network -->

<!-- Hello! I'm a comment. I won't appear in the README file in github. In this section we have to write something like "just run python3 main.py and amazing stuff will happen"-->


<!-- Cheatsheet: https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet -->
