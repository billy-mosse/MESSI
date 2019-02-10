#!/usr/bin sage -python
# -*- coding: utf-8 -*-
import os
os.environ['SAGE_ROOT'] = '/usr/share/sagemath'
os.environ['SAGE_SRC'] = '/usr/share/sagemath' 
os.environ['SAGE_DOC_SRC'] = '/usr/share/sagemath' 
os.environ['SAGE_LOCAL'] = '/usr/share/sagemath' 
os.environ['DOT_SAGE'] = '/usr/share/sagemath' 

from sage.all import *


from scipy import linalg, matrix


#Como instalar cosas nuevas:
#https://ask.sagemath.org/question/35457/importing-python-packages-into-sage-or-vice-versa/
#el que sirve es este: ./sage --python -m easy_install <package_name>

#import sage.all.VectorSpace


#No pude conectar sage con python por un problema de versiones
#asi que tuve que correr algunos comandos de sage a mano
#y usar numpy para las cosas que termine dejando en el programa.

def e(i,n,V):
	#Deje de mirarlo como vector de sage
	return [0]*(i-1) + [1] + [0]*(n-i)


def gs_cofficient(v1, v2):
    return numpy.dot(v2, v1) / numpy.dot(v1, v1)

def multiply(cofficient, v):
    return map((lambda x : x * cofficient), v)

def proj(v1, v2):
    return multiply(gs_cofficient(v1, v2) , v1)

def gs(X):
    Y = []
    for i in range(len(X)):
        temp_vec = X[i]
        for inY in Y :
            proj_vec = proj(inY, X[i])
            #print "i =", i, ", projection vector =", proj_vec
            temp_vec = map(lambda x, y : x - y, temp_vec, proj_vec)
            #print "i =", i, ", temporary vector =", temp_vec
        Y.append(temp_vec)
    return Y


def gram_schmidt_columns(X):
    Q, R = np.linalg.qr(X)
    return Q


n = 10
V = VectorSpace(Reals(),n)
T = V.subspace(
	[e(5, n, V)-e(1, n, V)-e(9, n, V),
	e(6, n, V)-e(2, n, V)-e(10, n, V),
	e(7, n, V)-e(2, n, V)-e(3, n, V),
	e(8, n, V)-e(10, n, V)-e(4, n, V),
	e(9, n, V)+e(1, n, V)-e(10, n, V)-e(2, n, V),
	e(2, n, V)+e(3, n, V)-e(10, n, V)-e(4, n, V)])

Bs = T.basis()

B_ort = Matrix(Bs).right_kernel()
#La matriz B_ort es el generador del ortogonal de T
B_ort = Matrix([[1, 0, 0, 0, 0, 0, 0, 0, -1, 0],
[0, 1, 0, 0, 2, 2, 1, 1, 2, 1],
[-0, 0, 1, 0, 1, 1, 1, 1, 1, 1],
[-0, 0, 0, 1, -1, -1, 0, -0, -1, -1]])



#B = Matrix([v for v in Bs]).transpose()

S = V.subspace(
	[e(1, n, V)+e(9, n, V)-e(5, n, V),
	e(2, n, V)+e(9, n, V)-e(5, n, V),
	e(1, n, V)+e(10, n, V)-e(6, n, V),
	e(1, n, V)+e(10, n, V)-e(6, n, V),
	e(3, n, V)+e(2, n, V)-e(7, n, V),
	e(4, n, V)+e(2, n, V)-e(7, n, V),
	e(4, n, V)+e(10, n, V)-e(8, n, V),
	e(3, n, V)+e(10, n, V)-e(8, n, V)])

Ms = S.basis()
Mt = Matrix([v for v in Ms])


#Computadas a mano con sage, porque no puedo correrlo directamente desde el programa.
#Igual esto esta mal...
B = [[1, 0, 0, 0, 0, 0],
	[-0, 1, 0, 0, 0, 0],
	[-0 -0, 1, 0, 0, 0],
	[-0 -0 -0, 1, 0, 0],
	[0 -0 -0 -0, 1, 0],
	[-1, -1, 1 -0, -1, 0],
	[0 -0, 0 -0, 0, 1],
	[ 0 -0, -1, -1, 0, -1],
	[1 -0, 0 -0, 0, 0],
	[ 0, 1, -1, 1, 0, 0]
]
'''
Bt=[
1 1 0 1 1 1 1 1 0 0;
0 0 1 1 0 0 1 1 0 0;
0 1 0 1 1 1 1 1 1 0;
1 0 1 0 1 1 1 1 0 1
];


B=[
 1  0  0  0  -1   0   0    0  1  0;
 0  1  0  0   0  -1   0    0  0  1;
 0  1  1  0   0   0  -1    0  0  0;
 0  0  0  1   0   0   0   -1  0  1;
-1  1  0  0   0   0   0    0 -1  1;
 0 -1 -1  1   0   0   0    0  0  1
];

'''
M = [
	[1, 0, 0, 0, 0, 0],
	[0, 1, 0, 0, 0, 0],
	[0, 0, 1, 0, 0, 0],
	[0, 0, 0, 1, 0, 0],
	[0, 0, 0, 0, 1, 0],
	[0, 0, 0, 0, 0, 1],
	[-1, -1, 0, 0, -1, -1],
	[1, 1, -1, -1, 1, 1],
	[0, 0, 0, 0, -1 -0],
	[-1, -1, 1, 1, -1, -2]
]

'''def getKappa2(x1, x2):
    M = Matrix([
    [-x1[1] * x9[1], x5[1], 0, 0, 0, 0, 0, 0, 0, -x4[1] * x1[1], -x8[1], 0],
    [0, 0, x5[1], -x2[1] * x1[1],  x6[1], -x3[1]*x2[1], x7[1], x7[1], 0, 0, 0],
    [0, 0, 0, 0, 0, 0, -x3[1] * x2[1], x7[1], 0, 0, 0, -x8[1]],
    [0, 0, 0, 0, 0, 0, 0, 0, x7[1], x4[1]*x1[1], x8[1], 0],
    [x1[1]*x9[1], -x5[1], -x5[1], 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, x2[1]*x1[1], -x6[1], -x6[1], 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, x3[1]*x2[1], -x7[1], -x7[1], 0, 0, 0],
    [0, 0, 0 0, 0, 0, 0, 0, 0, x4[1]*x1[1], -x8[1], x8[1]],
    [-x1[1]*x9[1], x5[1], x5[1], 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -x2[1]*x1[1], x6[1], x6[1], 0, 0, 0, 0, 0, x8[1]]
    ])
    print scipy.linalg.nullspace(M)


def getKappa(x1, x2):
    M = Matrix([
    [-x1[0]*x1[8], x1[4], 0, 0, 0, x1[5], 0, 0, 0, 0, 0, 0], #1
	[0, 0, x1[4], -x1[1]*x1[9], x1[5], -x1[2]*x1[1], x1[6], x1[6], 0, 0, 0], #2
	[0, 0, 0, 0, 0, 0, -x1[2]*x1[1], x1[6], 0, 0, 0, x1[7]], #3
	[0, 0, 0, 0, 0, 0, 0, 0, x1[6], -x1[3]*x1[9], x1[7], 0], #4
	[x1[0]*x1[8], -x1[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #5
	[0, 0, 0, x1[1]*x1[9], -x1[5], -x1[5], 0, 0, 0, 0, 0, 0],#6
	[0, 0, 0, 0, 0, 0, x1[2]*x1[1], -x1[6], -x1[6], 0, 0, 0],#7
	[0, 0, 0, 0, 0, 0, 0, 0, 0, x1[3]*x1[9], -x1[7], -x1[7]], #8
	[-x1[0]*x1[8], x1[4], x1[4], 0, 0, 0, 0, 0, 0, 0, 0, 0], #9
	[0, 0, 0, -x1[1]*x1[9], x1[5], x1[5], 0, 0, 0, -x1[3]*x1[9], x1[7], x1[7]]]#10
    sol = scipy.linalg.nullspace(M)
    k = sol[0]
    print k'''