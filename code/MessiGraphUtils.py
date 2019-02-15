class Complex:
	def __init__(self, species):
		"""species is a list. For example, ['S0', 'S1']
		"""
		self.species = species



def buildIncidenceMatrix(network):
    return None

def buildComplexesMatrix(network):
    return None

def buildStochiometricMatrix(network):
    return None


class EductVector:
    def __init__(self, complexes):
        self.complexes = complexes

    def doSomeAlgebraStuff(self):
        print("hello")


def buildEductVector(network):
    return None


class MessiNetwork:

	#network es un networkx digraph
    def __init__(self, network):
        """
        intermediates and cores are both lists of (the class) complexes

        This function should build:
        1. Matrix N.
        2. The monomial vector phi(x)?
        3. I, D, Y.
        4. diag(k)? No.
        3. Matrix M. For that I need Mer's input.


        """
        self.incidence_matrix = buildIncidenceMatrix(network)
        self.complexes_matrix = buildComplexesMatrix(network)
        self.stochiometric_matrix = buildStochiometricMatrix(network)


        #what whould the structure of the  educt(origin vector) be?
        #I think it suffices for it to be a list of complexes with some hidden methods...
        self.educt_vector = buildEductVector(network)