from utils import Utils

#TODO: maybe we could use sets.
#TODO: names of parameters could be nicer.
#TODO: change so that we use SIGMA.
def check_if_sigma_subperp_is_mixed(Mperp, Bt, s, d):

    """
    This function checks if SigmaSubperp is mixed.
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

    Sigma_supraperp = [["J", "menor M^\\perp", "Jc", "menor B^t", "Result"]]
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

        mB = Utils.minor(Mperp,JL)
        mM = Utils.minor(Bt, JcL)
        result = int(round(mB * mM * multiplier,0)) #son todos enteros, ver paper.
        if(result > 0): #and '+' not in signs:
            signs.append('+')
            witnesses.append(['+', J, Jc])
        elif(result < 0):# and '-' not in signs:
            signs.append('-')
            witnesses.append(['-', J, Jc])
        Sigma_supraperp.append([JL, mM, JcL, mB, result])

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