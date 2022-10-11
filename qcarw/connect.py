#! /usr/bin/env python
import time
import copy
from collections import OrderedDict
from .refine import equal

        
def compare_rxns(rxn1, rxn2):
    """
    Checking to see whether two reaction pathways are the same.
 
    parameters
    ----------
    rxn1, rxn2: [Molecule objects]
        Molecule objects reactant, TS, and product in a list.
 
    Return
    ----------
    "True" if they are identical.  
    "False" if they are different.
    """
    L1 = len(rxn1)
    L2 = len(rxn2)
    rxns = {L1: rxn1, L2: rxn2}

    if L1 == L2:
        equal_count = 0
        for M1, M2 in zip(rxn1, rxn2):
            if equal(M1, M2):
                equal_count += 1

        if equal_count == L1:
            return True

        equal_count = 0
        for M1, M2 in zip(rxn1[::-1], rxn2):
            if equal(M1, M2):
                equal_count += 1

        if equal_count == L1:
            return True
    else:
        shorter = rxns[min(rxns)]
        longer = rxns[max(rxns)]
        ite_num = int(abs(L1/3-L2/3) + 1)
        equal_count = 0
        for i in range(ite_num):
            for M1, M2 in zip(shorter, longer[i*3:]):
                if equal(M1, M2):
                    equal_count += 1
                    if equal_count == min(L1, L2):
                        return True

            for M1, M2 in zip(shorter, longer[i*3:][::-1]):
                if equal(M1, M2):
                    equal_count += 1
                    if equal_count == min(L1, L2):
                        return True

    return False

def check_repeat(M):
    if len(M) <= 6:
        return False

    for i in range(len(M)):
        if i%3 == 0:
            for j in range(int(len(M)/3)):            
                if j*3>i:
                    if equal(M[i],M[j*3]) or equal(M[i], M[-1]):
                        return True
    return False
                        


def filterTS(M_info, E1):
    temp = copy.deepcopy(M_info)
    for k, v in M_info.items():
        num = int(len(v)/3)
        Final_E = float(v[-1].qm_energies[0])
        if Final_E > E1:
            continue
        for i in range(num):
            Rct_E = float(v[i*3].qm_energies[0])

            TS_E = float(v[i*3+1].qm_energies[0])
                
            Prd_E = float(v[i*3+2].qm_energies[0])

            if TS_E > Rct_E and TS_E > Prd_E:
                if TS_E > E1:
                    del temp[k]
                    break 
            else:
                print("Transition state has lower energy than reactant/product.")

    return temp

def connect_rxns(M_info, outsiders = None, iteration=0):
    """
    This function will connect unit reactions.
    
    Parameters
    ----------
    M_info : OrderedDict
        OrderedDict with frames in key and Molcule objects in a list [reactant, TS, product]
 
    Return
    ----------
    final : OrderedDict
        Molecule objects consist of unit reactions 
    """
    print("-----------Iteration %i-----------"%iteration)
    rxns = OrderedDict()
    if outsiders is None:
        outsiders = OrderedDict()
    outsiders_temp = copy.deepcopy(M_info)
    filtered = copy.deepcopy(M_info)
    if iteration == 0:
        print("Filtering identical reaction pathways..")
        for i, (k1, v1) in enumerate(M_info.items()):
           for j, (k2, v2) in enumerate(M_info.items()):
               if j > i:
                   if compare_rxns(v1, v2):
                       if len(v1) >= len(v2) and k2 in filtered.keys():
                           del filtered[k2]
                       elif len(v1) < len(v2) and k1 in filtered.keys():
                           del filtered[k1]


        filtered_ratio=1-len(filtered)/len(M_info)
    print("{:.2f}% of unit reactions were filtered ({i} unit reactions are unique out of {i})".format(filtered_ratio*100, len(filtered), len(M_info)))

    connect = 0
    print('Detecting connection points..')
    for i, (k1, v1) in enumerate(filtered.items()):
        for j, (k2, v2) in enumerate(filtered.items()):
            if j > i:
                reac1 = v1[0]
                prod1 = v1[-1]
                reac2 = v2[0]
                prod2 = v2[-1]
                frm = k1 +"/"+ k2
                if equal(reac1, reac2):
                    M = v2[::-1] + v1
                    if check_repeat(M):
                        continue
                    else:
                        try:
                            del outsiders_temp[k1], outsiders_temp[k2]
                        except:
                            pass
                        rxns[frm]=M
                        connect += 1

                elif equal(reac1, prod2):
                    M = v2 + v1
                    if check_repeat(M):
                        continue
                    else:
                        try:
                            del outsiders_temp[k1], outsiders_temp[k2]
                        except:
                            pass
                        rxns[frm] = M
                        connect += 1

                elif equal(prod1, reac2):
                    M = v1 + v2
                    if check_repeat(M):
                        continue
                    else:
                        try:
                            del outsiders_temp[k1], outsiders_temp[k2]
                        except:
                            pass
                        rxns[frm] = M
                        connect += 1

                elif equal(prod1, prod2):
                    M = v1 + v2[::-1]
                    if check_repeat(M):
                        continue
                    else:
                        try:
                            del outsiders_temp[k1], outsiders_temp[k2]
                        except:
                            pass
                        rxns[frm] = M
                        connect += 1
    outsiders.update(outsiders_temp)
    print("Connections", connect)
    print("Number of connected rxns", len(rxns))
    print("Number of outsiders", len(outsiders))

    if connect == 0:
        print("Final filtering of the outsiders")
        filtered = copy.deepcopy(outsiders)
        for i, (k1, v1) in enumerate(outsiders.items()):
            for j, (k2, v2) in enumerate(outsiders.items()):
                if j > i:
                    if compare_rxns(v1, v2):
                        if len(v1) >= len(v2) and k2 in filtered.keys():
                            del filtered[k2]
                        elif len(v1) < len(v2) and k1 in filtered.keys():
                            del filtered[k1]

        print('Done! %i reactions were filtered from %i reactions.' %(len(filtered), len(outsiders)))
        return filtered#final
    else:
        iteration += 1
        return connect_rxns(rxns, outsiders, iteration)


if __name__=='__main__':
    pass    


