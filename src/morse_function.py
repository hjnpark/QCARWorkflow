""" Calculate the pairwise Morse potential between all pairs of atoms in a frame.  Used as an energy function to prevent atom clashes. """
from __future__ import print_function, division
import os, sys
import itertools
import numpy as np
from collections import namedtuple, OrderedDict 
from QCARWorkflow.molecule import Molecule

## The periodic table with atomic masses
AtomicMass = OrderedDict([('H', 1.0079), ('He', 4.0026), 
                          ('Li', 6.941), ('Be', 9.0122), ('B', 10.811), ('C', 12.0107), ('N', 14.0067), ('O', 15.9994), ('F', 18.9984), ('Ne', 20.1797),
                          ('Na', 22.9897), ('Mg', 24.305), ('Al', 26.9815), ('Si', 28.0855), ('P', 30.9738), ('S', 32.065), ('Cl', 35.453), ('Ar', 39.948), 
                          ('K', 39.0983), ('Ca', 40.078), ('Sc', 44.9559), ('Ti', 47.867), ('V', 50.9415), ('Cr', 51.9961), ('Mn', 54.938), ('Fe', 55.845), ('Co', 58.9332), 
                          ('Ni', 58.6934), ('Cu', 63.546), ('Zn', 65.39), ('Ga', 69.723), ('Ge', 72.64), ('As', 74.9216), ('Se', 78.96), ('Br', 79.904), ('Kr', 83.8), 
                          ('Rb', 85.4678), ('Sr', 87.62), ('Y', 88.9059), ('Zr', 91.224), ('Nb', 92.9064), ('Mo', 95.94), ('Tc', 98), ('Ru', 101.07), ('Rh', 102.9055), 
                          ('Pd', 106.42), ('Ag', 107.8682), ('Cd', 112.411), ('In', 114.818), ('Sn', 118.71), ('Sb', 121.76), ('Te', 127.6), ('I', 126.9045), ('Xe', 131.293), 
                          ('Cs', 132.9055), ('Ba', 137.327), ('La', 138.9055), ('Ce', 140.116), ('Pr', 140.9077), ('Nd', 144.24), ('Pm', 145), ('Sm', 150.36), 
                          ('Eu', 151.964), ('Gd', 157.25), ('Tb', 158.9253), ('Dy', 162.5), ('Ho', 164.9303), ('Er', 167.259), ('Tm', 168.9342), ('Yb', 173.04), 
                          ('Lu', 174.967), ('Hf', 178.49), ('Ta', 180.9479), ('W', 183.84), ('Re', 186.207), ('Os', 190.23), ('Ir', 192.217), ('Pt', 195.078), 
                          ('Au', 196.9665), ('Hg', 200.59), ('Tl', 204.3833), ('Pb', 207.2), ('Bi', 208.9804), ('Po', 209), ('At', 210), ('Rn', 222), 
                          ('Fr', 223), ('Ra', 226), ('Ac', 227), ('Th', 232.0381), ('Pa', 231.0359), ('U', 238.0289), ('Np', 237), ('Pu', 244), 
                          ('Am', 243), ('Cm', 247), ('Bk', 247), ('Cf', 251), ('Es', 252), ('Fm', 257), ('Md', 258), ('No', 259), 
                          ('Lr', 262), ('Rf', 261), ('Db', 262), ('Sg', 266), ('Bh', 264), ('Hs', 277), ('Mt', 268)])

## A list of the elements.
Elements = list(AtomicMass.keys())

## Covalent radii from Cordero et al. 'Covalent radii revisited' Dalton Transactions 2008, 2832-2838.
CovalentRadii = OrderedDict([(j, [0.31, 0.28,
                                  1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
                                  1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
                                  2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.61, 1.52, 1.50, 
                                  1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
                                  2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 
                                  1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
                                  2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 
                                  1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87,
                                  1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 
                                  1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50,
                                  2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69][i]) for i, j in list(enumerate(AtomicMass.keys()))[:96]])

## Pairs of known bond energies / lengths (in wavenumbers)
## Reference: Huheey, pps. A-21 to A-34; T.L. Cottrell, "The Strengths of Chemical Bonds," 2nd ed., Butterworths, London, 1958
## B. deB. Darwent, "National Standard Reference Data Series," National Bureau of Standards, No. 31, Washington, DC, 1970; 
## S.W. Benson, J. Chem. Educ., 42, 502 (1965).
## Web Site: http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
BondEnergyRadii = OrderedDict([(('H', 'H'), (432, 74)), (('H', 'B'), (389, 119)), (('H', 'C'), (411, 109)), (('H', 'Si'), (318, 148)), (('H', 'Ge'), (288, 153)), (('H', 'Sn'), (251, 170)), (('H', 'N'), (386, 101)),
                               (('H', 'P'), (322, 144)), (('H', 'As'), (247, 152)), (('H', 'O'), (459, 96)), (('H', 'S'), (363, 134)), (('H', 'Se'), (276, 146)), (('H', 'Te'), (238, 170)), (('H', 'F'), (565, 92)),
                               (('H', 'Cl'), (428, 127)), (('H', 'Br'), (362, 141)), (('H', 'I'), (295, 161)), (('B', 'B'), (293, 0)), (('B', 'O'), (536, 0)), (('B', 'F'), (613, 0)), (('B', 'Cl'), (456, 175)),
                               (('B', 'Br'), (377, 0)), (('C', 'C'), (346, 154)), (('C', 'Si'), (318, 185)), (('C', 'Ge'), (238, 195)), (('C', 'Sn'), (192, 216)), (('C', 'Pb'), (130, 230)), (('C', 'N'), (305, 147)),
                               (('C', 'P'), (264, 184)), (('C', 'O'), (358, 143)), (('C', 'B'), (356, 0)), (('C', 'S'), (272, 182)), (('C', 'F'), (485, 135)), (('C', 'Cl'), (327, 177)), (('C', 'Br'), (285, 194)),
                               (('C', 'I'), (213, 214)), (('Si', 'Si'), (222, 233)), (('Si', 'N'), (355, 0)), (('Si', 'O'), (452, 163)), (('Si', 'S'), (293, 200)), (('Si', 'F'), (565, 160)), (('Si', 'Cl'), (381, 202)),
                               (('Si', 'Br'), (310, 215)), (('Si', 'I'), (234, 243)), (('Ge', 'Ge'), (188, 241)), (('Ge', 'N'), (257, 0)), (('Ge', 'F'), (470, 168)), (('Ge', 'Cl'), (349, 210)), (('Ge', 'Br'), (276, 230)),
                               (('Ge', 'I'), (212, 0)), (('Sn', 'F'), (414, 0)), (('Sn', 'Cl'), (323, 233)), (('Sn', 'Br'), (273, 250)), (('Sn', 'I'), (205, 270)), (('Pb', 'F'), (331, 0)), (('Pb', 'Cl'), (243, 242)),
                               (('Pb', 'Br'), (201, 0)), (('Pb', 'I'), (142, 279)), (('N', 'N'), (167, 145)), (('N', 'O'), (201, 140)), (('N', 'F'), (283, 136)), (('N', 'Cl'), (313, 175)), (('P', 'P'), (201, 221)),
                               (('P', 'O'), (335, 163)), (('P', 'F'), (490, 154)), (('P', 'Cl'), (326, 203)), (('P', 'Br'), (264, 0)), (('P', 'I'), (184, 0)), (('As', 'As'), (146, 243)), (('As', 'O'), (301, 178)),
                               (('As', 'F'), (484, 171)), (('As', 'Cl'), (322, 216)), (('As', 'Br'), (458, 233)), (('As', 'I'), (200, 254)), (('Sb', 'Sb'), (121, 0)), (('Sb', 'F'), (440, 0)), (('Sb', 'Cl'), (248, 0)),
                               (('O', 'O'), (142, 148)), (('O', 'F'), (190, 142)), (('S', 'S'), (226, 205)), (('S', 'F'), (284, 156)), (('S', 'Cl'), (255, 207)), (('Se', 'Se'), (172, 0)), (('F', 'F'), (155, 142)),
                               (('Cl', 'Cl'), (240, 199)), (('Br', 'Br'), (190, 228)), (('I', 'I'), (148, 267)), (('At', 'At'), (116, 0)), (('I', 'O'), (201, 0)), (('I', 'F'), (273, 191)), (('I', 'Cl'), (208, 232)),
                               (('I', 'Br'), (175, 0)), (('Kr', 'F'), (50, 190)), (('Xe', 'O'), (84, 175)), (('Xe', 'F'), (130, 195))])
# Unused parameters (multiple bonds).
#C=C 602 134
#C#C 835 120
#C=N 615 129
#C#N 887 116
#C=O 799 120
#C#O 1072 113
#C=S 573 160
#N=N 418 125
#N#N 942 110
#N=O 607 121
#P=O 544 150
#P=S 335 186
#O=O 494 121
#S=O 522 143
#S=S 425 149
#Se=Se 272 215

## UFF Parameters
## Energy units are in Hartree
## Angles are in radians
## Parameter names: 
# r    = (Bohr) Bond Radius
# th0  = Equil. Angle in Radian
# x    = (Bohr) VdW Distance
# d    = (Hartree) Vdw Well Depth
# zeta = Scale Factor
# z    = Effective Charge
# vsp3 = Parameters for torsional barriers involving a pair of sp3 centers
# vsp2 = Parameters for torsional barriers involving a pair of sp2 centers
# chi  = (Hartree) Electronegativity
# J    = 
# rqeq = 
# qmin = 
# qmax = 
UFFParams = OrderedDict([("h_",    [ 0.668963,  3.14159,  5.45375,  7.01185e-05,  12.0,  7.12000e-01,  0.0,  3.18720e-03,  1.66399e-01,  4.77153e-01,  7.01088e-01, -1.0,  1.0]),
                         ("h_b",   [ 0.869274,  1.45735,  5.45375,  7.01185e-05,  12.0,  7.12000e-01,  0.0,  3.18720e-03,  1.66399e-01,  4.77153e-01,  7.01088e-01, -1.0,  1.0]),
                         ("he",    [  1.60438,  1.57080,  4.46353,  8.92417e-05,  15.2400,  9.80000e-02,  0.0,  3.18720e-03,  3.54995e-01,  1.09659,  2.45664,  0.0,  2.0]),
                         ("li",    [  2.52467,  3.14159,  4.63172,  3.98400e-05,  12.0,  1.026,  0.0,  3.18720e-03,  1.10468e-01,  1.75366e-01,  2.94230, -1.0,  1.0]),
                         ("be",    [  2.02957,  1.91061,  5.18730,  1.35456e-04,  12.0,  1.565,  0.0,  3.18720e-03,  1.79225e-01,  3.26552e-01,  2.34326, -1.0,  2.0]),
                         ("b_2",   [  1.56469,  2.09440,  7.71575,  2.86848e-04,  12.0520,  1.755,  0.0,  3.18720e-03,  1.87787e-01,  3.49116e-01,  1.55335, -5.0,  3.0]),
                         ("b_3",   [  1.58359,  1.91061,  7.71575,  2.86848e-04,  12.0520,  1.755,  0.0,  3.18720e-03,  1.87787e-01,  3.49116e-01,  1.55335, -5.0,  3.0]),
                         ("c_r",   [  1.37761,  2.09440,  7.27733,  1.67328e-04,  12.7300,  1.912,  3.37684e-03,  3.18720e-03,  1.96350e-01,  3.72120e-01,  1.43430, -4.0,  4.0]),
                         ("c_1",   [  1.33415,  3.14159,  7.27733,  1.67328e-04,  12.7300,  1.912,  3.37684e-03,  3.18720e-03,  1.96350e-01,  3.72120e-01,  1.43430, -4.0,  4.0]),
                         ("c_2",   [  1.38328,  2.09440,  7.27733,  1.67328e-04,  12.7300,  1.912,  3.37684e-03,  3.18720e-03,  1.96350e-01,  3.72120e-01,  1.43430, -4.0,  4.0]),
                         ("c_3",   [  1.43052,  1.91061,  7.27733,  1.67328e-04,  12.7300,  1.912,  3.37684e-03,  3.18720e-03,  1.96350e-01,  3.72120e-01,  1.43430, -4.0,  4.0]),
                         ("n_r",   [  1.32092,  2.09440,  6.91640,  1.09958e-04,  13.4070,  2.544,  7.17121e-04,  3.18720e-03,  2.53531e-01,  4.32168e-01,  1.35115, -3.0,  5.0]),
                         ("n_1",   [  1.23966,  3.14159,  6.91640,  1.09958e-04,  13.4070,  2.544,  7.17121e-04,  3.18720e-03,  2.53531e-01,  4.32168e-01,  1.35115, -3.0,  5.0]),
                         ("n_2",   [  1.29446,  1.94081,  6.91640,  1.09958e-04,  13.4070,  2.544,  7.17121e-04,  3.18720e-03,  2.53531e-01,  4.32168e-01,  1.35115, -3.0,  5.0]),
                         ("n_3",   [  1.32281,  1.86227,  6.91640,  1.09958e-04,  13.4070,  2.544,  7.17121e-04,  3.18720e-03,  2.53531e-01,  4.32168e-01,  1.35115, -3.0,  5.0]),
                         ("o_r",   [  1.28501,  1.91986,  6.61404,  9.56161e-05,  14.0850,  2.3,  2.86848e-05,  3.18720e-03,  3.21223e-01,  4.91114e-01,  1.26423, -2.0,  6.0]),
                         ("o_1",   [  1.20753,  3.14159,  6.61404,  9.56161e-05,  14.0850,  2.3,  2.86848e-05,  3.18720e-03,  3.21223e-01,  4.91114e-01,  1.26423, -2.0,  6.0]),
                         ("o_2",   [  1.19809,  2.09440,  6.61404,  9.56161e-05,  14.0850,  2.3,  2.86848e-05,  3.18720e-03,  3.21223e-01,  4.91114e-01,  1.26423, -2.0,  6.0]),
                         ("o_3",   [  1.24344,  1.82404,  6.61404,  9.56161e-05,  14.0850,  2.3,  2.86848e-05,  3.18720e-03,  3.21223e-01,  4.91114e-01,  1.26423, -2.0,  6.0]),
                         ("o_3_z", [ 0.997775,  2.54818,  6.61404,  9.56161e-05,  14.0850,  2.3,  2.86848e-05,  3.18720e-03,  3.21223e-01,  4.91114e-01,  1.26423, -2.0,  6.0]),
                         ("f_",    [  1.26234,  3.14159,  6.35704,  7.96801e-05,  14.7620,  1.735,  0.0,  3.18720e-03,  3.99609e-01,  5.49324e-01,  1.33415, -1.0,  7.0]),
                         ("ne4+4", [  1.73855,  1.57080,  6.12838,  6.69313e-05,  15.4400,  1.94000e-01,  0.0,  3.18720e-03,  4.05709e-01,  7.75404e-01,  3.34104,  0.0,  8.0]),
                         ("na",    [  2.90829,  3.14159,  5.63705,  4.78080e-05,  12.0,  1.081,  0.0,  1.99200e-03,  1.04477e-01,  1.68751e-01,  3.94008, -1.0,  1.0]),
                         ("mg3+2", [  2.68530,  1.91061,  5.70886,  1.76890e-04,  12.0,  1.787,  0.0,  1.99200e-03,  1.45195e-01,  2.71428e-01,  2.83459, -1.0,  2.0]),
                         ("al",    [  2.35082,  1.91061,  8.50188,  8.04769e-04,  11.2780,  1.792,  0.0,  1.99200e-03,  1.49201e-01,  2.63858e-01,  2.26956, -5.0,  3.0]),
                         ("si3",   [  2.11082,  1.91061,  8.11637,  6.40628e-04,  12.1750,  2.323,  1.95216e-03,  1.99200e-03,  1.53170e-01,  2.56288e-01,  2.22232, -4.0,  4.0]),
                         ("p_3+3", [  2.08059,  1.63712,  7.83669,  4.86048e-04,  13.0720,  2.863,  3.82464e-03,  1.99200e-03,  2.00760e-01,  2.93992e-01,  2.08248, -3.0,  5.0]),
                         ("p_3+5", [  1.99555,  1.91061,  7.83669,  4.86048e-04,  13.0720,  2.863,  3.82464e-03,  1.99200e-03,  2.00760e-01,  2.93992e-01,  2.08248, -3.0,  5.0]),
                         ("s_r",   [  2.03523,  1.60919,  7.62504,  4.36647e-04,  13.9690,  2.703,  7.71303e-04,  1.99200e-03,  2.54597e-01,  3.29712e-01,  1.97854, -2.0,  6.0]),
                         ("s_2",   [  1.61383,  2.09440,  7.62504,  4.36647e-04,  13.9690,  2.703,  7.71303e-04,  1.99200e-03,  2.54597e-01,  3.29712e-01,  1.97854, -2.0,  6.0]),
                         ("s_3+2", [  2.01067,  1.60745,  7.62504,  4.36647e-04,  13.9690,  2.703,  7.71303e-04,  1.99200e-03,  2.54597e-01,  3.29712e-01,  1.97854, -2.0,  6.0]),
                         ("s_3+4", [  1.98232,  1.80118,  7.62504,  4.36647e-04,  13.9690,  2.703,  7.71303e-04,  1.99200e-03,  2.54597e-01,  3.29712e-01,  1.97854, -2.0,  6.0]),
                         ("s_3+6", [  1.94075,  1.91061,  7.62504,  4.36647e-04,  13.9690,  2.703,  7.71303e-04,  1.99200e-03,  2.54597e-01,  3.29712e-01,  1.97854, -2.0,  6.0]),
                         ("cl",    [  1.97287,  3.14159,  7.45875,  3.61747e-04,  14.8660,  2.348,  0.0,  1.99200e-03,  3.14718e-01,  3.63521e-01,  1.87839, -1.0,  7.0]),
                         ("ar4+4", [  1.95020,  1.57080,  7.30946,  2.94816e-04,  15.7630,  3.00e-01,  0.0,  1.99200e-03,  3.47829e-01,  4.67080e-01,  3.98354,  0.0,  8.0]),
                         ("k_",    [  3.69063,  3.14159,  7.20364,  5.57760e-05,  12.0,  1.165,  0.0,  1.11552e-03,  8.89693e-02,  1.41116e-01,  4.88683, -1.0,  1.0]),
                         ("ca6+2", [  3.32781,  1.57080,  6.42318,  3.79277e-04,  12.0,  2.141,  0.0,  1.11552e-03,  1.18736e-01,  2.11674e-01,  3.77945, -1.0,  2.0]),
                         ("sc3+3", [  2.85916,  1.91061,  6.22665,  3.02784e-05,  12.0,  2.592,  0.0,  1.11552e-03,  1.24763e-01,  2.26374e-01,  3.30702, -1.0,  3.0]),
                         ("ti3+4", [  2.66829,  1.91061,  5.99988,  2.70912e-05,  12.0,  2.659,  0.0,  1.11552e-03,  1.27519e-01,  2.48423e-01,  3.03679, -1.0,  4.0]),
                         ("ti6+4", [  2.66829,  1.57080,  5.99988,  2.70912e-05,  12.0,  2.659,  0.0,  1.11552e-03,  1.27519e-01,  2.48423e-01,  3.03679, -1.0,  4.0]),
                         ("v_3+5", [  2.64940,  1.91061,  5.94130,  2.54976e-05,  12.0,  2.679,  0.0,  1.11552e-03,  1.34134e-01,  2.50628e-01,  2.77790, -1.0,  5.0]),
                         ("cr6+3", [  2.54168,  1.57080,  5.71264,  2.39040e-05,  12.0,  2.463,  0.0,  1.11552e-03,  1.25498e-01,  2.84070e-01,  2.64940, -1.0,  6.0]),
                         ("mn6+2", [  2.61160,  1.57080,  5.59548,  2.07168e-05,  12.0,  2.43,  0.0,  1.11552e-03,  1.22190e-01,  3.01709e-01,  2.89695, -1.0,  7.0]),
                         ("fe3+2", [  2.39995,  1.91061,  5.50288,  2.07168e-05,  12.0,  2.43,  0.0,  1.11552e-03,  1.38176e-01,  3.04282e-01,  2.63239, -1.0,  8.0]),
                         ("fe6+2", [  2.52278,  1.57080,  5.50288,  2.07168e-05,  12.0,  2.43,  0.0,  1.11552e-03,  1.38176e-01,  3.04282e-01,  2.63239, -1.0,  8.0]),
                         ("co6+3", [  2.34515,  1.57080,  5.42729,  2.23104e-05,  12.0,  2.43,  0.0,  1.11552e-03,  1.50855e-01,  3.06854e-01,  2.65695, -1.0,  9.0]),
                         ("ni4+2", [  2.19964,  1.57080,  5.35548,  2.39040e-05,  12.0,  2.43,  0.0,  1.11552e-03,  1.64084e-01,  3.09059e-01,  2.64184, -1.0,  10.0]),
                         ("cu3+1", [  2.46042,  1.91061,  6.60459,  7.96801e-06,  12.0,  2.43,  0.0,  1.11552e-03,  1.54346e-01,  3.10162e-01,  2.70987, -1.0,  10.0]),
                         ("zn3+2", [  2.25444,  1.91061,  5.22131,  1.97607e-04,  12.0,  1.308,  0.0,  1.11552e-03,  1.87640e-01,  3.14939e-01,  2.64562, -1.0,  2.0]),
                         ("ga3+3", [  2.38105,  1.91061,  8.28267,  6.61345e-04,  11.0,  1.821,  0.0,  1.11552e-03,  1.33803e-01,  2.32254e-01,  2.28846, -5.0,  3.0]),
                         ("ge3",   [  2.26200,  1.91061,  8.08803,  6.03975e-04,  12.0,  2.789,  1.11711e-03,  1.11552e-03,  1.48870e-01,  2.52686e-01,  2.24688, -4.0,  4.0]),
                         ("as3+3", [  2.28846,  1.60745,  7.99354,  4.92423e-04,  13.0,  2.864,  2.39040e-03,  1.11552e-03,  1.90654e-01,  2.79954e-01,  2.27523, -3.0,  5.0]),
                         ("se3+2", [  2.24877,  1.58127,  7.94630,  4.63738e-04,  14.0,  2.764,  5.33856e-04,  1.11552e-03,  2.36223e-01,  3.03620e-01,  2.31302, -2.0,  6.0]),
                         ("br",    [  2.25255,  3.14159,  7.91606,  3.99994e-04,  15.0,  2.519,  0.0,  1.11552e-03,  2.86275e-01,  3.25229e-01,  2.15618, -1.0,  7.0]),
                         ("kr4+4", [  2.16752,  1.57080,  7.82536,  3.50592e-04,  16.0,  4.52000e-01,  0.0,  1.11552e-03,  3.12550e-01,  4.20041e-01,  4.28968,  0.0,  8.0]),
                         ("rb",    [  4.27078,  3.14159,  7.77433,  6.37440e-05,  12.0,  1.592,  0.0,  3.18720e-04,  8.56619e-02,  1.35677e-01,  5.23454, -1.0,  1.0]),
                         ("sr6+2", [  3.87772,  1.57080,  6.88049,  3.74496e-04,  12.0,  2.449,  0.0,  3.18720e-04,  1.11129e-01,  1.79335e-01,  4.56369, -1.0,  2.0]),
                         ("y_3+3", [  3.20875,  1.91061,  6.32113,  1.14739e-04,  12.0,  3.257,  0.0,  3.18720e-04,  1.40749e-01,  2.06529e-01,  3.77567, -1.0,  3.0]),
                         ("zr3+4", [  2.95553,  1.91061,  5.90350,  1.09958e-04,  12.0,  3.667,  0.0,  3.18720e-04,  1.24947e-01,  2.60918e-01,  3.32214, -1.0,  4.0]),
                         ("nb3+5", [  2.78357,  1.91061,  5.98098,  9.40225e-05,  12.0,  3.618,  0.0,  3.18720e-04,  1.30459e-01,  2.48423e-01,  3.02923, -1.0,  5.0]),
                         ("mo3+6", [  2.80435,  1.91061,  5.76744,  8.92417e-05,  12.0,  3.4,  0.0,  3.18720e-04,  1.27335e-01,  2.75985e-01,  2.89128, -1.0,  6.0]),
                         ("mo6+6", [  2.77223,  1.57080,  5.76744,  8.92417e-05,  12.0,  3.4,  0.0,  3.18720e-04,  1.27335e-01,  2.75985e-01,  2.89128, -1.0,  6.0]),
                         ("tc6+5", [  2.49822,  1.57080,  5.66540,  7.64929e-05,  12.0,  3.4,  0.0,  3.18720e-04,  1.20904e-01,  2.93257e-01,  2.83459, -1.0,  7.0]),
                         ("ru6+2", [  2.79302,  1.57080,  5.59926,  8.92417e-05,  12.0,  3.4,  0.0,  3.18720e-04,  1.31378e-01,  2.95095e-01,  2.83459, -1.0,  8.0]),
                         ("rh6+3", [  2.51712,  1.57080,  5.53501,  8.44609e-05,  12.0,  3.508,  0.0,  3.18720e-04,  1.46077e-01,  2.94360e-01,  2.85160, -1.0,  9.0]),
                         ("pd4+2", [  2.52845,  1.57080,  5.47832,  7.64929e-05,  12.0,  3.21,  0.0,  3.18720e-04,  1.58756e-01,  2.93992e-01,  2.91774, -1.0,  10.0]),
                         ("ag1+1", [  2.61916,  3.14159,  5.94886,  5.73696e-05,  12.0,  1.956,  0.0,  3.18720e-04,  1.63019e-01,  2.30343e-01,  3.06514, -1.0,  10.0]),
                         ("cd3+2", [  2.65129,  1.91061,  5.38194,  3.63341e-04,  12.0,  1.65,  0.0,  3.18720e-04,  1.84994e-01,  2.90832e-01,  3.02356, -1.0,  2.0]),
                         ("in3+3", [  2.75711,  1.91061,  8.43385,  9.54567e-04,  11.0,  2.07,  0.0,  3.18720e-04,  1.28842e-01,  2.12850e-01,  2.65318, -5.0,  3.0]),
                         ("sn3",   [  2.64184,  1.91061,  8.29968,  9.03572e-04,  12.0,  2.961,  3.17127e-04,  3.18720e-04,  1.46518e-01,  2.29608e-01,  2.55869, -4.0,  4.0]),
                         ("sb3+3", [  2.65884,  1.59872,  8.35259,  7.15527e-04,  13.0,  2.704,  1.75296e-03,  3.18720e-04,  1.80033e-01,  2.45630e-01,  2.65318, -3.0,  5.0]),
                         ("te3+2", [  2.61916,  1.57516,  8.44708,  6.34253e-04,  14.0,  2.882,  4.78080e-04,  3.18720e-04,  2.13732e-01,  2.59154e-01,  2.60782, -2.0,  6.0]),
                         ("i_",    [  2.61160,  3.14159,  8.50377,  5.40231e-04,  15.0,  2.65,  0.0,  3.18720e-04,  2.50702e-01,  2.76500e-01,  2.519, -1.0,  7.0]),
                         ("xe4+4", [  2.39428,  1.57080,  8.32235,  5.29076e-04,  12.0,  5.56000e-01,  0.0,  3.18720e-04,  2.79109e-01,  3.65653e-01,  4.64684,  0.0,  8.0]),
                         ("cs",    [  4.85660,  3.14159,  8.53589,  7.17121e-05,  12.0,  1.573,  0.0,  1.59360e-04,  8.02231e-02,  1.25755e-01,  5.63894, -1.0,  1.0]),
                         ("ba6+2", [  4.30291,  1.57080,  6.99766,  5.80071e-04,  12.0,  2.727,  0.0,  1.59360e-04,  1.03412e-01,  1.76101e-01,  4.61471, -1.0,  2.0]),
                         ("la3+3", [  3.67174,  1.91061,  6.65561,  2.70912e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.04202e-01,  2.01495e-01,  3.91362, -1.0,  3.0]),
                         ("ce6+3", [  3.47899,  1.57080,  6.71987,  2.07168e-04,  12.0,  3.3,  0.0,  1.59360e-04,  1.01942e-01,  1.97857e-01,  3.63772, -1.0,  3.0]),
                         ("pr6+3", [  3.44497,  1.57080,  6.81435,  1.59360e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.05029e-01,  1.88449e-01,  3.79268, -1.0,  3.0]),
                         ("nd6+3", [  3.43174,  1.57080,  6.75577,  1.59360e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.05415e-01,  1.92602e-01,  3.79268, -1.0,  3.0]),
                         ("pm6+3", [  3.40340,  1.57080,  6.70286,  1.43424e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.05874e-01,  1.96460e-01,  3.77945, -1.0,  3.0]),
                         ("sm6+3", [  3.36371,  1.57080,  6.65184,  1.27488e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.06995e-01,  1.99878e-01,  3.73788, -1.0,  3.0]),
                         ("eu6+3", [  3.34670,  1.57080,  6.60081,  1.27488e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.05782e-01,  2.04876e-01,  4.20842, -1.0,  3.0]),
                         ("gd6+3", [  3.27867,  1.57080,  6.36460,  1.43424e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.16366e-01,  2.18620e-01,  3.71898, -1.0,  3.0]),
                         ("tb6+3", [  3.27301,  1.57080,  6.52144,  1.11552e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.10908e-01,  2.08293e-01,  3.69252, -1.0,  3.0]),
                         ("dy6+3", [  3.23143,  1.57080,  6.47798,  1.11552e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.12287e-01,  2.11050e-01,  3.65473, -1.0,  3.0]),
                         ("ho6+3", [  3.20498,  1.57080,  6.44208,  1.11552e-05,  12.0,  3.416,  0.0,  1.59360e-04,  1.14914e-01,  2.12483e-01,  3.63772, -1.0,  3.0]),
                         ("er6+3", [  3.16151,  1.57080,  6.40806,  1.11552e-05,  12.0,  3.3,  0.0,  1.59360e-04,  1.17101e-01,  2.14210e-01,  3.61883, -1.0,  3.0]),
                         ("tm6+3", [  3.13695,  1.57080,  6.37594,  9.56161e-06,  12.0,  3.3,  0.0,  1.59360e-04,  1.19486e-01,  2.15562e-01,  3.77945, -1.0,  3.0]),
                         ("yb6+3", [  3.09348,  1.57080,  6.34003,  3.63341e-04,  12.0,  2.618,  0.0,  1.59360e-04,  1.20864e-01,  2.17922e-01,  4.07803, -1.0,  3.0]),
                         ("lu6+3", [  3.15773,  1.57080,  6.87860,  6.53377e-05,  12.0,  3.271,  0.0,  1.59360e-04,  1.08884e-01,  1.81018e-01,  3.58292,  1.09470e+02,  3.0]),
                         ("hf3+4", [  3.04435,  1.91061,  5.93563,  1.14739e-04,  12.0,  3.921,  0.0,  1.59360e-04,  1.35971e-01,  2.49893e-01,  3.32403, -1.0,  4.0]),
                         ("ta3+5", [  2.85538,  1.91061,  5.99043,  1.29082e-04,  12.0,  4.075,  0.0,  1.59360e-04,  1.87420e-01,  2.09469e-01,  3.03301, -1.0,  5.0]),
                         ("w_3+4", [  2.88372,  1.91061,  5.79957,  1.06771e-04,  12.0,  3.7,  0.0,  1.59360e-04,  1.70148e-01,  2.43278e-01,  2.90640, -1.0,  6.0]),
                         ("w_3+6", [  2.60782,  1.91061,  5.79957,  1.06771e-04,  12.0,  3.7,  0.0,  1.59360e-04,  1.70148e-01,  2.43278e-01,  2.90640, -1.0,  6.0]),
                         ("w_6+6", [  2.63050,  1.57080,  5.79957,  1.06771e-04,  12.0,  3.7,  0.0,  1.59360e-04,  1.70148e-01,  2.43278e-01,  2.90640, -1.0,  6.0]),
                         ("re3+7", [  2.48310,  1.91061,  5.58225,  1.05178e-04,  12.0,  3.7,  0.0,  1.59360e-04,  1.45526e-01,  2.88112e-01,  3.02356, -1.0,  7.0]),
                         ("re6+5", [  2.59270,  1.57080,  5.58225,  1.05178e-04,  12.0,  3.7,  0.0,  1.59360e-04,  1.45526e-01,  2.88112e-01,  3.02356, -1.0,  7.0]),
                         ("os6+6", [  2.59270,  1.57080,  5.89595,  5.89632e-05,  12.0,  3.7,  0.0,  1.59360e-04,  1.88890e-01,  2.66798e-01,  3.21253, -1.0,  8.0]),
                         ("ir6+3", [  2.59081,  1.57080,  5.36682,  1.16333e-04,  12.0,  3.731,  0.0,  1.59360e-04,  1.83745e-01,  2.93992e-01,  3.52623, -1.0,  9.0]),
                         ("pt4+2", [  2.57759,  1.57080,  5.20431,  1.27488e-04,  12.0,  3.382,  0.0,  1.59360e-04,  1.76028e-01,  3.25596e-01,  2.94230, -1.0,  10.0]),
                         ("au4+3", [  2.38483,  1.57080,  6.22287,  6.21504e-05,  12.0,  2.625,  0.0,  1.59360e-04,  1.79850e-01,  1.90066e-01,  3.05758, -1.0,  10.0]),
                         ("hg1+2", [  2.53223,  3.14159,  5.11171,  6.13536e-04,  12.0,  1.75,  0.0,  1.59360e-04,  2.30416e-01,  3.05752e-01,  3.02356, -1.0,  2.0]),
                         ("tl3+3", [  2.86860,  2.09440,  8.21464,  1.08365e-03,  11.0,  2.068,  0.0,  1.59360e-04,  1.17597e-01,  2.13144e-01,  2.89128, -5.0,  3.0]),
                         ("pb3",   [  2.75711,  1.91061,  8.12015,  1.05656e-03,  12.0,  2.846,  1.59360e-04,  1.59360e-04,  1.43321e-01,  2.59448e-01,  2.72876, -4.0,  4.0]),
                         ("bi3+3", [  2.85727,  1.57080,  8.25810,  8.25485e-04,  13.0,  2.47,  1.59360e-03,  1.59360e-04,  1.72353e-01,  2.74883e-01,  2.86105, -3.0,  5.0]),
                         ("po3+2", [  2.83459,  1.57080,  8.89872,  5.17920e-04,  14.0,  2.33,  4.78080e-04,  1.59360e-04,  1.54713e-01,  3.09427e-01,  2.79679, -2.0,  6.0]),
                         ("at",    [  2.91963,  3.14159,  8.97620,  4.52583e-04,  15.0,  2.24,  0.0,  1.59360e-04,  1.74558e-01,  3.49116e-01,  2.77790, -1.0,  7.0]),
                         ("rn4+4", [  2.68341,  1.57080,  9.00454,  3.95213e-04,  16.0,  5.83000e-01,  0.0,  1.59360e-04,  1.97342e-01,  3.94684e-01,  4.15740,  0.0,  8.0]),
                         ("fr",    [  5.44241,  3.14159,  9.25966,  7.96801e-05,  12.0,  1.847,  0.0,  1.59360e-04,  7.34980e-02,  1.46996e-01,  4.34637, -1.0,  1.0]),
                         ("ra6+2", [  4.74699,  1.57080,  6.94852,  6.43815e-04,  12.0,  2.92,  0.0,  1.59360e-04,  1.04477e-01,  1.78894e-01,  4.15740, -1.0,  2.0]),
                         ("ac6+3", [  3.74733,  1.57080,  6.57247,  5.25888e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.04183e-01,  2.08367e-01,  3.98354, -1.0,  3.0]),
                         ("th6+4", [  3.25222,  1.57080,  6.41751,  4.14336e-05,  12.0,  4.202,  0.0,  1.59360e-04,  1.16678e-01,  2.13512e-01,  3.81347, -1.0,  3.0]),
                         ("pa6+4", [  3.23332,  1.57080,  6.47042,  3.50592e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.09696e-01,  2.13512e-01,  3.40151, -1.0,  3.0]),
                         ("u_6+4", [  3.18230,  1.57080,  6.41562,  3.50592e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.22778e-01,  2.09690e-01,  3.23710, -1.0,  3.0]),
                         ("np6+4", [  3.14828,  1.57080,  6.47042,  3.02784e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.30422e-01,  1.99694e-01,  3.40151, -1.0,  3.0]),
                         ("pu6+4", [  3.13128,  1.57080,  6.47042,  2.54976e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.19177e-01,  2.07191e-01,  3.47710, -1.0,  3.0]),
                         ("am6+4", [  3.13695,  1.57080,  6.38916,  2.23104e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.09861e-01,  2.20751e-01,  3.66985, -1.0,  3.0]),
                         ("cm6+3", [  3.40340,  1.57080,  6.28523,  2.07168e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.04055e-01,  2.34422e-01,  3.59048, -1.0,  3.0]),
                         ("bk6+3", [  3.32781,  1.57080,  6.30980,  2.07168e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.17358e-01,  2.23103e-01,  3.59048, -1.0,  3.0]),
                         ("cf6+3", [  3.30702,  1.57080,  6.26066,  2.07168e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.17487e-01,  2.27917e-01,  3.59048, -1.0,  3.0]),
                         ("es6+3", [  3.25789,  1.57080,  6.23421,  1.91232e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.22484e-01,  2.27035e-01,  3.59048, -1.0,  3.0]),
                         ("fm6+3", [  3.23521,  1.57080,  6.20964,  1.91232e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.24947e-01,  2.27844e-01,  3.59048, -1.0,  3.0]),
                         ("md6+3", [  3.19175,  1.57080,  6.18696,  1.75296e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.27519e-01,  2.28579e-01,  3.59048, -1.0,  3.0]),
                         ("no6+3", [  3.17285,  1.57080,  6.13783,  1.75296e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.27703e-01,  2.33356e-01,  3.59048, -1.0,  3.0]),
                         ("lr6+3", [  3.20875,  1.57080,  6.11515,  1.75296e-05,  12.0,  3.9,  0.0,  1.59360e-04,  1.28622e-01,  2.35194e-01,  3.59048,  1.09470e+02,  3.0])])

UFFSymbols = OrderedDict([(j, ['h_', 'he', 
                               'li', 'be', 'b_3', 'c_3', 'n_3', 'o_3', 'f_', 'ne4+4', 
                               'na', 'mg3+2', 'al', 'si3', 'p_3+3', 's_3+2', 'cl', 'ar4+4', 
                               'k_', 'ca6+2', 'sc3+3', 'ti6+4', 'v_3+5', 'cr6+3', 'mn6+2', 'fe6+2', 'co6+3', 'ni4+2', 'cu3+1', 'zn3+2', 'ga3+3', 'ge3', 'as3+3', 'se3+2', 'br', 'kr4+4', 
                               'rb', 'sr6+2', 'y_3+3', 'zr3+4', 'nb3+5', 'mo6+6', 'tc6+5', 'ru6+2', 'rh6+3', 'pd4+2', 'ag1+1', 'cd3+2', 'in3+3', 'sn3', 'sb3+3', 'te3+2', 'i_', 'xe4+4', 
                               'cs', 'ba6+2', 'la3+3', 'ce6+3', 'pr6+3', 'nd6+3', 'pm6+3', 'sm6+3', 'eu6+3', 'gd6+3', 'tb6+3', 'dy6+3', 'ho6+3', 'er6+3', 'tm6+3', 'yb6+3', 'lu6+3', 
                               'hf3+4', 'ta3+5', 'w_6+6', 're6+5', 'os6+6', 'ir6+3', 'pt4+2', 'au4+3', 'hg1+2', 'tl3+3', 'pb3', 'bi3+3', 'po3+2', 'at', 'rn4+4', 
                               'fr', 'ra6+2', 'ac6+3', 'th6+4', 'pa6+4', 'u_6+4', 'np6+4', 'pu6+4', 'am6+4', 'cm6+3', 'bk6+3', 'cf6+3', 'es6+3', 'fm6+3', 'md6+3', 'no6+3', 'lr6+3'][i]) for i, j in list(enumerate(AtomicMass.keys()))[:103]])

## The following is taken from my old calc-uff-parms.py
## It is intended to give us the force constants for bonds.
auangs = 0.529177
aukcal = 627.5096
auev   = 27.2114
l      = 0.1332
b      = 664.12

## UFF parameters named tuple.
UFFP = namedtuple('UFFP', ['r', 'th0', 'x', 'd', 'zeta', 'z', 'vsp3', 'vsp2', 'chi', 'J', 'rqeq', 'qmin', 'qmax'])

## Build a dictionary of elements to UFF parameters.
## In cases where multiple UFF parameters exist for an element,
## I went for the highest coordination.
UFFElements = OrderedDict()
for s in UFFSymbols:
    d = UFFParams[UFFSymbols[s]]
    A = UFFP(r = d[0] * auangs, th0 = d[1], x = d[2] * auangs,
             d = d[3] * aukcal,
             zeta = d[4],
             z = d[5],
             vsp3 = d[6],
             vsp2 = d[7],
             chi = d[8] * auev,
             J = d[9],
             rqeq = d[10],
             qmin = d[11],
             qmax = d[12])
    UFFElements[s] = A

def uffbond(i,j,bo):
    """
    Return the UFF parameters for a bonded interaction in Gromacs units (nm, kJ mol^-1 nm^-2).
    i = Element symbol (string)
    j = Element symbol (string)
    bo = Bond order (integer)
    """
    ri = UFFElements[i].r
    rj = UFFElements[j].r
    chii = UFFElements[i].chi
    chij = UFFElements[j].chi
    rbo = -1 * l * (ri + rj) * np.log(bo)
    ren = ri*rj*(np.sqrt(chii)-np.sqrt(chij))**2/(chii*ri+chij*rj) #Division here needs to be reviewed
    # print "Symbols = %s %s Base bond length = %.4f Bond Order Correction = %.4f Electronegativity Correction = %.4f" % (i, j, ri+rj, rbo, ren)
    rij = ri + rj + rbo - ren
    kij = b * UFFElements[i].z * UFFElements[j].z / rij**3 #Division here needs to be reviewed 
    # print "rij = %.4e nm kij = %.4e kJ/mol/nm**2" % (rij/10,kij*418.4)
    return rij/10,kij*418.4 #Division here needs to be reviewed 

def uffangle(i,j,k,boij,bojk,theta):
    """
    Return the UFF parameters for an angle interaction in Gromacs units (degrees, kJ mol^-1 rad^-2).
    Not used for the nebterpolator but I decided to keep this code.
    
    i = Element symbol (string)
    j = Element symbol for the middle atom (string)
    k = Element symbol (string)
    boij = Bond order for i-j bond
    bojk = Bond order for j-k bond
    theta = Angle (in radians)
    """
    rij = uffbond(i,j,boij)[0]
    rjk = uffbond(j,k,bojk)[0]
    if abs(theta) < 1.0e-3:
        theta = UFFElements[j].th0
    rik = sqrt(rij**2 + rjk**2 - 2*rij*rjk*np.cos(theta))
    zi = UFFElements[i].z
    zk = UFFElements[k].z
    kijk = b * zi * zk * (3*rij*rjk*(1-np.cos(theta)**2) - rik**2*np.cos(theta)) / rik**5 #Division here needs to be reviewed 
    # print "theta = %.4f nm kij = %.4e kJ/mol/rad**2" % (theta*180./pi,kijk*4.184)
    return theta*180./np.pi, kijk*4.184

def uffLJ(i):
    """
    Return the Lennard-Jones parameters for an element in Gromacs units (nm, kJ/mol)
    
    i = Element symbol (string)
    """
    sigi = (UFFElements[i].x)*(0.5**(1./6))
    epsi = UFFElements[i].d
    # print "UFF uses geometric combination rules"
    # print "sig = %.4e eps = %.4e" % (sigi/10,epsi*4.184)
    return sigi,epsi

MorseParams = OrderedDict()
for i, ei in enumerate(Elements[:96]):
    for j, ej in enumerate(Elements[:96]):
        if 'H' in [ei, ej]:
            r0, k = uffbond(ei, ej, 1)
        else:
            r0, k = uffbond(ei, ej, 2)
        Dij = 0
        # Look up the pair of elements in the list of known bond energies
        if (ei, ej) in BondEnergyRadii:
            Dij = float(BondEnergyRadii[(ei, ej)][0])
        elif (ej, ei) in BondEnergyRadii:
            Dij = float(BondEnergyRadii[(ej, ei)][0])
        # If not present, try to construct the bond energies of the heterobond using homobonds
        elif (ei, ei) in BondEnergyRadii:
            if (ej, ej) in BondEnergyRadii:
                Dij = 0.5*(float(BondEnergyRadii[(ei, ei)][0]) + float(BondEnergyRadii[(ej, ej)][0]))
            else:
                Dij = float(BondEnergyRadii[(ei, ei)][0])
        elif (ej, ej) in BondEnergyRadii:
            Dij = float(BondEnergyRadii[(ej, ej)][0])
        # Finally use the last ditch arbitrary guess of 100.0 kJ/mol
        else:
            Dij = 100.0
        a = np.sqrt(k/(2*Dij)) #Division here needs to be reviewed 
        MorseParams[(ei, ej)] = (r0, Dij, a)
        # if i < 18 and j < 18:
        #     print ei, ej, r0, Dij, a, Dij*(1.0 - np.exp(-a*-r0))**2
            
def PairwiseMorse(M, dzero=True, repulsive=False):
    """
    Calculate the sum of pairwise Morse potentials (and atomic gradients) for all atoms in a Molecule object.
    dzero = Define the zero energy at the dissociation limit (default False)
    repulsive = Purely repulsive interaction (default True)
    """
    apairs, dists = M.distance_matrix()
    R0 = []
    D = []
    A = []
    for i, j in apairs:
        mp = MorseParams[(M.elem[i], M.elem[j])]
        r0 = mp[0]
        d = mp[1]
        a = mp[2]
        if (not repulsive) and (M.elem[i] == 'H' and M.elem[j] == 'H'):
            # Take care of the H-clumping problem
            d *= 0.1
        R0.append(r0)
        D.append(d)
        A.append(a)
    R0 = np.array(R0)
    D = np.array(D)
    A = np.array(A)
    Energies = []
    Gradients = []
    for sn in range(len(M)):
        R = dists[sn]*0.1 # Nanometers
        # Second derivative will be discontinuous if repulsive.
        dar0 = -1.0*A*(R-R0) * (1 if not repulsive else R<R0)
        mexp = np.exp(dar0)
        Gradient = np.zeros((M.na, 3), dtype=float)
        iij = 0
        for ij, rij in zip(apairs, R):
            i = ij[0]
            j = ij[1]
            xij = (M.xyzs[sn][i] - M.xyzs[sn][j]) / 10 #Division here needs to be reviewed 
            grd = 2*D[iij]*(1.0-mexp[iij])*mexp[iij]*A[iij]*xij/rij #Division here needs to be reviewed 
            grd *= (1 if not repulsive else rij < R0[iij])
            Gradient[i, :] += grd
            Gradient[j, :] -= grd
            iij += 1
        # if repulsive:
        #     Energies.append(np.sum(D*mexp**2))
        if (repulsive or not dzero):
            Energies.append(np.sum(D*(1.0-mexp)**2))
        else:
            Energies.append(np.sum(D*((1.0-mexp)**2 - 1.0)))
        Gradients.append(Gradient)
    print(Energies)
    return Energies, Gradients

def PairwiseMorse2(M):
    """
    Calculate the pairwise Morse potentials (and atomic gradients) for all atoms in a Molecule object.
    dzero = Define the zero energy at the dissociation limit (default False)
    repulsive = Purely repulsive interaction (default True)
    """
    apairs, dists, disps = M.distance_displacement()
    R0 = []
    D = []
    A = []
    for i, j in apairs:
        r0, d, a = MorseParams[(M.elem[i], M.elem[j])]
        R0.append(r0)
        D.append(d)
        A.append(a)
    R0 = np.array(R0)
    D = np.array(D)
    A = np.array(A)
    dsqt = D ** 0.5
    Energies = []
    Gradients = []
    for sn in range(len(M)):
        R = dists[sn]*0.1 # Nanometers
        X = disps[sn]*0.1 # Nanometers
        # Second derivative will be discontinuous if repulsive.
        dar0 = -1.0*A*(R-R0) * (R<R0)
        # dar0 = -1.0*A*(R-R0)
        mexp = np.exp(dar0)
        Energies.append(dsqt*(1.0 - mexp))
        Gradient = np.zeros((len(apairs), M.na, 3), dtype=float)
        GrdFac = dsqt*mexp*A
        for iij in range(len(apairs)):
            if R[iij] > R0[iij]: continue
            i = apairs[iij][0]
            j = apairs[iij][1]
            rij = R[iij]
            xij = X[iij] # xj - xi
            grdj = GrdFac[iij]*xij/rij #Division here needs to be reviewed 
            Gradient[iij, i, :] = -grdj
            Gradient[iij, j, :] = +grdj
        Gradients.append(Gradient)

    return Energies, Gradients
