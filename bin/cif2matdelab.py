#!/usr/bin/env python
# A Python script to take .CIF files and generate inputs for MatDeLab programs
#
# - We use the CifFile library
#
import re
import os
import itertools
import numpy
#
from math    import ceil, floor, acos, sqrt
import cif2matdelab.structure
import pymatgen

versionstr = '%(prog)s version 0.0\npymatgen version '+pymatgen.__version__

#
# Note on units:
#
# In this script we convert all quantities to Rydberg atomic units
# (rather than Hartree atomic units).
#
rmt_reduce    = 0.97 # scale the muffin tin radius down to avoid overlaps

chemical_name = {} # element name by atomic number
chemical_symb = {} # chemical symbol by atomic number
chemical_chrg = {} # nuclear charge by chemical symbol
element_rmt   = {} # muffin tin radius by atomic number
element_lmb   = {} # maximum l-quantum number for the muffin tin basis by atomic number
                   # incidentally this is also the minimum l-quantum number required to 
                   # describe all occupied orbitals
element_ntle  = {} # by atomic number a list of the number of basis functions for each l-quantum number
element_augm  = {} # by atomic number a list of nature of the basis functions LOC/APW
element_atocc = {} # by atomic number a list of the orbital occupation
element_ptnl  = {} # by atomic number a list of logarithmic derivatives 
element_idmd  = {} # by atomic number a list of core/valence(=0) or unoccupied(=1) status of basis function

chemical_name[1]   = "hydrogen"
chemical_name[2]   = "helium"
chemical_name[3]   = "lithium"
chemical_name[4]   = "berylium"
chemical_name[5]   = "boron"
chemical_name[6]   = "carbon"
chemical_name[7]   = "nitrogen"
chemical_name[8]   = "oxygen"
chemical_name[9]   = "fluorine"
chemical_name[10]  = "neon"
chemical_name[11]  = "sodium"
chemical_name[12]  = "magnesium"
chemical_name[13]  = "aluminium"
chemical_name[14]  = "silicon"
chemical_name[15]  = "phosphorus"
chemical_name[16]  = "sulfur"
chemical_name[17]  = "chlorine"
chemical_name[18]  = "argon"
chemical_name[19]  = "potassium"
chemical_name[20]  = "calcium"
chemical_name[21]  = "scandium"
chemical_name[22]  = "titanium"
chemical_name[23]  = "vanadium"
chemical_name[24]  = "chromium"
chemical_name[25]  = "manganese"
chemical_name[26]  = "iron"
chemical_name[27]  = "cobalt"
chemical_name[28]  = "nickel"
chemical_name[29]  = "copper"
chemical_name[30]  = "zinc"
chemical_name[31]  = "gallium"
chemical_name[32]  = "germanium"
chemical_name[33]  = "arsenic"
chemical_name[34]  = "selenium"
chemical_name[35]  = "bromine"
chemical_name[36]  = "krypton"
chemical_name[37]  = "rubidium"
chemical_name[38]  = "strontium"
chemical_name[39]  = "yttrium"
chemical_name[40]  = "zirconium"
chemical_name[41]  = "niobium"
chemical_name[42]  = "molybdenum"
chemical_name[43]  = "technetium"
chemical_name[44]  = "ruthenium"
chemical_name[45]  = "rhodium"
chemical_name[46]  = "paladium"
chemical_name[47]  = "silver"
chemical_name[48]  = "cadmium"
chemical_name[49]  = "indium"
chemical_name[50]  = "tin"
chemical_name[51]  = "antimony"
chemical_name[52]  = "tellurium"
chemical_name[53]  = "iodine"
chemical_name[54]  = "xenon"
chemical_name[55]  = "cesium"
chemical_name[56]  = "barium"
chemical_name[57]  = "lanthanum"
chemical_name[58]  = "cerium"
chemical_name[59]  = "praseodium"
chemical_name[60]  = "neodium"
chemical_name[61]  = "promethium"
chemical_name[62]  = "samarium"
chemical_name[63]  = "europium"
chemical_name[64]  = "gadolinium"
chemical_name[65]  = "terbium"
chemical_name[66]  = "dysprosium"
chemical_name[67]  = "holmium"
chemical_name[68]  = "erbium"
chemical_name[69]  = "thulium"
chemical_name[70]  = "ytterbium"
chemical_name[71]  = "lutetium"
chemical_name[72]  = "hafnium"
chemical_name[73]  = "tantalum"
chemical_name[74]  = "tungsten"
chemical_name[75]  = "rhenium"
chemical_name[76]  = "osmium"
chemical_name[77]  = "iridium"
chemical_name[78]  = "platinum"
chemical_name[79]  = "gold"
chemical_name[80]  = "mercury"
chemical_name[81]  = "thallium"
chemical_name[82]  = "lead"
chemical_name[83]  = "bismuth"
chemical_name[84]  = "polonium"
chemical_name[85]  = "astatine"
chemical_name[86]  = "radon"
chemical_name[87]  = "francium"
chemical_name[88]  = "radium"
chemical_name[89]  = "actinium"
chemical_name[90]  = "thorium"
chemical_name[91]  = "protactinium"
chemical_name[92]  = "uranium"
chemical_name[93]  = "neptunium"
chemical_name[94]  = "plutonium"
chemical_name[95]  = "americium"
chemical_name[96]  = "curium"
chemical_name[97]  = "berkelium"
chemical_name[98]  = "californium"
chemical_name[99]  = "einsteinium"
chemical_name[100] = "fermium"
chemical_name[101] = "mendelevium"
chemical_name[102] = "nobelium"
chemical_name[103] = "lawrencium"
chemical_name[104] = "rutherfordium"
chemical_name[105] = "dubnium"
chemical_name[106] = "seaborgium"
chemical_name[107] = "bohrium"
chemical_name[108] = "hassium"
chemical_name[109] = "meitnerium"
chemical_name[110] = "darmstadtium"
chemical_name[111] = "rontgenium"
chemical_name[112] = "copernicium"
chemical_name[113] = "nihonium"
chemical_name[114] = "flerovium"
chemical_name[115] = "moscovium"
chemical_name[116] = "livermorium"
chemical_name[117] = "tennessine"
chemical_name[118] = "oganesson"

chemical_symb[1]   = "h"
chemical_symb[2]   = "he"
chemical_symb[3]   = "li"
chemical_symb[4]   = "be"
chemical_symb[5]   = "b"
chemical_symb[6]   = "c"
chemical_symb[7]   = "n"
chemical_symb[8]   = "o"
chemical_symb[9]   = "f"
chemical_symb[10]  = "ne"
chemical_symb[11]  = "na"
chemical_symb[12]  = "mg"
chemical_symb[13]  = "al"
chemical_symb[14]  = "si"
chemical_symb[15]  = "p"
chemical_symb[16]  = "s"
chemical_symb[17]  = "cl"
chemical_symb[18]  = "ar"
chemical_symb[19]  = "k"
chemical_symb[20]  = "ca"
chemical_symb[21]  = "sc"
chemical_symb[22]  = "ti"
chemical_symb[23]  = "v"
chemical_symb[24]  = "cr"
chemical_symb[25]  = "mn"
chemical_symb[26]  = "fe"
chemical_symb[27]  = "co"
chemical_symb[28]  = "ni"
chemical_symb[29]  = "cu"
chemical_symb[30]  = "zn"
chemical_symb[31]  = "ga"
chemical_symb[32]  = "ge"
chemical_symb[33]  = "as"
chemical_symb[34]  = "se"
chemical_symb[35]  = "br"
chemical_symb[36]  = "kr"
chemical_symb[37]  = "rb"
chemical_symb[38]  = "sr"
chemical_symb[39]  = "y"
chemical_symb[40]  = "zr"
chemical_symb[41]  = "nb"
chemical_symb[42]  = "mo"
chemical_symb[43]  = "tc"
chemical_symb[44]  = "ru"
chemical_symb[45]  = "rh"
chemical_symb[46]  = "pd"
chemical_symb[47]  = "ag"
chemical_symb[48]  = "cd"
chemical_symb[49]  = "in"
chemical_symb[50]  = "sn"
chemical_symb[51]  = "sb"
chemical_symb[52]  = "te"
chemical_symb[53]  = "i"
chemical_symb[54]  = "xe"
chemical_symb[55]  = "cs"
chemical_symb[56]  = "ba"
chemical_symb[57]  = "la"
chemical_symb[58]  = "ce"
chemical_symb[59]  = "pr"
chemical_symb[60]  = "nd"
chemical_symb[61]  = "pm"
chemical_symb[62]  = "sm"
chemical_symb[63]  = "eu"
chemical_symb[64]  = "gd"
chemical_symb[65]  = "tb"
chemical_symb[66]  = "dy"
chemical_symb[67]  = "ho"
chemical_symb[68]  = "er"
chemical_symb[69]  = "tm"
chemical_symb[70]  = "yb"
chemical_symb[71]  = "lu"
chemical_symb[72]  = "hf"
chemical_symb[73]  = "ta"
chemical_symb[74]  = "w"
chemical_symb[75]  = "re"
chemical_symb[76]  = "os"
chemical_symb[77]  = "ir"
chemical_symb[78]  = "pt"
chemical_symb[79]  = "au"
chemical_symb[80]  = "hg"
chemical_symb[81]  = "tl"
chemical_symb[82]  = "pb"
chemical_symb[83]  = "bi"
chemical_symb[84]  = "po"
chemical_symb[85]  = "at"
chemical_symb[86]  = "rn"
chemical_symb[87]  = "fr"
chemical_symb[88]  = "ra"
chemical_symb[89]  = "ac"
chemical_symb[90]  = "th"
chemical_symb[91]  = "pa"
chemical_symb[92]  = "u"
chemical_symb[93]  = "np"
chemical_symb[94]  = "pu"
chemical_symb[95]  = "am"
chemical_symb[96]  = "cm"
chemical_symb[97]  = "bk"
chemical_symb[98]  = "cf"
chemical_symb[99]  = "es"
chemical_symb[100] = "fm"
chemical_symb[101] = "md"
chemical_symb[102] = "no"
chemical_symb[103] = "lr"
chemical_symb[104] = "rf"
chemical_symb[105] = "db"
chemical_symb[106] = "sg"
chemical_symb[107] = "bh"
chemical_symb[108] = "hs"
chemical_symb[109] = "mt"
chemical_symb[110] = "ds"
chemical_symb[111] = "rg"
chemical_symb[112] = "cn"
chemical_symb[113] = "nh"
chemical_symb[114] = "fl"
chemical_symb[115] = "mc"
chemical_symb[116] = "lv"
chemical_symb[117] = "ts"
chemical_symb[118] = "og"

for ii in range(1,119):
    chemical_chrg[chemical_symb[ii]] = int(ii)

# Data for element_rmt taken from the species files of Elk 3.1.12
# for elements 1-104. The data for elements 105-118 are set based
# on similarity to the lower atomic number elements.
# The rmt are given in Bohr.
element_rmt[1]   = 1.4
element_rmt[2]   = 1.4
element_rmt[3]   = 1.8
element_rmt[4]   = 1.8
element_rmt[5]   = 1.8
element_rmt[6]   = 1.8
element_rmt[7]   = 1.8
element_rmt[8]   = 1.8
element_rmt[9]   = 2.0
element_rmt[10]  = 1.6
element_rmt[11]  = 2.2
element_rmt[12]  = 2.2
element_rmt[13]  = 2.2
element_rmt[14]  = 2.2
element_rmt[15]  = 2.2
element_rmt[16]  = 2.2
element_rmt[17]  = 2.2
element_rmt[18]  = 2.0
element_rmt[19]  = 2.4
element_rmt[20]  = 2.4
element_rmt[21]  = 2.4
element_rmt[22]  = 2.4
element_rmt[23]  = 2.4
element_rmt[24]  = 2.4
element_rmt[25]  = 2.4
element_rmt[26]  = 2.4
element_rmt[27]  = 2.4
element_rmt[28]  = 2.4
element_rmt[29]  = 2.4
element_rmt[30]  = 2.4
element_rmt[31]  = 2.4
element_rmt[32]  = 2.4
element_rmt[33]  = 2.4
element_rmt[34]  = 2.4
element_rmt[35]  = 2.4
element_rmt[36]  = 2.2
element_rmt[37]  = 2.6
element_rmt[38]  = 2.6
element_rmt[39]  = 2.6
element_rmt[40]  = 2.6
element_rmt[41]  = 2.6
element_rmt[42]  = 2.6
element_rmt[43]  = 2.6
element_rmt[44]  = 2.6
element_rmt[45]  = 2.6
element_rmt[46]  = 2.6
element_rmt[47]  = 2.6
element_rmt[48]  = 2.6
element_rmt[49]  = 2.6
element_rmt[50]  = 2.6
element_rmt[51]  = 2.6
element_rmt[52]  = 2.6
element_rmt[53]  = 2.6
element_rmt[54]  = 2.4
element_rmt[55]  = 2.8
element_rmt[56]  = 2.8
element_rmt[57]  = 2.8
element_rmt[58]  = 2.8
element_rmt[59]  = 2.8
element_rmt[60]  = 2.8
element_rmt[61]  = 2.8
element_rmt[62]  = 2.8
element_rmt[63]  = 2.8
element_rmt[64]  = 2.8
element_rmt[65]  = 2.8
element_rmt[66]  = 2.8
element_rmt[67]  = 2.8
element_rmt[68]  = 2.8
element_rmt[69]  = 2.8
element_rmt[70]  = 2.8
element_rmt[71]  = 2.8
element_rmt[72]  = 2.8
element_rmt[73]  = 2.8
element_rmt[74]  = 2.8
element_rmt[75]  = 2.8
element_rmt[76]  = 2.8
element_rmt[77]  = 2.8
element_rmt[78]  = 2.8
element_rmt[79]  = 2.8
element_rmt[80]  = 2.8
element_rmt[81]  = 2.8
element_rmt[82]  = 2.8
element_rmt[83]  = 2.8
element_rmt[84]  = 2.8
element_rmt[85]  = 2.8
element_rmt[86]  = 2.6
element_rmt[87]  = 3.0
element_rmt[88]  = 3.0
element_rmt[89]  = 3.0
element_rmt[90]  = 3.0
element_rmt[91]  = 3.0
element_rmt[92]  = 3.0
element_rmt[93]  = 3.0
element_rmt[94]  = 3.0
element_rmt[95]  = 3.0
element_rmt[96]  = 3.0
element_rmt[97]  = 3.0
element_rmt[98]  = 3.0
element_rmt[99]  = 3.0
element_rmt[100] = 3.0
element_rmt[101] = 3.0
element_rmt[102] = 3.0
element_rmt[103] = 3.0
element_rmt[104] = 3.0

element_rmt[105] = 3.0
element_rmt[106] = 3.0
element_rmt[107] = 3.0
element_rmt[108] = 3.0
element_rmt[109] = 3.0
element_rmt[110] = 3.0
element_rmt[111] = 3.0
element_rmt[112] = 3.0
element_rmt[113] = 3.0
element_rmt[114] = 3.0
element_rmt[115] = 3.0
element_rmt[116] = 3.0
element_rmt[117] = 3.0
element_rmt[118] = 2.8

# Below we initialize element_lmb, element_ntle, element_augm, element_atocc, element_ptnl, and
# element_idmd. 
# Because these fields for a given element are related it makes sense to group the initialization
# of all these fields by element. 
# Data from WIEN2k *.outputst files

# Rules of thumb:
# 1. Core orbitals are so compact that they are alway local (LOC)
# 2. Non-core s-, and p-orbitals are always so diffuse that they are APW
# 3. F-orbitals are always so compact that they are LOC
# 4. D-orbitals are kind of in between and can be LOC or APW depending on the element
# See also "Introductory Solid State Physics", H.P. Myers, Publisher Taylor & Francis, 1997.
# ISBN: 0-203-21255-X. Pp. 224 and further.
# These rules of thumb have not yet been applied to the table below.

element_lmb[1]   = 0
element_ntle[1]  = [2]
element_augm[1]  = ["APW","LOC"]
element_atocc[1] = [1.0,0.0]
element_ptnl[1]  = [1.8,2.8]
element_idmd[1]  = [0,1]

element_lmb[2]   = 0
element_ntle[2]  = [2]
element_augm[2]  = ["APW","LOC"]
element_atocc[2] = [2.0,0.0]
element_ptnl[2]  = [1.8,2.8]
element_idmd[2]  = [0,1]

element_lmb[3]   = 1
element_ntle[3]  = [3,1]
element_augm[3]  = ["LOC","APW","LOC", "APW"]
element_atocc[3] = [2.0,1.0,0.0, 0.0]
element_ptnl[3]  = [1.8,2.8,3.8, 2.8]
element_idmd[3]  = [0,0,1, 0]

element_lmb[4]   = 1
element_ntle[4]  = [2,1]
element_augm[4]  = ["APW","LOC", "APW"]
element_atocc[4] = [2.0,0.0, 0.0]
element_ptnl[4]  = [2.8,3.8, 2.8]
element_idmd[4]  = [0,1, 0]

element_lmb[5]   = 1
element_ntle[5]  = [2,2]
element_augm[5]  = ["APW","LOC", "APW","LOC"]
element_atocc[5] = [2.0,0.0, 1.0,0.0]
element_ptnl[5]  = [2.8,3.8, 2.8,3.8]
element_idmd[5]  = [0,1, 0,1]

element_lmb[6]   = 1
element_ntle[6]  = [2,2]
element_augm[6]  = ["APW","LOC", "APW","LOC"]
element_atocc[6] = [2.0,0.0, 2.0,0.0]
element_ptnl[6]  = [2.8,3.8, 2.8,3.8]
element_idmd[6]  = [0,1, 0,1]

element_lmb[7]   = 1
element_ntle[7]  = [2,2]
element_augm[7]  = ["APW","LOC", "APW","LOC"]
element_atocc[7] = [2.0,0.0, 3.0,0.0]
element_ptnl[7]  = [2.8,3.8, 2.8,3.8]
element_idmd[7]  = [0,1, 0,1]

element_lmb[8]   = 1
element_ntle[8]  = [2,2]
element_augm[8]  = ["APW","LOC", "APW","LOC"]
element_atocc[8] = [2.0,0.0, 4.0,0.0]
element_ptnl[8]  = [2.8,3.8, 2.8,3.8]
element_idmd[8]  = [0,1, 0,1]

element_lmb[9]   = 1
element_ntle[9]  = [2,2]
element_augm[9]  = ["APW","LOC", "APW","LOC"]
element_atocc[9] = [2.0,0.0, 5.0,0.0]
element_ptnl[9]  = [2.8,3.8, 2.8,3.8]
element_idmd[9]  = [0,1, 0,1]

element_lmb[10]   = 1
element_ntle[10]  = [2,2]
element_augm[10]  = ["APW","LOC", "APW","LOC"]
element_atocc[10] = [2.0,0.0, 6.0,0.0]
element_ptnl[10]  = [2.8,3.8, 2.8,3.8]
element_idmd[10]  = [0,1, 0,1]

element_lmb[11]   = 1
element_ntle[11]  = [3,2]
element_augm[11]  = ["LOC","APW","LOC","APW","LOC"]
element_atocc[11] = [2.0,1.0,0.0,6.0,0.0]
element_ptnl[11]  = [2.8,3.8,4.8,2.8,3.8]
element_idmd[11]  = [0,0,1,0,1]

element_lmb[12]   = 1
element_ntle[12]  = [3,2]
element_augm[12]  = ["LOC","APW","LOC","APW","LOC"]
element_atocc[12] = [2.0,2.0,0.0,6.0,0.0]
element_ptnl[12]  = [2.8,3.8,4.8,2.8,3.8]
element_idmd[12]  = [0,0,1,0,1]

element_lmb[13]   = 1
element_ntle[13]  = [3,3]
element_augm[13]  = ["LOC","APW","LOC","LOC","APW","LOC"]
element_atocc[13] = [2.0,2.0,0.0,6.0,1.0,0.0]
element_ptnl[13]  = [2.8,3.8,4.8,2.8,3.8,4.8]
element_idmd[13]  = [0,0,1,0,0,1]

element_lmb[14]   = 1
element_ntle[14]  = [3,3]
element_augm[14]  = ["LOC","APW","LOC", "LOC","APW","LOC"]
element_atocc[14] = [2.0,2.0,0.0, 6.0,2.0,0.0]
element_ptnl[14]  = [2.8,3.8,4.8, 2.8,3.8,4.8]
element_idmd[14]  = [0,0,1, 0,0,1]

element_lmb[15]   = 1
element_ntle[15]  = [3,3]
element_augm[15]  = ["LOC","APW","LOC", "LOC","APW","LOC"]
element_atocc[15] = [2.0,2.0,0.0, 6.0,3.0,0.0]
element_ptnl[15]  = [2.8,3.8,4.8, 2.8,3.8,4.8]
element_idmd[15]  = [0,0,1, 0,0,1]

element_lmb[16]   = 1
element_ntle[16]  = [3,3]
element_augm[16]  = ["LOC","APW","LOC", "LOC","APW","LOC"]
element_atocc[16] = [2.0,2.0,0.0, 6.0,4.0,0.0]
element_ptnl[16]  = [2.8,3.8,4.8, 2.8,3.8,4.8]
element_idmd[16]  = [0,0,1, 0,0,1]

element_lmb[17]   = 1
element_ntle[17]  = [2,2]
element_augm[17]  = ["APW","LOC","APW","LOC"]
element_atocc[17] = [2.0,0.0,5.0,0.0]
element_ptnl[17]  = [3.8,4.8,3.8,4.8]
element_idmd[17]  = [0,1,0,1]

element_lmb[18]   = 1
element_ntle[18]  = [2,2]
element_augm[18]  = ["APW","LOC","APW","LOC"]
element_atocc[18] = [2.0,0.0,6.0,0.0]
element_ptnl[18]  = [3.8,4.8,3.8,4.8]
element_idmd[18]  = [0,1,0,1]

element_lmb[19]   = 1
element_ntle[19]  = [3,2]
element_augm[19]  = ["LOC","APW","LOC","APW","LOC"]
element_atocc[19] = [2.0,1.0,0.0,6.0,0.0]
element_ptnl[19]  = [3.8,4.8,5.8,3.8,4.8]
element_idmd[19]  = [0,0,1,0,1]

element_lmb[20]   = 1
element_ntle[20]  = [3,2]
element_augm[20]  = ["LOC","APW","LOC","APW","LOC"]
element_atocc[20] = [2.0,2.0,0.0,6.0,0.0]
element_ptnl[20]  = [3.8,4.8,5.8,3.8,4.8]
element_idmd[20]  = [0,0,1,0,1]

element_lmb[21]   = 2
element_ntle[21]  = [3,2,2]
element_augm[21]  = ["LOC","APW","LOC","APW","LOC","APW","LOC"]
element_atocc[21] = [2.0,2.0,0.0,6.0,0.0,1.0,0.0]
element_ptnl[21]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[21]  = [0,0,1,0,1,0,1]

element_lmb[22]   = 2
element_ntle[22]  = [3,2,2]
element_augm[22]  = ["LOC","APW","LOC","APW","LOC","APW","LOC"]
element_atocc[22] = [2.0,2.0,0.0,6.0,0.0,2.0,0.0]
element_ptnl[22]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[22]  = [0,0,1,0,1,0,1]

element_lmb[23]   = 2
element_ntle[23]  = [3,2,2]
element_augm[23]  = ["LOC","APW","LOC","APW","LOC","APW","LOC"]
element_atocc[23] = [2.0,2.0,0.0,6.0,0.0,3.0,0.0]
element_ptnl[23]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[23]  = [0,0,1,0,1,0,1]

element_lmb[24]   = 2
element_ntle[24]  = [3,2,2]
element_augm[24]  = ["LOC","APW","LOC","APW","LOC","APW","LOC"]
element_atocc[24] = [2.0,1.0,0.0,6.0,0.0,5.0,0.0]
element_ptnl[24]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[24]  = [0,0,1,0,1,0,1]

element_lmb[25]   = 2
element_ntle[25]  = [3,2,2]
element_augm[25]  = ["LOC","APW","LOC","LOC","APW","APW","LOC"]
element_atocc[25] = [2.0,2.0,0.0,6.0,0.0,5.0,0.0]
element_ptnl[25]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[25]  = [0,0,1,0,1,0,1]

element_lmb[26]   = 2
element_ntle[26]  = [3,2,2]
element_augm[26]  = ["LOC","APW","LOC","LOC","APW","APW","LOC"]
element_atocc[26] = [2.0,2.0,0.0,6.0,0.0,6.0,0.0]
element_ptnl[26]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[26]  = [0,0,1,0,0,0,1]

element_lmb[27]   = 2
element_ntle[27]  = [3,2,2]
element_augm[27]  = ["LOC","APW","LOC","APW","LOC","APW","LOC"]
element_atocc[27] = [2.0,2.0,0.0,6.0,0.0,7.0,0.0]
element_ptnl[27]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[27]  = [0,0,1,0,1,0,1]

element_lmb[28]   = 2
element_ntle[28]  = [3,2,2]
element_augm[28]  = ["LOC","APW","LOC","APW","LOC","APW","LOC"]
element_atocc[28] = [2.0,2.0,0.0,6.0,0.0,8.0,0.0]
element_ptnl[28]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[28]  = [0,0,1,0,1,0,1]

element_lmb[29]   = 2
element_ntle[29]  = [3,2,2]
element_augm[29]  = ["LOC","APW","LOC","APW","LOC","APW","LOC"]
element_atocc[29] = [2.0,2.0,0.0,6.0,0.0,9.0,0.0]
element_ptnl[29]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[29]  = [0,0,1,0,1,0,1]

element_lmb[30]   = 2
element_ntle[30]  = [3,2,2]
element_augm[30]  = ["LOC","APW","LOC","APW","LOC","APW","LOC"]
element_atocc[30] = [2.0,2.0,0.0,6.0,0.0,10.0,0.0]
element_ptnl[30]  = [3.8,4.8,5.8,3.8,4.8,3.8,4.8]
element_idmd[30]  = [0,0,1,0,1,0,1]

element_lmb[31]   = 2
element_ntle[31]  = [2,2,2]
element_augm[31]  = ["APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[31] = [2.0,0.0, 1.0,0.0, 10.0,0.0]
element_ptnl[31]  = [4.8,5.8, 4.8,5.8, 3.8,4.8]
element_idmd[31]  = [0,1, 0,1, 0,1]

element_lmb[32]   = 2
element_ntle[32]  = [3,3,2]
element_augm[32]  = ["LOC","APW","LOC", "LOC","APW","LOC", "APW","LOC"]
element_atocc[32] = [2.0,2.0,0.0, 6.0,2.0,0.0, 10.0,0.0]
element_ptnl[32]  = [3.8,4.8,5.8, 3.8,4.8,5.8, 3.8,4.8]
element_idmd[32]  = [0,0,1, 0,0,1, 0,1]

element_lmb[33]   = 2
element_ntle[33]  = [3,3,2]
element_augm[33]  = ["LOC","APW","LOC", "LOC","APW","LOC", "LOC","APW"]
element_atocc[33] = [2.0,2.0,0.0, 6.0,3.0,0.0, 10.0,0.0]
element_ptnl[33]  = [3.8,4.8,5.8, 3.8,4.8,5.8, 3.8,4.8]
element_idmd[33]  = [0,0,1, 0,0,1, 0,0]

element_lmb[34]   = 2
element_ntle[34]  = [3,3,2]
element_augm[34]  = ["LOC","APW","LOC", "LOC","APW","LOC", "APW","LOC"]
element_atocc[34] = [2.0,2.0,0.0, 6.0,4.0,0.0, 10.0,0.0]
element_ptnl[34]  = [3.8,4.8,5.8, 3.8,4.8,5.8, 3.8,4.8]
element_idmd[34]  = [0,0,1, 0,0,1, 0,1]

element_lmb[35]   = 2
element_ntle[35]  = [2,2,2]
element_augm[35]  = ["APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[35] = [2.0,0.0, 5.0,0.0, 10.0,0.0]
element_ptnl[35]  = [4.8,5.8, 4.8,5.8, 3.8,4.8]
element_idmd[35]  = [0,1, 0,1, 0,1]

element_lmb[36]   = 2
element_ntle[36]  = [2,2,2]
element_augm[36]  = ["APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[36] = [2.0,0.0, 6.0,0.0, 10.0,0.0]
element_ptnl[36]  = [4.8,5.8, 4.8,5.8, 3.8,4.8]
element_idmd[36]  = [0,1, 0,1, 0,1]

element_lmb[37]   = 2
element_ntle[37]  = [3,2,2]
element_augm[37]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[37] = [2.0,1.0,0.0, 6.0,0.0, 10.0,0.0]
element_ptnl[37]  = [4.8,5.8,6.8, 4.8,5.8, 3.8,4.8]
element_idmd[37]  = [0,0,1, 0,1, 0,1]

element_lmb[38]   = 2
element_ntle[38]  = [3,2,2]
element_augm[38]  = ["LOC","APW","LOC", "LOC","APW", "LOC","APW"]
element_atocc[38] = [2.0,2.0,0.0, 6.0,0.0, 10.0,0.0]
element_ptnl[38]  = [4.8,5.8,6.8, 4.8,5.8, 3.8,4.8]
element_idmd[38]  = [0,0,1, 0,0, 0,0]

element_lmb[39]   = 2
element_ntle[39]  = [3,2,2]
element_augm[39]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[39] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0]
element_ptnl[39]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[39]  = [0,0,1, 0,1, 0,1]

element_lmb[40]   = 2
element_ntle[40]  = [3,2,2]
element_augm[40]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[40] = [2.0,2.0,0.0, 6.0,0.0, 2.0,0.0]
element_ptnl[40]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[40]  = [0,0,1, 0,1, 0,1]

element_lmb[41]   = 2
element_ntle[41]  = [3,2,2]
element_augm[41]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[41] = [2.0,1.0,0.0, 6.0,0.0, 4.0,0.0]
element_ptnl[41]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[41]  = [0,0,1, 0,1, 0,1]

element_lmb[42]   = 2
element_ntle[42]  = [3,2,2]
element_augm[42]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[42] = [2.0,1.0,0.0, 6.0,0.0, 5.0,0.0]
element_ptnl[42]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[42]  = [0,0,1, 0,1, 0,1]

element_lmb[43]   = 2
element_ntle[43]  = [3,2,2]
element_augm[43]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[43] = [2.0,1.0,0.0, 6.0,0.0, 6.0,0.0]
element_ptnl[43]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[43]  = [0,0,1, 0,1, 0,1]

element_lmb[44]   = 2
element_ntle[44]  = [3,2,2]
element_augm[44]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[44] = [2.0,1.0,0.0, 6.0,0.0, 7.0,0.0]
element_ptnl[44]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[44]  = [0,0,1, 0,1, 0,1]

element_lmb[45]   = 2
element_ntle[45]  = [3,2,2]
element_augm[45]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[45] = [2.0,1.0,0.0, 6.0,0.0, 8.0,0.0]
element_ptnl[45]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[45]  = [0,0,1, 0,1, 0,1]

element_lmb[46]   = 2
element_ntle[46]  = [3,2,2]
element_augm[46]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[46] = [2.0,1.0,0.0, 6.0,0.0, 9.0,0.0]
element_ptnl[46]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[46]  = [0,0,1, 0,1, 0,1]

#Maybe the basis below is better for Ag
#element_lmb[47]   = 2
#element_ntle[47]  = [3,2,2]
#element_augm[47]  = ["LOC","LOC","APW", "LOC","APW", "LOC","APW"]
#element_atocc[47] = [2.0,1.0,0.0, 6.0,0.0, 10.0,0.0]
#element_ptnl[47]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
#element_idmd[47]  = [0,0,0, 0,0, 0,0]

element_lmb[47]   = 2
element_ntle[47]  = [3,2,2]
element_augm[47]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[47] = [2.0,1.0,0.0, 6.0,0.0, 10.0,0.0]
element_ptnl[47]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[47]  = [0,0,1, 0,1, 0,1]

element_lmb[48]   = 2
element_ntle[48]  = [3,2,2]
element_augm[48]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[48] = [2.0,2.0,0.0, 6.0,0.0, 10.0,0.0]
element_ptnl[48]  = [4.8,5.8,6.8, 4.8,5.8, 4.8,5.8]
element_idmd[48]  = [0,0,1, 0,1, 0,1]

element_lmb[49]   = 2
element_ntle[49]  = [2,3,2]
element_augm[49]  = ["APW","LOC", "LOC","APW","LOC", "APW","LOC"]
element_atocc[49] = [2.0,0.0, 6.0,1.0,0.0, 10.0,0.0]
element_ptnl[49]  = [5.8,6.8, 4.8,5.8,6.8, 4.8,5.8]
element_idmd[49]  = [0,1, 0,0,1, 0,1]

element_lmb[50]   = 2
element_ntle[50]  = [3,3,2]
element_augm[50]  = ["LOC","APW","LOC", "LOC","APW","LOC", "APW","LOC"]
element_atocc[50] = [2.0,2.0,0.0, 6.0,2.0,0.0, 10.0,0.0]
element_ptnl[50]  = [4.8,5.8,6.8, 4.8,5.8,6.8, 4.8,5.8]
element_idmd[50]  = [0,0,1, 0,0,1, 0,1]

element_lmb[51]   = 2
element_ntle[51]  = [3,3,2]
element_augm[51]  = ["LOC","APW","LOC", "LOC","APW","LOC", "APW","LOC"]
element_atocc[51] = [2.0,2.0,0.0, 6.0,3.0,0.0, 10.0,0.0]
element_ptnl[51]  = [4.8,5.8,6.8, 4.8,5.8,6.8, 4.8,5.8]
element_idmd[51]  = [0,0,1, 0,0,1, 0,1]

element_lmb[52]   = 2
element_ntle[52]  = [3,3,2]
element_augm[52]  = ["LOC","APW","LOC", "LOC","APW","LOC", "APW","LOC"]
element_atocc[52] = [2.0,2.0,0.0, 6.0,4.0,0.0, 10.0,0.0]
element_ptnl[52]  = [4.8,5.8,6.8, 4.8,5.8,6.8, 4.8,5.8]
element_idmd[52]  = [0,0,1, 0,0,1, 0,1]

element_lmb[53]   = 2
element_ntle[53]  = [2,2,2]
element_augm[53]  = ["APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[53] = [2.0,0.0, 5.0,0.0, 10.0,0.0]
element_ptnl[53]  = [5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[53]  = [0,1, 0,1, 0,1]

element_lmb[54]   = 2
element_ntle[54]  = [2,2,2]
element_augm[54]  = ["APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[54] = [2.0,0.0, 6.0,0.0, 10.0,0.0]
element_ptnl[54]  = [5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[54]  = [0,1, 0,1, 0,1]

element_lmb[55]   = 2
element_ntle[55]  = [3,2,2]
element_augm[55]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[55] = [2.0,1.0,0.0, 6.0,0.0, 10.0,0.0]
element_ptnl[55]  = [5.8,6.8,7.8, 5.8,6.8, 4.8,5.8]
element_idmd[55]  = [0,0,1, 0,1, 0,1]

element_lmb[56]   = 2
element_ntle[56]  = [3,2,2]
element_augm[56]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[56] = [2.0,2.0,0.0, 6.0,0.0, 10.0,0.0]
element_ptnl[56]  = [5.8,6.8,7.8, 5.8,6.8, 4.8,5.8]
element_idmd[56]  = [0,0,1, 0,1, 0,1]

element_lmb[57]   = 3
element_ntle[57]  = [3,2,2,1]
element_augm[57]  = ["LOC","APW","LOC", "LOC","APW", "APW","LOC", "APW"]
element_atocc[57] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 0.0]
element_ptnl[57]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8]
element_idmd[57]  = [0,0,1, 0,0, 0,1, 0]

element_lmb[58]   = 3
element_ntle[58]  = [3,2,2,2]
element_augm[58]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[58] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 1.0,0.0]
element_ptnl[58]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[58]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[59]   = 3
element_ntle[59]  = [3,2,2,2]
element_augm[59]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[59] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 2.0,0.0]
element_ptnl[59]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[59]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[60]   = 3
element_ntle[60]  = [3,2,2,2]
element_augm[60]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[60] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 3.0,0.0]
element_ptnl[60]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[60]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[61]   = 3
element_ntle[61]  = [3,2,2,2]
element_augm[61]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[61] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 4.0,0.0]
element_ptnl[61]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[61]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[62]   = 3
element_ntle[62]  = [3,2,2,2]
element_augm[62]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[62] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 5.0,0.0]
element_ptnl[62]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[62]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[63]   = 3
element_ntle[63]  = [3,2,2,2]
element_augm[63]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[63] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 6.0,0.0]
element_ptnl[63]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[63]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[64]   = 3
element_ntle[64]  = [3,2,2,2]
element_augm[64]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[64] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 7.0,0.0]
element_ptnl[64]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[64]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[65]   = 3
element_ntle[65]  = [3,2,2,2]
element_augm[65]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[65] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 8.0,0.0]
element_ptnl[65]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[65]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[66]   = 3
element_ntle[66]  = [3,2,2,2]
element_augm[66]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[66] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 9.0,0.0]
element_ptnl[66]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[66]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[67]   = 3
element_ntle[67]  = [3,2,2,2]
element_augm[67]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[67] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 10.0,0.0]
element_ptnl[67]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[67]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[68]   = 3
element_ntle[68]  = [3,2,2,2]
element_augm[68]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "LOC","APW"]
element_atocc[68] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 11.0,0.0]
element_ptnl[68]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[68]  = [0,0,1, 0,1, 0,1, 0,0]

element_lmb[69]   = 3
element_ntle[69]  = [3,2,2,2]
element_augm[69]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[69] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 12.0,0.0]
element_ptnl[69]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[69]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[70]   = 3
element_ntle[70]  = [3,2,2,2]
element_augm[70]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[70] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 13.0,0.0]
element_ptnl[70]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[70]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[71]   = 3
element_ntle[71]  = [3,2,2,2]
element_augm[71]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[71] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 14.0,0.0]
element_ptnl[71]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[71]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[72]   = 3
element_ntle[72]  = [3,2,2,2]
element_augm[72]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[72] = [2.0,2.0,0.0, 6.0,0.0, 2.0,0.0, 14.0,0.0]
element_ptnl[72]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[72]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[73]   = 3
element_ntle[73]  = [3,2,2,2]
element_augm[73]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[73] = [2.0,2.0,0.0, 6.0,0.0, 3.0,0.0, 14.0,0.0]
element_ptnl[73]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[73]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[74]   = 3
element_ntle[74]  = [3,2,2,2]
element_augm[74]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[74] = [2.0,2.0,0.0, 6.0,0.0, 4.0,0.0, 14.0,0.0]
element_ptnl[74]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[74]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[75]   = 3
element_ntle[75]  = [3,2,2,2]
element_augm[75]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[75] = [2.0,2.0,0.0, 6.0,0.0, 5.0,0.0, 14.0,0.0]
element_ptnl[75]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[75]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[76]   = 3
element_ntle[76]  = [3,2,2,2]
element_augm[76]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[76] = [2.0,2.0,0.0, 6.0,0.0, 6.0,0.0, 14.0,0.0]
element_ptnl[76]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[76]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[77]   = 3
element_ntle[77]  = [3,2,2,2]
element_augm[77]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[77] = [2.0,2.0,0.0, 6.0,0.0, 7.0,0.0, 14.0,0.0]
element_ptnl[77]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[77]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[78]   = 3
element_ntle[78]  = [3,2,2,2]
element_augm[78]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[78] = [2.0,1.0,0.0, 6.0,0.0, 9.0,0.0, 14.0,0.0]
element_ptnl[78]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[78]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[79]   = 3
element_ntle[79]  = [3,2,2,2]
element_augm[79]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[79] = [2.0,1.0,0.0, 6.0,0.0, 10.0,0.0, 14.0,0.0]
element_ptnl[79]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[79]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[80]   = 3
element_ntle[80]  = [3,2,2,2]
element_augm[80]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[80] = [2.0,2.0,0.0, 6.0,0.0, 10.0,0.0, 14.0,0.0]
element_ptnl[80]  = [5.8,6.8,7.8, 5.8,6.8, 5.8,6.8, 4.8,5.8]
element_idmd[80]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[81]   = 3
element_ntle[81]  = [2,2,2,1]
element_augm[81]  = ["APW","LOC", "APW","LOC", "APW","LOC", "APW"]
element_atocc[81] = [2.0,0.0, 1.0,0.0, 10.0,0.0, 0.0]
element_ptnl[81]  = [6.8,7.8, 6.8,7.8, 5.8,6.8, 5.8]
element_idmd[81]  = [0,1, 0,1, 0,1, 0]

element_lmb[82]   = 3
element_ntle[82]  = [2,2,2,1]
element_augm[82]  = ["APW","LOC", "APW","LOC", "APW","LOC", "APW"]
element_atocc[82] = [2.0,0.0, 2.0,0.0, 10.0,0.0, 0.0]
element_ptnl[82]  = [6.8,7.8, 6.8,7.8, 5.8,6.8, 5.8]
element_idmd[82]  = [0,1, 0,1, 0,1, 0]

element_lmb[83]   = 3
element_ntle[83]  = [2,2,2,1]
element_augm[83]  = ["APW","LOC", "APW","LOC", "APW","LOC", "APW"]
element_atocc[83] = [2.0,0.0, 3.0,0.0, 10.0,0.0, 0.0]
element_ptnl[83]  = [6.8,7.8, 6.8,7.8, 5.8,6.8, 5.8]
element_idmd[83]  = [0,1, 0,1, 0,1, 0]

element_lmb[84]   = 3
element_ntle[84]  = [2,2,2,1]
element_augm[84]  = ["APW","LOC", "APW","LOC", "APW","LOC", "APW"]
element_atocc[84] = [2.0,0.0, 4.0,0.0, 10.0,0.0, 0.0]
element_ptnl[84]  = [6.8,7.8, 6.8,7.8, 5.8,6.8, 5.8]
element_idmd[84]  = [0,1, 0,1, 0,1, 0]

element_lmb[85]   = 3
element_ntle[85]  = [2,2,2,1]
element_augm[85]  = ["APW","LOC", "APW","LOC", "APW","LOC", "APW"]
element_atocc[85] = [2.0,0.0, 5.0,0.0, 10.0,0.0, 0.0]
element_ptnl[85]  = [6.8,7.8, 6.8,7.8, 5.8,6.8, 5.8]
element_idmd[85]  = [0,1, 0,1, 0,1, 0]

element_lmb[86]   = 3
element_ntle[86]  = [2,2,2,1]
element_augm[86]  = ["APW","LOC", "APW","LOC", "APW","LOC", "APW"]
element_atocc[86] = [2.0,0.0, 6.0,0.0, 10.0,0.0, 0.0]
element_ptnl[86]  = [6.8,7.8, 6.8,7.8, 5.8,6.8, 5.8]
element_idmd[86]  = [0,1, 0,1, 0,1, 0]

element_lmb[87]   = 3
element_ntle[87]  = [3,2,2,1]
element_augm[87]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW"]
element_atocc[87] = [2.0,1.0,0.0, 6.0,0.0, 10.0,0.0, 0.0]
element_ptnl[87]  = [6.8,7.8,8.8, 6.8,7.8, 5.8,6.8, 5.8]
element_idmd[87]  = [0,0,1, 0,1, 0,1, 0]

element_lmb[88]   = 3
element_ntle[88]  = [3,2,2,1]
element_augm[88]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW"]
element_atocc[88] = [2.0,2.0,0.0, 6.0,0.0, 10.0,0.0, 0.0]
element_ptnl[88]  = [6.8,7.8,8.8, 6.8,7.8, 5.8,6.8, 5.8]
element_idmd[88]  = [0,0,1, 0,1, 0,1, 0]

element_lmb[89]   = 3
element_ntle[89]  = [3,2,2,2]
element_augm[89]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[89] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 0.0,0.0]
element_ptnl[89]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[89]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[90]   = 3
element_ntle[90]  = [3,2,2,2]
element_augm[90]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[90] = [2.0,2.0,0.0, 6.0,0.0, 2.0,0.0, 0.0,0.0]
element_ptnl[90]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[90]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[91]   = 3
element_ntle[91]  = [3,2,2,2]
element_augm[91]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[91] = [2.0,2.0,0.0, 6.0,0.0, 2.0,0.0, 1.0,0.0]
element_ptnl[91]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[91]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[92]   = 3
element_ntle[92]  = [3,2,2,2]
element_augm[92]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[92] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 3.0,0.0]
element_ptnl[92]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[92]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[93]   = 3
element_ntle[93]  = [3,2,2,2]
element_augm[93]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[93] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 4.0,0.0]
element_ptnl[93]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[93]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[94]   = 3
element_ntle[94]  = [3,2,2,2]
element_augm[94]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[94] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 5.0,0.0]
element_ptnl[94]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[94]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[95]   = 3
element_ntle[95]  = [3,2,2,2]
element_augm[95]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[95] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 6.0,0.0]
element_ptnl[95]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[95]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[96]   = 3
element_ntle[96]  = [3,2,2,2]
element_augm[96]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[96] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 7.0,0.0]
element_ptnl[96]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[96]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[97]   = 3
element_ntle[97]  = [3,2,2,2]
element_augm[97]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[97] = [2.0,2.0,0.0, 6.0,0.0, 0.0,0.0, 9.0,0.0]
element_ptnl[97]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[97]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[98]   = 3
element_ntle[98]  = [3,2,2,2]
element_augm[98]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[98] = [2.0,2.0,0.0, 6.0,0.0, 0.0,0.0, 10.0,0.0]
element_ptnl[98]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[98]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[99]   = 3
element_ntle[99]  = [3,2,2,2]
element_augm[99]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[99] = [2.0,2.0,0.0, 6.0,0.0, 0.0,0.0, 11.0,0.0]
element_ptnl[99]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[99]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[100]   = 3
element_ntle[100]  = [3,2,2,2]
element_augm[100]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[100] = [2.0,2.0,0.0, 6.0,0.0, 0.0,0.0, 12.0,0.0]
element_ptnl[100]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[100]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[101]   = 3
element_ntle[101]  = [3,2,2,2]
element_augm[101]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[101] = [2.0,2.0,0.0, 6.0,0.0, 0.0,0.0, 13.0,0.0]
element_ptnl[101]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[101]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[102]   = 3
element_ntle[102]  = [3,2,2,2]
element_augm[102]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[102] = [2.0,2.0,0.0, 6.0,0.0, 0.0,0.0, 14.0,0.0]
element_ptnl[102]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[102]  = [0,0,1, 0,1, 0,1, 0,1]

element_lmb[103]   = 3
element_ntle[103]  = [3,2,2,2]
element_augm[103]  = ["LOC","APW","LOC", "APW","LOC", "APW","LOC", "APW","LOC"]
element_atocc[103] = [2.0,2.0,0.0, 6.0,0.0, 1.0,0.0, 14.0,0.0]
element_ptnl[103]  = [6.8,7.8,8.8, 6.8,7.8, 6.8,7.8, 5.8,6.8]
element_idmd[103]  = [0,0,1, 0,1, 0,1, 0,1]

regname = re.compile("[a-z,A-Z]+")

error = None


def parse_arguments():
    """
    Parse command line arguments for the cif2matdelab.py script when run from
    the command line.
    """
    from argparse import ArgumentParser
    prs = ArgumentParser()
    prs.add_argument("--version",action='version',
                     version=versionstr)
    prs.add_argument("-m","--method",dest="method",help="The method to run, default is qp",
                     default="qp",choices=["dft","hf","gw","qp"])
    prs.add_argument("--code",dest="code",help="The code suite to generate inputs for",
                     default="comsuite",choices=["comsuite","elk","wien2k"])
    prs.add_argument("--cell",dest="cell",help="The kind of unit cell to use",
                     default="primitive",choices=["primitive","conventional"])
    prs.add_argument("--dft-it",dest="dftit",type=int,
                     help="The modified number of DFT iterations")
    prs.add_argument("--hf-it",dest="hfit",type=int,
                     help="The modified number of HF iterations")
    prs.add_argument("--gw-it",dest="gwit",type=int,
                     help="The modified number of GW iterations")
    prs.add_argument("--qp-it",dest="qpit",type=int,
                     help="The modified number of QP iterations")
    prs.add_argument("-k","--Kmax",dest="Kmax",help="The planewave cutoff",
                     default="0")
    prs.add_argument("-r","--restart",help="Restart an old calculation",
                     action="store_true",default=False)
    prs.add_argument('--band',dest="band",default='Y',choices=['Y','N'],
                     help="Generate bandstructure data")
    prs.add_argument('--dos',dest="dos",default='Y',choices=['Y','N'],
                     help="Generate DOS data")
    prs.add_argument('-T','--temperature',dest="temperature",default=1000.0,type=float,
                     help="The temperature (K), the default is 1000.0 K")
    prs.add_argument("ciffile",help="the ciffile of the input structure")
    args = prs.parse_args()
    return args

def cif2float(cifnum):
    """
    Convert a cif-floating point number that may include an uncertainty
    indication to a proper floating point number. 
    In a .cif file the value "0.4254(4)" is a floating point number where
    the digit in brackets gives the uncertainty. To convert this number to
    a regular Python floating point number the uncertainty needs to be 
    eliminated and the resulting string converted.
    """
    ii = cifnum.find("(")
    if ii >= 0:
      pnum = float(cifnum[:ii])
    else:
      pnum = float(cifnum)
    return pnum

def cif2element(label):
    """
    Convert a label for an atom to the corresponding chemical symbol
    of the element. Examples of labels are "Al1", "Al3+", "boron2a",
    "H*251", etc.
    This algorithm could certainly be implemented more efficiently.
    """
    tag = label.lower()
    tag = re.match(regname,tag)
    symbol = None
    if tag:
        tag = tag.string[:tag.end()]
        for ii in range(1,119):
            if tag == chemical_name[ii] or tag == chemical_symb[ii]:
                symbol = chemical_symb[ii]
                break
    if symbol == None:
        error = "Unknown element: "+tag
        print(error)
    return symbol

def translate_elements(ini_struct):
    """
    CIF files may specify elements in funny ways. The structure 
    component just extracts whatever the CIF file contains but this
    may not be suitable for any program. This routine translates the
    elements from the way the CIF file specifies them into chemical
    symbols for the elements.
    The updated dictionary is returned.
    """
    elem_in = ini_struct["symbol"]
    elem_out = []
    for elem in elem_in:
        elem_out.append(cif2element(elem))
    ini_struct["symbol"] = elem_out
    return ini_struct

def establish_method(ini_struct,method):
    """
    Takes the method string and effects the method through suitable
    choices of the different iteration numbers. The different choices
    are effected as follows:
    - dft: iter_dft=40, iter_hf=0 , iter_gw=0 , iter_qp=0
    - hf : iter_dft=40, iter_hf=40, iter_gw=0 , iter_qp=0
    - gw : iter_dft=40, iter_hf=0 , iter_gw=20, iter_qp=0
    - qp : iter_dft=40, iter_hf=0 , iter_gw=0 , iter_qp=20
    """
    if method == "hf":
        ini_struct["iter_dft"] = 40
        ini_struct["iter_hf"]  = 40
        ini_struct["iter_gw"]  = 0
        ini_struct["iter_qp"]  = 0
    elif method == "dft":
        ini_struct["iter_dft"] = 40
        ini_struct["iter_hf"]  = 0
        ini_struct["iter_gw"]  = 0
        ini_struct["iter_qp"]  = 0
    elif method == "gw":
        ini_struct["iter_dft"] = 40
        ini_struct["iter_hf"]  = 0
        ini_struct["iter_gw"]  = 20
        ini_struct["iter_qp"]  = 0
    elif method == "qp":
        ini_struct["iter_dft"] = 40
        ini_struct["iter_hf"]  = 0
        ini_struct["iter_gw"]  = 0
        ini_struct["iter_qp"]  = 20
    else:
        ini_struct["iter_dft"] = 0
        ini_struct["iter_hf"]  = 0
        ini_struct["iter_gw"]  = 0
        ini_struct["iter_qp"]  = 0
        error = "Unknown method: "+method
        print(error)
    return ini_struct

def rkm_fact(ielm):
    """
    Following the Wien2K scheme for the muffin-tin radius adjustment.
    These factors are used to compute the relative size of two atoms.
    I.e. for a pair of atoms of ielm and jelm separated by a distance
    D (and Radius(ielm)+Radius(jelm) > D) then the new radii are 
    computed as

      Radius(ielm) = 0.5*(1+(rkm_fact(ielm)-rkm_fact(jelm))*0.005)*D
      Radius(jelm) = 0.5*(1-(rkm_fact(ielm)-rkm_fact(jelm))*0.005)*D

    This function returns the rkm_fact factors that Wien2K uses.
    See Wien2K setrmt_lapw for details.
    """
    if ielm == 3 or ielm == 13 or ielm == 14:
      # Li, Al, Si
      return 45.0
    elif ielm == 4 or ielm == 5:
      # Be, B
      return 50.0
    elif ielm == 6 or ielm == 15:
      # C, P
      return 55.0
    elif ielm == 7 or ielm == 16:
      # N, S
      return 60.0
    elif (ielm == 8 or (ielm >= 11 and ielm <= 13) or ielm == 17 or 
          ielm == 19 or ielm == 20 or ielm == 37 or ielm == 38 or
          ielm == 55 or ielm == 56):
      # O, Na, Mg, Cl, K, Ca, Rb, Sr, Cs, Ba
      return 65.0
    elif ielm == 9:
      # F
      return 70.0
    elif ((ielm >= 21 and ielm <= 24) or (ielm >= 39 and ielm <= 42) or
          (ielm >= 31 and ielm <= 35)):
      # Sc-Cr, Ga-Br, Y-Mo
      return 75.0
    elif ((ielm >= 25 and ielm <= 30) or (ielm >= 44 and ielm <= 53) or 
          ielm == 57 or ielm == 58 or (ielm >= 72 and ielm <= 75)):
      # Mn-Zn, Ru-I, La, Ce, Hf-Re
      return 80.0
    elif ((ielm >= 59 and ielm <= 71) or (ielm >= 76 and ielm <= 85) or
          (ielm >= 87 and ielm <= 118)):
      # Pr-Lu, Os-At, Fr-Og
      return 85.0
    else:
      return 60.0

def establish_mt_radii(ini_struct):
    """
    Takes the elements of the sites, the radii stored in element_rmt,
    and the distance matrix and generates muffin tin radii that are
    compatible with the current material structure.
    """
    element_rad = {}
    elements = ini_struct["symbol"]
    distmat  = ini_struct["dist_matrix"]
    a        = ini_struct["a"]
    b        = ini_struct["b"]
    c        = ini_struct["c"]
    lena     = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
    lenb     = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])
    lenc     = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])
    minlen   = min(lena,lenb,lenc)
    #
    # Initialize the muffin tin radii taking the lattice vectors into account
    #
    for ii in range(0,len(elements)):
        eli = elements[ii]
        nmi = chemical_chrg[eli]
        rmi = element_rmt[nmi]
        # the muffin tin radius must be smaller than half the shortest lattice vector
        element_rad[nmi] = min(rmi,0.5*minlen)
    #
    # Adjust the muffin tin radii based on the interatomic distances
    #
    for ii in range(0,len(elements)):
        row = distmat[ii]
        eli = elements[ii]
        nmi = chemical_chrg[eli]
        for jj in range(0,len(elements)):
            if ii < jj and row[jj] < 1.0e-6:
                print("ERROR: atoms ",ii+1, " and ",jj+1, " are in the same position!")
                print("Atom ",ii+1,":")
                print("  Element    : ",eli)
                print("  Coordinates: ",ini_struct["xcoord"][ii],
                                        ini_struct["ycoord"][ii],
                                        ini_struct["zcoord"][ii])
                print("Atom ",jj+1,":")
                print("  Element    : ",elj)
                print("  Coordinates: ",ini_struct["xcoord"][jj],
                                        ini_struct["ycoord"][jj],
                                        ini_struct["zcoord"][jj])
            if ii != jj and row[jj] > 0.0:
                rmi = element_rmt[nmi]
                rr  = row[jj]
                elj = elements[jj]
                nmj = chemical_chrg[elj]
                rmj = element_rmt[nmj]
                if rmi+rmj > rr:
                    fi = rkm_fact(nmi)
                    fj = rkm_fact(nmj)
                    rmi = 0.5*(1.0+(fi-fj)*0.005)*rr
                    rmj = 0.5*(1.0-(fi-fj)*0.005)*rr
                    scale = rr/(rmi+rmj)
                if nmi in element_rad:
                    element_rad[nmi] = min(element_rad[nmi],rmi)
                else:
                    element_rad[nmi] = rmi
                if nmj in element_rad:
                    element_rad[nmj] = min(element_rad[nmj],rmj)
                else:
                    element_rad[nmj] = rmj
    ini_struct["element_rad"] = element_rad
    return ini_struct

def volume_sphere(radius):
    """
    Return the volume of a sphere given its radius

    The volume of a sphere is 4/3*pi*(r^3)
    """
    pi = acos(-1.0)
    vol = 4.0/3.0*pi*(radius**3)
    return vol

def establish_atoms_volume(ini_struct):
    """
    Calculate the total volume of the atoms within the muffin tin radii

    For all atoms in the unit cell calculate the volume of each atom
    and add all volumes up. Add the total volume to the ini_struct 
    dictionary.
    """
    elements = ini_struct["symbol"]
    radii    = ini_struct["element_rad"]
    atoms_volume = 0.0
    for element in elements:
        number = chemical_chrg[element]
        radius = radii[number]
        atoms_volume += volume_sphere(radius)
    ini_struct["atoms_volume"] = atoms_volume
    return ini_struct

def establish_Kmin(ini_struct):
    """
    Establish the minimum plane wave cut-off

    From the relation RK=l [1] where R is the muffin tin radius, K
    is the plane wave cut-off and l is the highest angular momentum
    of the occupied orbitals of an atom we can find K=l/R and use
    that to establish the smallest plane wave cut-off that still
    provides a qualitatively correct description of the system.

    In addition Klim is established which derives from the fact that 
    the GW code currently cannot handle l quantum numbers larger than
    10. Hence Klim is the minimum value of 10/R across all elements.

    The values of Kmin and Klim are returned in the ini_struct
    dictionary.

    #References#

    [1] D.J. Singh, "Planewaves, pseudo-potentials and the LAPW method",
        Springer (1994), ISBN: 978-1-4757-2314-4, 
        DOI: 10.1007/978-1-4757-2312-0, pp. 62-63.
    """
    element_rad = ini_struct["element_rad"]
    elements    = ini_struct["symbol"]
    Kmin        = 0.00
    Klim        = 1000.0
    for ii in range(0,len(elements)):
        eli  = elements[ii]
        nmi  = chemical_chrg[eli]
        ll   = 2*element_lmb[nmi]
        rr   = element_rad[nmi]
        Kmin = max(Kmin,ll/rr)
        Klim = min(Klim,10.0/rr) # the GW code can handle max l upto 10
    Kmin = int(ceil(Kmin))
    Klim = int(floor(Klim))
    ini_struct["Kmin"] = Kmin
    ini_struct["Klim"] = Klim
    return ini_struct

def establish_Kmax(ini_struct,Kmax=0):
    """
    Establish the maximum plane wave cut off 

    The maximum plane wave cut off must satisfy the cut off requirement
    for the worst calculation that is in principle correct. It should
    also be small enough that the maximum l-quantum number never 
    exceeds 10 as is encapsulated in Klim. Otherwise it should adopt a
    user specified cut off.

    The value of Kmax is returned in the ini_struct dictionary.
    """
    Kmin = ini_struct["Kmin"]
    Klim = ini_struct["Klim"]
    Kfin = min(max(Kmax,Kmin,1),Klim)
    if Klim < Kmin:
        print("Trouble with the Muffin Tin radii:")
        print("Kmin = %d Kmax = %d" %(Kmin,Klim))
        print("No satisfactory solution available!")
        print("Generating a low quality but working input instead!")
    ini_struct["Kmax"] = Kfin
    return ini_struct

def establish_r_grid(ini_struct):
    """
    Establish the real space grid

    Given the K_max and the lattice vectors work out how many grid
    points are needed in each dimension.

    The number of grid points in all dimensions is returned in
    ini_struct under key "nrdiv". Also returned is "mdiv" which
    is about 4/3*nrdiv.

    In addition this function takes the packing factor of the
    material into account. For packing factors larger than 0.5
    the regular approach to calculating the grid sizes is sufficient.
    For packing factors smaller than 0.3 about twice as many plane
    waves are needed to represent the interstitial region well. 
    Between 0.3 and 0.5 the grid sizes are scaled linearly with
    a factor 2 to 1.
    """
    Vcel  = ini_struct["cell_volume"]
    Vatm  = ini_struct["atoms_volume"]
    SpGr  = ini_struct["spacegroup"]
    pfac  = Vatm/Vcel  # packing factor
    scale = 1.0
    if pfac < 0.3:
        scale = 2.0
    elif pfac < 0.5:
        scale = 1.0+(0.5-pfac)*5
    pi    = acos(-1.0)
    r43   = 4.0/3.0
    Kmax  = ini_struct["Kmax"]
    a     = ini_struct["a"]
    b     = ini_struct["b"]
    c     = ini_struct["c"]
    a     = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
    b     = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])
    c     = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])
    mdiv  = []
    nrdiv = []
    mra   = 2*int(ceil(scale*r43*Kmax/pi*a))
    mrb   = 2*int(ceil(scale*r43*Kmax/pi*b))
    mrc   = 2*int(ceil(scale*r43*Kmax/pi*c))
    nra   = 2*int(ceil(scale*Kmax/pi*a))
    nrb   = 2*int(ceil(scale*Kmax/pi*b))
    nrc   = 2*int(ceil(scale*Kmax/pi*c))
    #
    maxplw = 400000
    l1     = 2*nra+8
    l2     = 2*nrb+8
    l3     = 2*nrc+8
    pi     = acos(-1.0)
    # pi/6 comes from dividing the volume of a sphere by the volume of a
    # cube with a side equal to the sphere diameter
    numplw = (2*l1+1)*(2*l2+1)*(2*l3+1)*(pi/6.0)
    while numplw > maxplw:
        #
        # The maxplw parameter in GW set_control.F sets an upper limit 
        # on the number of plane waves. If the mdiv number of grid points
        # exceed this limit the calculation will crash with the error
        # messages:
        #
        #    ERROR from GETPLW : In sphere    0.0000 generated      0 vectors
        #
        # one should either increase maxplw in the code or reduce the 
        # the number of grid points. In this script the only intervention
        # we can enforce is to reduce the number of grid points, which is 
        # done here.
        #
        fac    = (float(maxplw)/float(numplw))**(1.0/3.0)
        nra    = int(nra*fac)
        nrb    = int(nrb*fac)
        nrc    = int(nrc*fac)
        l1     = 2*nra+8
        l2     = 2*nrb+8
        l3     = 2*nrc+8
        numplw = (2*l1+1)*(2*l2+1)*(2*l3+1)*(pi/6.0)
    #
    if SpGr == 166:
        #
        # The number of grid points along all lattice vectors has to
        # be the same to ensure that the realspace grid is 
        # commensurate with the space group symmetry.
        #
        mra = min(mra,mrb,mrc)
        mrb = mra
        mrc = mra
        nra = min(nra,nrb,nrc)
        nrb = nra
        nrc = nra
    mdiv.append(mra)
    mdiv.append(mrb)
    mdiv.append(mrc)
    nrdiv.append(nra)
    nrdiv.append(nrb)
    nrdiv.append(nrc)
    ini_struct["mdiv"]  = mdiv
    ini_struct["nrdiv"] = nrdiv
    return ini_struct

def establish_lmax(ini_struct):
    """
    For every element establish the maximum angular momentum to be used

    From the relation RK=l_max [1] where R is the muffin tin radius, K
    is the plane wave cut-off and l_max is the highest angular momentum
    of the orbitals of an atom we can find l_max for each element.
    The l_max values are stored in a dictionary.

    The values of l_max are returned in the ini_struct dictionary.

    #References#

    [1] D.J. Singh, L. Nordstrom, "Planewaves, pseudo-potentials and
        the LAPW method", Springer (1994), ISBN: 978-0-387-29684-5,
        pp. 62-63.
    """
    element_rad = ini_struct["element_rad"]
    elements    = ini_struct["symbol"]
    Kmax        = ini_struct["Kmax"]
    l_max       = {}
    for ii in range(0,len(elements)):
        eli        = elements[ii]
        nmi        = chemical_chrg[eli]
        rr         = element_rad[nmi]
        l_max[nmi] = int(rr*Kmax)
    ini_struct["l_max"] = l_max
    return ini_struct

def expand_atomic_basis(ini_struct):
    """
    For every element (atomic type) expand the basis up to the given l_max

    For every atomic type, using the valence basis set information and the
    maximum angular momentum, expand the stored basis set information to
    the full basis set. This includes adding missing angular momenta but
    also for the lower angular momenta add more radial functions up to 
    the appropriate N-quantum number.

    The corresponding dictionaries are returned in ini_struct. The names
    of the dictionaries are elmnt_exp_* corresponding to the global names
    element_* from which they were derived.
    """
    elements        = ini_struct["symbol"]
    l_max           = ini_struct["l_max"]
    elmnt_exp_lmb   = {}
    elmnt_exp_ntle  = {}
    elmnt_exp_augm  = {}
    elmnt_exp_atocc = {}
    elmnt_exp_ptnl  = {}
    elmnt_exp_idmd  = {}
    for ii in range(0,len(elements)):
        eli        = elements[ii]
        nmi        = chemical_chrg[eli]
        lmax       = l_max[nmi]
        Nmax       = lmax+1
        lstored    = len(element_ntle[nmi])-1
        augm       = []
        atocc      = []
        ptnl       = []
        idmd       = []
        ntle       = []
        iread      = 0
        for lcur in range(0,lstored+1): # lcur is the current angular momentum
            nfunc = element_ntle[nmi][lcur]
            Nhomo = 0
            for ifunc in range(0,nfunc):
                augm.append(element_augm[nmi][iread])
                atocc.append(element_atocc[nmi][iread])
                ptnl.append(element_ptnl[nmi][iread])
                idmd.append(element_idmd[nmi][iread])
                Ncur = int(element_ptnl[nmi][iread])
                if element_atocc[nmi][iread] > 0.0:
                    Nhomo = Ncur
                iread += 1
            Ncur += 1
            while Ncur <= min(Nmax,Nhomo+2):
                augm.append("LOC")
                atocc.append(0.0)
                ptnl.append(Ncur+0.8)
                idmd.append(1.0)
                nfunc += 1
                Ncur  += 1
            ntle.append(nfunc)
        for lcur in range(lstored+1,lmax+1):
            Ncur  = lcur+1
            nfunc = 0
            while Ncur <= lcur+1:
                augm.append("APW")
                atocc.append(0.0)
                ptnl.append(Ncur+0.8)
                idmd.append(0.0)
                nfunc += 1
                Ncur  += 1
            ntle.append(nfunc)
        elmnt_exp_lmb[nmi]   = lmax
        elmnt_exp_ntle[nmi]  = ntle
        elmnt_exp_augm[nmi]  = augm
        elmnt_exp_atocc[nmi] = atocc
        elmnt_exp_ptnl[nmi]  = ptnl
        elmnt_exp_idmd[nmi]  = idmd
    ini_struct["element_lmb"]   = elmnt_exp_lmb
    ini_struct["element_ntle"]  = elmnt_exp_ntle
    ini_struct["element_augm"]  = elmnt_exp_augm
    ini_struct["element_atocc"] = elmnt_exp_atocc
    ini_struct["element_ptnl"]  = elmnt_exp_ptnl
    ini_struct["element_idmd"]  = elmnt_exp_idmd
    return ini_struct

def write_inifile(ini_struct,inifile):
    """
    Take a dictionary with all the relevant ini settings and write an
    input file for the GW+DMFT code.
    """
    inifile.write("TEXT band structure calculation\n")
    inifile.write("CONTROL   iter_dft=%3i  iter_hf=%3i  iter_gw=%3i iter_qp=%3i\n" % 
                  (ini_struct["iter_dft"],ini_struct["iter_hf"],
                   ini_struct["iter_gw"],ini_struct["iter_qp"]))
    inifile.write("          admix=0.200  adspin=0.600 adm_gw=0.100 acc_it_gw=0.15\n")
    inifile.write("          iexch=005 scal_spin= 1.0000\n")
    inifile.write("          nproc_tau=   1 nproc_k=   1\n")
    inifile.write("          irel=1 clight=274.074e+00 rel_interst=F irel_core=1\n")
    if ini_struct["restart"]:
        restart="T"
    else:
        restart="F"
    inifile.write("          temperature=%10.2f  restart=%s\n" % (ini_struct["temperature"],restart))
    inifile.write("FILES\n")
    inifile.write("  allfile=mdl\n")
    if True:
        symmops = ini_struct["symgen"]
        str_inp=""
        ngroup=1
        for symmop in symmops:
            str_op=""
            if symmop.kind == "I":
                ngroup*=2
                if sqrt(numpy.dot(symmop.translation_vector,
                                  symmop.translation_vector))<1.0e-10:
                    str_op   = "I_"
                else:
                    str_op   = "I:T["
                    str_op  += str(symmop.translation_vector[0])+","
                    str_op  += str(symmop.translation_vector[1])+","
                    str_op  += str(symmop.translation_vector[2])+"]_"
            elif symmop.kind == "M":
                ngroup*=2
                if sqrt(numpy.dot(symmop.translation_vector,
                                  symmop.translation_vector))<1.0e-10:
                    str_op   = "M["
                    str_op  += str(symmop.axis[0])+","
                    str_op  += str(symmop.axis[1])+","
                    str_op  += str(symmop.axis[2])+"]_"
                else:
                    str_op   = "M["
                    str_op  += str(symmop.axis[0])+","
                    str_op  += str(symmop.axis[1])+","
                    str_op  += str(symmop.axis[2])+"]:T["
                    str_op  += str(symmop.translation_vector[0])+","
                    str_op  += str(symmop.translation_vector[1])+","
                    str_op  += str(symmop.translation_vector[2])+"]_"
            elif symmop.kind == "R":
                ngroup*=symmop.n
                if sqrt(numpy.dot(symmop.translation_vector,
                                  symmop.translation_vector))<1.0e-10:
                    str_op   = "R"+str(symmop.n)+"["
                    str_op  += str(symmop.axis[0])+","
                    str_op  += str(symmop.axis[1])+","
                    str_op  += str(symmop.axis[2])+"]_"
                else:
                    str_op   = "R"+str(symmop.n)+"["
                    str_op  += str(symmop.axis[0])+","
                    str_op  += str(symmop.axis[1])+","
                    str_op  += str(symmop.axis[2])+"]:T["
                    if abs(int(symmop.translation_vector[0])-
                               symmop.translation_vector[0]) < 1.0e-10:
                        str_op  += str(symmop.translation_vector[0])+","
                    else:
                        str_op  += '{:.16f}'.format(symmop.translation_vector[0])+","
                    if abs(int(symmop.translation_vector[1])-
                               symmop.translation_vector[1]) < 1.0e-10:
                        str_op  += str(symmop.translation_vector[1])+","
                    else:
                        str_op  += '{:.16f}'.format(symmop.translation_vector[1])+","
                    if abs(int(symmop.translation_vector[2])-
                               symmop.translation_vector[2]) < 1.0e-10:
                        str_op  += str(symmop.translation_vector[2])+"]_"
                    else:
                        str_op  += '{:.16f}'.format(symmop.translation_vector[2])+"]_"
            if ngroup <= 48:
                str_inp += str_op
        inifile.write("SYM symgen=%s\n" % str_inp)
    else:
        inifile.write("SYM symgen=input\n")
        nsymmops = len(ini_struct["symgen"])
        # The GW code can handle at most 48 symmetry operations
        # It seems to be common to skip symmetry generators to limit
        # the number of symmetry operations. Here we just cut the 
        # symmetry operations off at number 48.
        if nsymmops > 48:
           nsymmops = 48
        inifile.write("    number of symmetry operations=%3d\n" % (nsymmops))
        for ii in range(0,nsymmops):
           inifile.write("    symmetry operation %3d\n" % (ii+1))
           symmop = ini_struct["symgen"][ii]
           rotation = symmop.rotation_matrix
           translation = symmop.translation_vector
           for jj in range(0,3):
              line = (rotation[jj][0],rotation[jj][1],rotation[jj][2],translation[jj])
              inifile.write("    (%14.10f,%14.10f,%14.10f)    (%14.10f)\n" % line)
    nsort = len(ini_struct["symeqv"])
    inifile.write("STRUCTURE  par=%11.7f  natom=%3d nsort=%3d istruct=%3d\n" %
                  (ini_struct["par"],ini_struct["natom"],nsort,ini_struct["istruc"]))
    inifile.write("      is=")
    islist = ini_struct["islist"]
    for isnum in islist:
        inifile.write("%3d" % isnum)
    inifile.write("\n")
    inifile.write("      b/a=%9.6f  c/a=%9.6f\n" % (ini_struct["b_a"],ini_struct["c_a"]))
    inifile.write("      a=%21.16f%21.16f%21.16f\n" % ini_struct["a"])
    inifile.write("      b=%21.16f%21.16f%21.16f\n" % ini_struct["b"])
    inifile.write("      c=%21.16f%21.16f%21.16f\n" % ini_struct["c"])
    natom = ini_struct["natom"]
    for ii in range(0,natom):
        inifile.write("    tau=%21.16f%21.16f%21.16f\n" %
                      (ini_struct["xcoord"][ii],ini_struct["ycoord"][ii],ini_struct["zcoord"][ii]))
    mdiv  = ini_struct["mdiv"]
    nrdiv = ini_struct["nrdiv"]
    inifile.write("REAL SPACE MESHES mdiv=%4d %4d %4d\n" % (mdiv[0], mdiv[1], mdiv[2]))
    inifile.write("                 nrdiv=%4d %4d %4d\n" % (nrdiv[0],nrdiv[1],nrdiv[2]))
    inifile.write("BASIS  cut_lapw_ratio=0.610 cut_pb_ratio=0.980\n")
    inifile.write("       eps_pb=1.e-03\n")
    inifile.write("ZONES nbndf=   0\n")
    if ini_struct["band"]:
        band = ' T'
    else:
        band = ' F'
    if ini_struct["dos"]:
        dos = ' T'
    else:
        dos = ' F'
    inifile.write("DOS   emindos=-15.000  emaxdos= 15.000   ndos=  800\n")
    inifile.write("      n_cont_frac=  30 e_small=2.e-02\n")
    inifile.write("      dos=%s           bandstructure=%s\n" % (dos,band))
    inifile.write("K_POINT  ndiv=   4   4   4  metal=T n_k_div= 27 k_line=010\n")
    inifile.write("MULTI_SCF vv0=  1.00\n")
    inifile.write("MAGNET  b_extval=   0.000000 iter_h_ext=0000100\n")
    inifile.write("        b_ext=  0.000  0.000  1.000\n")
    inifile.write("TAU MESH n_tau=   46 n_tau_int= 1200\n")
    inifile.write("OMEGA MESH n_omega_exa=   29 n_omega_asy=   18 omega_max=  200.00 \n")
    inifile.write("           interp_omega_d= 2\n")
    inifile.write("NU MESH n_nu_exa=   29 n_nu_asy=   18 nu_max=  200.00\n")
    inifile.write("        interp_nu_d= 2\n")
    inifile.write("ATOMIC DATA --------------------------------------------------------\n")
    element_rad = ini_struct["element_rad"]
    isold = 0
    for ii in range(0,natom):
        isnew = ini_struct["islist"][ii]
        if isold != isnew:
            symbol = ini_struct["symbol"][ii]
            if len(symbol) == 1:
                symb2 = symbol+"  "
            elif len(symbol) == 2:
                symb2 = symbol+" "
            elif len(symbol) == 3:
                symb2 = symbol
            else:
                error = "Strange chemical symbol:"+symbol
                print(error)
            number = chemical_chrg[symbol]
            smt    = element_rad[number]*rmt_reduce
            inifile.write("  txtel=%s   z=%5.1f magn_shift= 0.050\n" %
                          (symb2,float(number)))
            inifile.write("  smt=%8.5f h= 0.0120 nrad= 1216 z_dop=0.000\n" % 
                          smt)
            lmb = ini_struct["element_lmb"][number]
            # lmpb is the maximum l-quantum number for the product basis
            # we set this to the same value as lmb for now...
            inifile.write("  lmb=%2d  lmpb=%2d\n" % (lmb,lmb))
            ntle = ini_struct["element_ntle"][number]
            inifile.write("  lim_pb_mt=")
            for ii in range(0,len(ntle)):
                inifile.write("%3d" % 30)
            inifile.write("\n")
            inifile.write("  ntle=")
            for ii in range(0,len(ntle)):
                inifile.write("%3d" % ntle[ii])
            inifile.write("\n")
            inifile.write("  l  augm  atocc   ptnl  corr idmd\n")
            kk = 0
            for ii in range(0,len(ntle)):
                ntlen = ntle[ii]
                for jj in range(0,ntlen):
                    inifile.write("%3d   %s%7.3f%7.3f     %s    %1d\n" %
                                  (ii,ini_struct["element_augm"][number][kk],
                                   ini_struct["element_atocc"][number][kk],
                                   ini_struct["element_ptnl"][number][kk],"N",
                                   ini_struct["element_idmd"][number][kk]))
                    kk+=1
        isold = isnew

def write_kpathfile(ini_struct,kpathfile):
    """
    Take a dictionary with all the relevant information for the
    structure, extract the Kpath, and write the data to the
    kpathfile.
    """
    kpath = ini_struct["kpath"]
    if not kpath:
        raise ValueError("None is an invalid value")
    kpoints = kpath["kpoints"]
    temp = itertools.chain.from_iterable(kpath["path"])
    path = []
    for point in temp:
        path.append(point)
    length = len(path)
    for ii in range(0,length-1,2):
        point0 = str(path[ii])
        point1 = str(path[ii+1])
        kpoint0 = kpoints[point0]
        kpoint1 = kpoints[point1]
        point0 = str.replace(point0,"\\Gamma","G")
        point1 = str.replace(point1,"\\Gamma","G")
        # Cannot replace Sigma with S because S is already taken
        point0 = str.replace(point0,"\\Sigma","s")
        point1 = str.replace(point1,"\\Sigma","s")
        kpathfile.write("%s %12.8f %12.8f %12.8f       %s %12.8f %12.8f %12.8f\n" %
                        (point0,kpoint0[0],kpoint0[1],kpoint0[2],
                         point1,kpoint1[0],kpoint1[1],kpoint1[2]))
    oddlength = ((length % 2) == 1)
    if oddlength:
        point0 = str(path[length-1])
        kpoint0 = kpoints[point0]
        point0 = str.replace(point0,"\\Gamma","G")
        # Cannot replace Sigma with S because S is already taken
        point0 = str.replace(point0,"\\Sigma","s")
        kpathfile.write("%s %12.8f %12.8f %12.8f\n" %
                        (point0,kpoint0[0],kpoint0[1],kpoint0[2]))

def write_plot1dfile(ini_struct,plot1dfile):
    """
    Take a dictionary with all the relevant information for the
    structure, extract the Kpath, and write the data to the
    plot1dfile for the Elk code.
    """
    kpath = ini_struct["kpath"]
    if not kpath:
        raise ValueError("None is an invalid value")
    kpoints = kpath["kpoints"]
    temp = itertools.chain.from_iterable(kpath["path"])
    path = []
    for point in temp:
        path.append(point)
    length = len(path)
    plot1dfile.write("plot1d\n")
    plot1dfile.write("  %d\n" % length)
    plot1dfile.write("  200\n")
    for ii in range(0,length):
        point = str(path[ii])
        kpoint = kpoints[point]
        plot1dfile.write("%12.8f %12.8f %12.8f\n" %
                        (kpoint[0],kpoint[1],kpoint[2]))

def write_kpointfile(ini_struct,kpointfile):
    """
    Take a dictionary with all the relevant information for the
    structure, extract the Kpoints, and write the data to the
    kpointfile.
    """
    if not ini_struct["kpoints"]:
        raise ValueError("None is an invalid value")
    (kpoints,kpath) = ini_struct["kpoints"]
    length = len(kpath)
    for ii in range(0,length):
        point  = str(kpath[ii])
        kpoint = kpoints[ii]
        point  = str.replace(point,"\\Gamma","G")
        # Cannot replace Sigma with S because S is already taken
        point  = str.replace(point,"\\Sigma","s")
        kpointfile.write("  %12.8f %12.8f %12.8f   %s\n" %
                         (kpoint[0],kpoint[1],kpoint[2],point))

def write_klistfile(ini_struct,kpointfile):
    """
    Take a dictionary with all the relevant information for the
    structure, extract the Kpoints, and write the data to the
    Wien2k klistfile.

    The Wien2K klist format is special in that the k-points are
    given in terms of integers instead of fractional coordinates.
    So we need to convert all fractional coordinates of each
    k-point to integers first.
    """
    import fractions
    if not ini_struct["kpoints"]:
        raise ValueError("None is an invalid value")
    (kpoints,kpath) = ini_struct["kpoints"]
    length = len(kpath)
    for ii in range(0,length):
        point  = str(kpath[ii])
        kpoint = kpoints[ii]
        fracx  = fractions.Fraction.from_float(kpoint[0]).limit_denominator(10000)
        fracy  = fractions.Fraction.from_float(kpoint[1]).limit_denominator(10000)
        fracz  = fractions.Fraction.from_float(kpoint[2]).limit_denominator(10000)
        idv    = fracx.denominator
        if fracy.denominator != fracx.denominator:
            idv *= fracy.denominator
        if (fracz.denominator != fracy.denominator and 
            fracz.denominator != fracx.denominator):
            idv *= fracz.denominator
        point = str.replace(point,"\\Gamma","G")
        # Cannot replace Sigma with S because S is already taken
        point = str.replace(point,"\\Sigma","s")
        isx   = int(kpoint[0]*idv+0.5)
        isy   = int(kpoint[1]*idv+0.5)
        isz   = int(kpoint[2]*idv+0.5)
        kpointfile.write("%-10s%5d%5d%5d%5d%5.2f%5.2f%5.2f%3s\n" %
                         (point,isx,isy,isz,idv,2.0,0.0,0.0,"   "))
    kpointfile.write("END\n")

def any2utf8(infile):
    """
    Convert the encoding of a text file to UTF8 encoding

    CIF files are supposed to contain only ASCII characters but some sources
    generate CIF files with other characters. This may cause problems if 
    these files are encoded in anything else but UTF8. Python and therefore
    Pymatgen expect CIF files to be UTF8 encoded. If they are not the code
    will throw an exception. 

    To avoid that problem this function takes the name of a CIF file. It
    detects the encoding of that file, reads it in using that encoding,
    writes out the contents to a temporary file using UTF8 encoding,
    and returns the name of the temporary file.
    """
    import codecs
    import os
    from chardet.universaldetector import UniversalDetector
    #
    # Generate a unique filename
    #
    pid = os.getpid()
    outfile = "/tmp/tmp"+str(pid)+".cif"
    #
    # Detect the encoding of the input file
    #
    detector = UniversalDetector()
    detector.reset()
    file = open(infile,'rb')
    for line in file:
        detector.feed(line)
        if detector.done:
            break
    detector.close()
    file.close()
    #
    # Read the input file and decode the contents
    #
    codec = detector.result['encoding']
    f = codecs.open(infile,mode='r',encoding=codec)
    contents_in = f.readlines()
    f.close()
    #
    # Map all characters to ASCII characters (non ASCII characters are skipped)
    #
    contents_out = []
    for line in contents_in:
        line_out = ""
        for char in line:
            if ord(char) <= 127:
                line_out = line_out + char
        contents_out.append(line_out)
    #
    # Write the contents to the temporary file UTF8 encoded
    #
    f = codecs.open(outfile,mode='w',encoding='utf8')
    f.writelines(contents_out)
    f.close()
    #
    # Return the temporary filename
    #
    return outfile

def write_comsuite(ini):
    """
    Take the input data and generate the input files for the Comsuite of
    programs.
    """
    #
    # Create an .ini file for data of <key> structure
    #
    filename = "ini"
    inifile = open(filename,'w')
    #
    # Write MatDeLab.ini
    #
    write_inifile(ini,inifile)
    #
    # Close .ini file
    #
    inifile.close()    
    #
    # Create an .kpath file for data of <key> structure
    #
    filename = "kpath"
    kpathfile = open(filename,'w')
    #
    # Write MatDeLab.kpath
    #
    try:
        write_kpathfile(ini,kpathfile)
        kpathfile.close()
    except ValueError:
        kpathfile.close()
        os.remove(filename)
    #
    # Create an .kpoint file for data of <key> structure
    #
    filename = "kpoints"
    kpointfile = open(filename,'w')
    #
    # Write MatDeLab.kpath
    #
    try:
        write_kpointfile(ini,kpointfile)
        kpointfile.close()
    except ValueError:
        kpointfile.close()
        os.remove(filename)

def unique_species(ini_struct):
    """
    Return the list of different chemical elements there are in the
    current structure.
    """
    natom = ini_struct["natom"]
    elmlist = []
    for ii in range(0,natom):
        symbol = ini_struct["symbol"][ii]
        if not symbol in elmlist:
            elmlist.append(symbol)
    return elmlist

def write_elkfile(ini_struct,elkfile):
    """
    Take a dictionary with all the relevant input settings and write an
    input file for the Elk program.
    """
    elkfile.write("tasks\n")
    elkfile.write("  0\n")
    if ini_struct["dos"]:
        elkfile.write("  10\n")
    if ini_struct["band"]:
        elkfile.write("  20\n")
        elkfile.write("  21\n")
    elkfile.write("\n")
    elkfile.write("isgkmax\n")
    elkfile.write("  -3\n")
    elkfile.write("\n")
    elkfile.write("spinpol\n")
    elkfile.write("  .true.\n")
    elkfile.write("\n")
    if ini_struct["dos"] or ini_struct["band"]:
        # vhighq seems rather expensive to run, maybe highq is good enough
        elkfile.write("highq\n")
        elkfile.write("  .true.\n")
        elkfile.write("\n")
    elkfile.write("tempk\n")
    elkfile.write("  %s\n" % ini_struct["temperature"])
    elkfile.write("\n")
    elkfile.write("scale\n")
    elkfile.write("  %11.7f\n" % ini_struct["par"])
    elkfile.write("\n")
    elkfile.write("avec\n")
    elkfile.write("  %21.16f %21.16f %21.16f\n" % ini_struct["a"])
    elkfile.write("  %21.16f %21.16f %21.16f\n" % ini_struct["b"])
    elkfile.write("  %21.16f %21.16f %21.16f\n" % ini_struct["c"])
    elk_species_path = os.environ.get('ELK_SPECIES_PATH')
    if not elk_species_path:
        error = "Environment variable ELK_SPECIES_PATH not set"
        elk_species_path = "."
        print(error)
    elkfile.write("\n")
    elkfile.write("atoms\n")
    natom = ini_struct["natom"]
    elmlist = unique_species(ini_struct)
    nelm = len(elmlist)
    elkfile.write("  %d\n" % nelm)
    for element in elmlist:
        elmname = element.strip()
        elmname = elmname.capitalize()
        elkfile.write("  '%s.in'\n" % (elk_species_path + "/" + elmname) )
        elkfile.write("  %d\n" % ini_struct["symbol"].count(element) )
        for ii in range(0,natom):
            symbol = ini_struct["symbol"][ii]
            if element == symbol:
                elkfile.write("  %21.16f %21.16f %21.16f 0.0 0.0 0.0\n" % (ini_struct["xcoord"][ii],ini_struct["ycoord"][ii],ini_struct["zcoord"][ii]))
    elkfile.write("\n")
    elkfile.write("nempty\n")
    if ini_struct["dos"] or ini_struct["band"]:
        elkfile.write("  30\n")
    else:
        elkfile.write("  5\n")
    elkfile.write("\n")
    elkfile.write("ngridk\n")
    elkfile.write("  %d %d %d\n" % (ini_struct["mdiv"][0],ini_struct["mdiv"][1],ini_struct["mdiv"][2]))
    elkfile.write("\n")
    try:
        write_plot1dfile(ini_struct,elkfile)
    except ValueError:
        pass

def write_elk(ini):
    """
    Take the input data and generate the input files for the Elk
    program.
    """
    #
    # Create an .ini file for data of <key> structure
    #
    filename = "elk.in"
    elkfile = open(filename,'w')
    #
    # Write elk.in
    #
    write_elkfile(ini,elkfile)
    #
    # Close elk.in file
    #
    elkfile.close()    
    #
    # Create an .kpoint file for data of <key> structure
    #
    filename = "kpoints"
    kpointfile = open(filename,'w')
    #
    # Write MatDeLab.kpath
    #
    try:
        write_kpointfile(ini,kpointfile)
        kpointfile.close()
    except ValueError:
        kpointfile.close()
        os.remove(filename)

def write_wien2k(ini):
    """
    Take the input data and generate the input files for Wien2K 
    package.
    """
    import sys
    import subprocess
    case = ini["cif"].split(".")[0]
    case_st2 = case+".struct"
    case_st1 = case_st2+"_st"
    case_st3 = "/tmp/"+case_st2
    make_struct = []
    make_struct.append("cif2struct")
    make_struct.append("/tmp/"+os.path.basename(ini["cif"]))
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_st3)
    make_struct.append(case_st2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("x")
    make_struct.append("symmetry")
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_st1)
    make_struct.append(case_st2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    # Run setrmt_lapw to choose the muffin tin radii.
    # The -r flag specifies the percentage reduction of the radii from
    # the just-touching radii. This is a REQUIRED flag because 
    # setrmt_lapw rounds the radii after calculating them to 2 decimal
    # places. The test whether the spheres are overlapping uses at
    # least 5 decimal places. Hence, if the spheres are not reduced
    # the rounding of the radii may cause the non-overlapping requirement
    # to be violated, and the calculation will abort!
    #
    # - setrmt_lapw case -r 3
    #
    case_st2 = case+".struct"
    case_st1 = case_st2+"_setrmt"
    make_struct = []
    make_struct.append("setrmt_lapw")
    make_struct.append(case)
    make_struct.append("-r")
    make_struct.append("3")
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_st1)
    make_struct.append(case_st2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()

    case_ins = case+".inst"
    make_struct = []
    make_struct.append("rm")
    make_struct.append("-f")
    make_struct.append(case_ins)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("instgen_lapw")
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("x")
    make_struct.append("lstart")
    process = subprocess.Popen(make_struct,stdin=subprocess.PIPE)
    outs,errs = process.communicate(bytes("5\n-6.0\n","utf-8"))
    case_in2 = case+".in0"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    case_in2 = case+".in1"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    case_in2 = case+".vsp"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    case_in3 = case+".in2"
    case_in1 = case_in3+"_ls"
    case_in2 = case_in3+"_sy"
    #make_struct = []
    #make_struct.append("cat")
    #make_struct.append(case_in1)
    #make_struct.append(case_in2)
    #make_struct.append(">")
    #make_struct.append(case_in3)
    line = "cat "+str(case_in1)+" "+str(case_in2)+" > "+str(case_in3)
    result = subprocess.run(line,shell=True)
    #if sys.version_info.major==3:
    #    result = result.decode()
    #
    # If Wien2K thinks there is no inversion symmetry the code needs 
    # .in1c and .in2c files instead of .in1 and .in2 files. 
    # To generate the former files just copy the latter.
    #
    # - cp case.in1 case.in1c
    #
    case_in1 = case+".in1"
    case_in2 = case_in1+"c"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    # - cp case.in2 case.in2c
    #
    case_in1 = case+".in2"
    case_in2 = case_in1+"c"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    make_struct = []
    make_struct.append("x")
    make_struct.append("kgen")
    process = subprocess.Popen(make_struct,stdin=subprocess.PIPE)
    line = "0\n"+str(ini["mdiv"][0])+" "+str(ini["mdiv"][1])+" "+str(ini["mdiv"][2])+"\n0"
    print("line = %s\n" % line)
    outs,errs = process.communicate(bytes(line,"utf-8"))
    make_struct = []
    make_struct.append("x")
    make_struct.append("dstart")
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    case_in2 = case+".inm"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    case_in2 = case+".inc"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    # Create the .klist_band file
    #
    case_in1 = case+".klist_band"
    klist_file = open(case_in1,"w")
    try:
        write_klistfile(ini,klist_file)
        klist_file.close()
    except ValueError:
        klist_file.close()
        os.remove(case_in1)

    
def execute_with_arguments(args):
    """
    Run the input generator for the method and ciffile as specified by the 
    command line arguments.
    """
    #
    # Convert data from CifFile to MatDeLab.ini
    #
    cell    = str(args.cell)
    tmpfile = any2utf8(args.ciffile)
    struct  = cif2matdelab.structure(tmpfile,cell)
    struct.struct.to(fmt="cif",filename="/tmp/"+os.path.basename(str(args.ciffile)))
    struct.prints()
    Kmax   = float(args.Kmax)
    method = str(args.method)
    code   = str(args.code)
    ini = {}
    ini["cif"] = args.ciffile
    ini["par"] = 1.0
    ini["b_a"] = 1.0
    ini["c_a"] = 1.0
    ini["restart"] = args.restart
    ini["band"]    = (args.band == 'Y')
    ini["dos"]     = (args.dos  == 'Y')
    ini = cif2matdelab.retr_cell_volume(struct,ini)
    ini = cif2matdelab.retr_lattice_vecs(struct,ini)
    ini = cif2matdelab.retr_sites(struct,ini)
    ini = cif2matdelab.retr_distance_matrix(struct,ini)
    ini = cif2matdelab.retr_lattice_type(struct,ini)
    if True:
        ini = cif2matdelab.retr_symmetry_generators(struct,ini)
    else:
        ini = cif2matdelab.retr_symmetry_operations(struct,ini)
    ini = cif2matdelab.retr_kpath(struct,ini)
    ini = cif2matdelab.retr_kpoints(struct,ini)
    ini = cif2matdelab.retr_spacegroup_number(struct,ini)
    #
    # Convert some data items into useful information
    #
    ini = translate_elements(ini)
    ini = establish_mt_radii(ini)
    ini = establish_atoms_volume(ini)
    ini = establish_Kmin(ini)
    ini = establish_Kmax(ini,Kmax=Kmax)
    ini = establish_r_grid(ini)
    ini = establish_lmax(ini)
    ini = establish_method(ini,method)
    if args.dftit!=None:
        ini["iter_dft"] = args.dftit
    if args.hfit!=None:
        ini["iter_hf"]  = args.hfit
    if args.gwit!=None:
        ini["iter_gw"]  = args.gwit
    if args.qpit!=None:
        ini["iter_qp"]  = args.qpit
    if args.temperature:
        ini["temperature"]  = args.temperature
    ini = expand_atomic_basis(ini)
    #
    # Create input files
    #
    if code == "comsuite":
        write_comsuite(ini)
    elif code == "elk":
        write_elk(ini)
    elif code == "wien2k":
        write_wien2k(ini)
    else:
        error = "Unknown code suite: "+code+" No input files generated"
        print(error)
        


def main():
    execute_with_arguments(parse_arguments())

if __name__ == "__main__":
    main()
