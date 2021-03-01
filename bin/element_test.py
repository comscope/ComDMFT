import string as string
import matplotlib as mat
import re as re
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import matplotlib.font_manager as fm
from collections import OrderedDict
import json, os, shutil

    

chemical_name = {} # element name by atomic number
chemical_symb = {} # chemical symbol by atomic number
chemical_chrg = {} # nuclear charge by chemical symbol
rmt   = {} # muffin tin radius by atomic number
lmb   = {} # maximum l-quantum number for the muffin tin basis by atomic number
                   # incidentally this is also the minimum l-quantum number required to 
                   # describe all occupied orbitals
ntle  = {} # by atomic number a list of the number of basis functions for each l-quantum number
augm  = {} # by atomic number a list of nature of the basis functions LOC/APW
atocc = {} # by atomic number a list of the orbital occupation
ptnl  = {} # by atomic number a list of logarithmic derivatives 
idmd  = {} # by atomic number a list of core/valence(=0) or unoccupied(=1) status of basis function

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

# Data for rmt taken from the species files of Elk 3.1.12
# for elements 1-104. The data for elements 105-118 are set based
# on similarity to the lower atomic number elements.
# The rmt are given in Bohr.
rmt[1]   = 1.4
rmt[2]   = 1.4
rmt[3]   = 1.8
rmt[4]   = 1.8
rmt[5]   = 1.8
rmt[6]   = 1.8
rmt[7]   = 1.8
rmt[8]   = 1.8
rmt[9]   = 2.0
rmt[10]  = 1.6
rmt[11]  = 2.2
rmt[12]  = 2.2
rmt[13]  = 2.2
rmt[14]  = 2.2
rmt[15]  = 2.2
rmt[16]  = 2.2
rmt[17]  = 2.2
rmt[18]  = 2.0
rmt[19]  = 2.4
rmt[20]  = 2.4
rmt[21]  = 2.4
rmt[22]  = 2.4
rmt[23]  = 2.4
rmt[24]  = 2.4
rmt[25]  = 2.4
rmt[26]  = 2.4
rmt[27]  = 2.4
rmt[28]  = 2.4
rmt[29]  = 2.4
rmt[30]  = 2.4
rmt[31]  = 2.4
rmt[32]  = 2.4
rmt[33]  = 2.4
rmt[34]  = 2.4
rmt[35]  = 2.4
rmt[36]  = 2.2
rmt[37]  = 2.6
rmt[38]  = 2.6
rmt[39]  = 2.6
rmt[40]  = 2.6
rmt[41]  = 2.6
rmt[42]  = 2.6
rmt[43]  = 2.6
rmt[44]  = 2.6
rmt[45]  = 2.6
rmt[46]  = 2.6
rmt[47]  = 2.6
rmt[48]  = 2.6
rmt[49]  = 2.6
rmt[50]  = 2.6
rmt[51]  = 2.6
rmt[52]  = 2.6
rmt[53]  = 2.6
rmt[54]  = 2.4
rmt[55]  = 2.8
rmt[56]  = 2.8
rmt[57]  = 2.8
rmt[58]  = 2.8
rmt[59]  = 2.8
rmt[60]  = 2.8
rmt[61]  = 2.8
rmt[62]  = 2.8
rmt[63]  = 2.8
rmt[64]  = 2.8
rmt[65]  = 2.8
rmt[66]  = 2.8
rmt[67]  = 2.8
rmt[68]  = 2.8
rmt[69]  = 2.8
rmt[70]  = 2.8
rmt[71]  = 2.8
rmt[72]  = 2.8
rmt[73]  = 2.8
rmt[74]  = 2.8
rmt[75]  = 2.8
rmt[76]  = 2.8
rmt[77]  = 2.8
rmt[78]  = 2.8
rmt[79]  = 2.8
rmt[80]  = 2.8
rmt[81]  = 2.8
rmt[82]  = 2.8
rmt[83]  = 2.8
rmt[84]  = 2.8
rmt[85]  = 2.8
rmt[86]  = 2.6
rmt[87]  = 3.0
rmt[88]  = 3.0
rmt[89]  = 3.0
rmt[90]  = 3.0
rmt[91]  = 3.0
rmt[92]  = 3.0
rmt[93]  = 3.0
rmt[94]  = 3.0
rmt[95]  = 3.0
rmt[96]  = 3.0
rmt[97]  = 3.0
rmt[98]  = 3.0
rmt[99]  = 3.0
rmt[100] = 3.0
rmt[101] = 3.0
rmt[102] = 3.0
rmt[103] = 3.0
rmt[104] = 3.0

rmt[105] = 3.0
rmt[106] = 3.0
rmt[107] = 3.0
rmt[108] = 3.0
rmt[109] = 3.0
rmt[110] = 3.0
rmt[111] = 3.0
rmt[112] = 3.0
rmt[113] = 3.0
rmt[114] = 3.0
rmt[115] = 3.0
rmt[116] = 3.0
rmt[117] = 3.0
rmt[118] = 2.8

# Below we initialize lmb, ntle, augm, atocc, ptnl, and
# idmd. 
# Because these fields for a given element are related it makes sense to group the initialization
# of all these fields by element. 
# Data from WIEN2k *.outputst files

lmb[1]   = 1
ntle[1]  = [2,1]
augm[1]  = ["APW","LOC",    "APW"]
atocc[1] = [ 1.0,  0.0,      0.0]
ptnl[1]  = [ 1.8,  2.8,      2.8]
idmd[1]  = [  0,    1,        0]


lmb[2]   = 1
ntle[2]  = [2,1]
augm[2]  = ["APW","LOC",    "APW"]
atocc[2] = [ 2.0,  0.0,      0.0]
ptnl[2]  = [ 1.8,  2.8,      2.8]
idmd[2]  = [  0,    1,        0]


lmb[3]   = 2
ntle[3]  = [3,2,1]
augm[3]  = ["LOC","APW","LOC",    "APW","LOC",     "APW"]
atocc[3] = [ 2.0,  1.0,  0.0,      0.0,  0.0,       0.0]
ptnl[3]  = [ 1.8,  2.8,  3.8,      2.8,  3.8,       3.8]
idmd[3]  = [  0,    0,    1,        0,     1,        0]  


lmb[4]   = 2
ntle[4]  = [2,2,1]
augm[4]  = ["APW","LOC",    "APW","LOC",     "APW"]
atocc[4] = [ 2.0,  0.0,      0.0,  0.0,       0.0]
ptnl[4]  = [ 2.8,  3.8,      2.8,  3.8,       3.8]
idmd[4]  = [  0,    1,        0,     1,        0 ] 

lmb[5]   = 2                                       
ntle[5]  = [2,2,1]                                 
augm[5]  = ["APW","LOC",    "APW","LOC",     "APW"]
atocc[5] = [ 2.0,  0.0,      1.0,  0.0,       0.0] 
ptnl[5]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
idmd[5]  = [  0,    1,        0,     1,        0 ] 


lmb[6]   = 2                                       
ntle[6]  = [2,2,1]                                 
augm[6]  = ["APW","LOC",    "APW","LOC",     "APW"]
atocc[6] = [ 2.0,  0.0,      2.0,  0.0,       0.0] 
ptnl[6]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
idmd[6]  = [  0,    1,        0,     1,        0 ] 


lmb[7]   = 2                                       
ntle[7]  = [2,2,1]                                 
augm[7]  = ["APW","LOC",    "APW","LOC",     "APW"]
atocc[7] = [ 2.0,  0.0,      3.0,  0.0,       0.0] 
ptnl[7]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
idmd[7]  = [  0,    1,        0,     1,        0 ] 


lmb[8]   = 2                                       
ntle[8]  = [2,2,1]                                 
augm[8]  = ["APW","LOC",    "APW","LOC",     "APW"]
atocc[8] = [ 2.0,  0.0,      4.0,  0.0,       0.0] 
ptnl[8]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
idmd[8]  = [  0,    1,        0,     1,        0 ] 


lmb[9]   = 2                                       
ntle[9]  = [2,2,1]                                 
augm[9]  = ["APW","LOC",    "APW","LOC",     "APW"]
atocc[9] = [ 2.0,  0.0,      5.0,  0.0,       0.0] 
ptnl[9]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
idmd[9]  = [  0,    1,        0,     1,        0 ] 


lmb[10]   = 2                                       
ntle[10]  = [2,2,1]                                 
augm[10]  = ["APW","LOC",    "APW","LOC",     "APW"]
atocc[10] = [ 2.0,  0.0,      6.0,  0.0,       0.0] 
ptnl[10]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
idmd[10]  = [  0,    1,        0,     1,        0 ] 

lmb[11]   = 2
ntle[11]  = [3,3,1]
augm[11]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[11] = [ 2.0,  1.0,  0.0,       6.0,  0.0,  0.0,       0.0]
ptnl[11]  = [ 2.8,  3.8,  4.8,       2.8,  3.8,  4.8,       3.8]
idmd[11]  = [  0,    0,    1,         0,    0,    1,         0] 


lmb[12]   = 2                                                    
ntle[12]  = [2,3,1]                                              
augm[12]  = ["APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[12] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       0.0] 
ptnl[12]  = [ 3.8,  4.8,       2.8,  3.8,  4.8,       3.8] 
idmd[12]  = [  0,    1,         0,    0,    1,         0 ] 


lmb[13]   = 2                                              
ntle[13]  = [2,2,1]                                        
augm[13]  = ["APW","LOC",     "APW","LOC",     "APW"]
atocc[13] = [ 2.0,  0.0,       1.0,  0.0,       0.0] 
ptnl[13]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
idmd[13]  = [  0,    1,         0,    1,         0 ] 


lmb[14]   = 2                                              
ntle[14]  = [2,2,1]                                        
augm[14]  = ["APW","LOC",     "APW","LOC",     "APW"]
atocc[14] = [ 2.0,  0.0,       2.0,  0.0,       0.0] 
ptnl[14]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
idmd[14]  = [  0,    1,         0,    1,         0] 


lmb[15]   = 2                                              
ntle[15]  = [2,2,1]                                        
augm[15]  = ["APW","LOC",     "APW","LOC",     "APW"]
atocc[15] = [ 2.0,  0.0,       3.0,  0.0,       0.0] 
ptnl[15]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
idmd[15]  = [  0,    1,         0,    1,         0 ] 

lmb[16]   = 2                                              
ntle[16]  = [2,2,1]                                        
augm[16]  = ["APW","LOC",     "APW","LOC",     "APW"]
atocc[16] = [ 2.0,  0.0,       4.0,  0.0,       0.0] 
ptnl[16]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
idmd[16]  = [  0,    1,         0,    1,         0] 


lmb[17]   = 2                                              
ntle[17]  = [2,2,1]                                        
augm[17]  = ["APW","LOC",     "APW","LOC",     "APW"]
atocc[17] = [ 2.0,  0.0,       5.0,  0.0,       0.0] 
ptnl[17]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
idmd[17]  = [  0,    1,         0,    1,         0] 


lmb[18]   = 2                                              
ntle[18]  = [2,2,1]                                        
augm[18]  = ["APW","LOC",     "APW","LOC",     "APW"]
atocc[18] = [ 2.0,  0.0,       6.0,  0.0,       0.0] 
ptnl[18]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
idmd[18]  = [  0,    1,         0,    1,         0 ] 

lmb[19]   = 2
ntle[19]  = [3,3,1]
augm[19]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[19] = [ 2.0,  1.0,  0.0,       6.0,  0.0,  0.0,       0.0]
ptnl[19]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8]
idmd[19]  = [  0,    0,    1,         0,    0,    1,         0]

lmb[20]   = 2
ntle[20]  = [3,3,1]
augm[20]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[20] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0]
ptnl[20]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8]
idmd[20]  = [  0,    0,    1,         0,    0,    1,         0]

lmb[21]   = 4
ntle[21]  = [3,3,2,1,1]
augm[21]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[21] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,      0.0,      0.0]
ptnl[21]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8]
idmd[21]  = [  0,    0,    1,         0,    0,    1,         0,    1,        0,        0]  


lmb[22]   = 4                                                                              
ntle[22]  = [3,3,2,1,1]                                                                    
augm[22]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[22] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       2.0,  0.0,      0.0,      0.0] 
ptnl[22]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[22]  = [  0,    0,    1,         0,    0,    1,         0,    1,        0,        0]  


lmb[23]   = 4                                                                              
ntle[23]  = [3,3,2,1,1]                                                                    
augm[23]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[23] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       3.0,  0.0,      0.0,      0.0] 
ptnl[23]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[23]  = [  0,    0,    1,         0,    0,    1,         0,    1,        0,        0]  

lmb[24]   = 4                                                                        
ntle[24]  = [2,3,2,1,1]                                                              
augm[24]  = ["APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[24] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       4.0,  0.0,      0.0,      0.0] 
ptnl[24]  = [ 4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[24]  = [  0,    1,         0,    0,    1,         0,    1,        0,        0]  


lmb[25]   = 4                                                                        
ntle[25]  = [2,3,2,1,1]                                                              
augm[25]  = ["APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[25] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       5.0,  0.0,      0.0,      0.0] 
ptnl[25]  = [ 4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[25]  = [  0,    1,         0,    0,    1,         0,    1,        0,        0]  

lmb[26]   = 4                                                                        
ntle[26]  = [2,3,2,1,1]                                                              
augm[26]  = ["APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[26] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       6.0,  0.0,      0.0,      0.0] 
ptnl[26]  = [ 4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[26]  = [  0,    1,         0,    0,    1,         0,    1,        0,        0]  


lmb[27]   = 4                                                                        
ntle[27]  = [2,3,2,1,1]                                                              
augm[27]  = ["APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[27] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       7.0,  0.0,      0.0,      0.0] 
ptnl[27]  = [ 4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[27]  = [  0,    1,         0,    0,    1,         0,    1,        0,        0]  


lmb[28]   = 4                                                                  
ntle[28]  = [2,2,2,1,1]                                                        
augm[28]  = ["APW","LOC",     "APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[28] = [ 2.0,  0.0,       0.0,  0.0,       8.0,  0.0,      0.0,      0.0] 
ptnl[28]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[28]  = [  0,    1,         0,    1,         0,    1,        0,        0]  


lmb[29]   = 4                                                                  
ntle[29]  = [2,2,2,1,1]                                                        
augm[29]  = ["APW","LOC",     "APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[29] = [ 2.0,  0.0,       0.0,  0.0,       9.0,  0.0,      0.0,      0.0] 
ptnl[29]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[29]  = [  0,    1,         0,    1,         0,    1,        0,        0]  

lmb[30]   = 4                                                                  
ntle[30]  = [2,2,2,1,1]                                                        
augm[30]  = ["APW","LOC",     "APW","LOC",     "APW","LOC",    "APW",    "APW"]
atocc[30] = [ 2.0,  0.0,       0.0,  0.0,      10.0,  0.0,      0.0,      0.0] 
ptnl[30]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
idmd[30]  = [  0,    1,         0,    1,         0,    1,        0,        0]  


lmb[31]   = 2
ntle[31]  = [2,2,3]
augm[31]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[31] = [ 2.0,  0.0,       1.0,  0.0,       10.0, 0.0, 0.0]
ptnl[31]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]
idmd[31]  = [  0,    1,         0,    1,         0,    0,   1]

lmb[32]   = 2
ntle[32]  = [2,2,3]
augm[32]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[32] = [ 2.0,  0.0,       2.0,  0.0,       10.0, 0.0, 0.0]
ptnl[32]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]
idmd[32]  = [  0,    1,         0,    1,         0,    0,   1]


lmb[33]   = 2
ntle[33]  = [2,2,3]
augm[33]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[33] = [ 2.0,  0.0,       3.0,  0.0,       10.0, 0.0, 0.0]
ptnl[33]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]
idmd[33]  = [  0,    1,         0,    1,         0,    0,   1]


lmb[34]   = 2
ntle[34]  = [2,2,3]
augm[34]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[34] = [ 2.0,  0.0,       4.0,  0.0,       10.0, 0.0, 0.0]
ptnl[34]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]
idmd[34]  = [  0,    1,         0,    1,         0,    0,   1]   


lmb[35]   = 2                                                    
ntle[35]  = [2,2,3]                                              
augm[35]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[35] = [ 2.0,  0.0,       5.0,  0.0,       10.0, 0.0, 0.0]  
ptnl[35]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]  
idmd[35]  = [  0,    1,         0,    1,         0,    0,   1]   

lmb[36]   = 2                                                    
ntle[36]  = [2,2,3]                                              
augm[36]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[36] = [ 2.0,  0.0,       6.0,  0.0,       10.0, 0.0, 0.0]  
ptnl[36]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]  
idmd[36]  = [  0,    1,         0,    1,         0,    0,   1]   

lmb[37]   = 2
ntle[37]  = [3,3,1]
augm[37]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[37] = [ 2.0,  1.0,  0.0,       6.0,  0.0,  0.0,       0.0 ]
ptnl[37]  = [ 4.8,  5.8,  6.8,       4.8,  5.8,  6.8,       4.8]
idmd[37]  = [  0,    0,    1,         0,    0,    1,         0]  

lmb[38]   = 2                                                    
ntle[38]  = [3,3,1]                                              
augm[38]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[38] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0 ]
ptnl[38]  = [ 4.8,  5.8,  6.8,       4.8,  5.8,  6.8,       4.8] 
idmd[38]  = [  0,    0,    1,         0,    0,    1,         0]  

lmb[39]   = 4
ntle[39]  = [3,3,2,1,1]
augm[39]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
atocc[39] = [ 2.0,  2.0,  0.0,      6.0,  0.0,  0.0,        1.0,  0.0,        0.0,        0.0]
ptnl[39]  = [ 4.8,  5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]
idmd[39]  = [  0,    0,    1,        0,    0,    1,          0,    1,          0,          0]   

lmb[40]   = 4                                                                                   
ntle[40]  = [3,3,2,1,1]                                                                         
augm[40]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
atocc[40] = [ 2.0,  2.0,  0.0,      6.0,  0.0,  0.0,        2.0,  0.0,        0.0,        0.0]  
ptnl[40]  = [ 4.8,  5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]  
idmd[40]  = [  0,    0,    1,        0,    0,    1,          0,    1,          0,          0]   

lmb[41]   = 4                                                                                   
ntle[41]  = [3,3,2,1,1]                                                                         
augm[41]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
atocc[41] = [ 2.0,  1.0,  0.0,      6.0,  0.0,  0.0,        4.0,  0.0,        0.0,        0.0]  
ptnl[41]  = [ 4.8,  5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]  
idmd[41]  = [  0,    0,    1,        0,    0,    1,          0,    1,          0,          0]   

lmb[42]   = 4                                                                                   
ntle[42]  = [3,3,2,1,1]                                                                         
augm[42]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
atocc[42] = [ 2.0,  1.0,  0.0,      6.0,  0.0,  0.0,        5.0,  0.0,        0.0,        0.0]  
ptnl[42]  = [ 4.8,  5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]  
idmd[42]  = [  0,    0,    1,        0,    0,    1,          0,    1,          0,          0]   

lmb[43]   = 4                                                                            
ntle[43]  = [2,3,2,1,1]                                                                   
augm[43]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
atocc[43] = [ 2.0,  0.0,      6.0,  0.0,  0.0,        5.0,  0.0,        0.0,        0.0]  
ptnl[43]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]  
idmd[43]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]   

lmb[44]   = 4                                                                             
ntle[44]  = [2,3,2,1,1]                                                                   
augm[44]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
atocc[44] = [ 1.0,  0.0,      6.0,  0.0,  0.0,        7.0,  0.0,        0.0,        0.0] 
ptnl[44]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
idmd[44]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

lmb[45]   = 4                                                                            
ntle[45]  = [2,3,2,1,1]                                                                  
augm[45]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
atocc[45] = [ 1.0,  0.0,      6.0,  0.0,  0.0,        8.0,  0.0,        0.0,        0.0] 
ptnl[45]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
idmd[45]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

lmb[46]   = 4                                                                            
ntle[46]  = [2,3,2,1,1]                                                                  
augm[46]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
atocc[46] = [ 0.0,  0.0,      6.0,  0.0,  0.0,       10.0,  0.0,        0.0,        0.0] 
ptnl[46]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
idmd[46]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

lmb[47]   = 4                                                                            
ntle[47]  = [2,3,2,1,1]                                                                  
augm[47]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
atocc[47] = [ 1.0,  0.0,      6.0,  0.0,  0.0,       10.0,  0.0,        0.0,        0.0] 
ptnl[47]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
idmd[47]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

lmb[48]   = 4                                                                            
ntle[48]  = [2,3,2,1,1]                                                                  
augm[48]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
atocc[48] = [ 2.0,  0.0,      6.0,  0.0,  0.0,       10.0,  0.0,        0.0,        0.0] 
ptnl[48]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
idmd[48]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

lmb[49]   = 2
ntle[49]  = [2,2,3]
augm[49]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[49] = [ 2.0,  0.0,       1.0,  0.0,      10.0,  0.0,  0.0 ]
ptnl[49]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
idmd[49]  = [  0,    1,         0,    1,         0,    0,    1  ]

lmb[50]   = 2                                                    
ntle[50]  = [2,2,3]                                              
augm[50]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[50] = [ 2.0,  0.0,       2.0,  0.0,      10.0,  0.0,  0.0 ]
ptnl[50]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
idmd[50]  = [  0,    1,         0,    1,         0,    0,    1  ]

lmb[51]   = 2                                                    
ntle[51]  = [2,2,3]                                              
augm[51]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[51] = [ 2.0,  0.0,       3.0,  0.0,      10.0,  0.0,  0.0 ]
ptnl[51]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
idmd[51]  = [  0,    1,         0,    1,         0,    0,    1  ]

lmb[52]   = 2                                                    
ntle[52]  = [2,2,3]                                              
augm[52]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[52] = [ 2.0,  0.0,       4.0,  0.0,      10.0,  0.0,  0.0 ]
ptnl[52]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
idmd[52]  = [  0,    1,         0,    1,         0,    0,    1  ]

lmb[53]   = 2                                                    
ntle[53]  = [2,2,3]                                              
augm[53]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[53] = [ 2.0,  0.0,       5.0,  0.0,      10.0,  0.0,  0.0 ]
ptnl[53]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
idmd[53]  = [  0,    1,         0,    1,         0,    0,    1  ]

lmb[54]   = 2                                                    
ntle[54]  = [2,2,3]                                              
augm[54]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
atocc[54] = [ 2.0,  0.0,       6.0,  0.0,      10.0,  0.0,  0.0 ]
ptnl[54]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
idmd[54]  = [  0,    1,         0,    1,         0,    0,    1  ]

lmb[55]   = 2
ntle[55]  = [3,3,1]
augm[55]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[55] = [ 2.0,  1.0,  0.0,       6.0,  0.0,  0.0,       0.0 ]
ptnl[55]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8 ]
idmd[55]  = [  0,    0,    1,         0,    0,    1,         0  ]

lmb[56]   = 2                                                    
ntle[56]  = [3,3,1]                                              
augm[56]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[56] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0 ]
ptnl[56]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8 ]
idmd[56]  = [  0,    0,    1,         0,    0,    1,         0  ]


lmb[57]   = 6
ntle[57]  = [3,3,2,2,1,1,1]
augm[57]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",      "APW",     "APW",     "APW"]
atocc[57] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,       0.0,  0.0,        0.0,       0.0,       0.0 ]
ptnl[57]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,        5.8,       6.8,       7.8 ]
idmd[57]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,          0,         0,         0  ]

lmb[58]   = 6
ntle[58]  = [3,3,2,2,1,1,1]
augm[58]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[58] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,       1.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[58]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[58]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[59]   = 6                                                                                                            
ntle[59]  = [3,3,2,2,1,1,1]                                                                                              
augm[59]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[59] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       3.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[59]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[59]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[60]   = 6                                                                                                            
ntle[60]  = [3,3,2,2,1,1,1]                                                                                              
augm[60]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[60] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       4.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[60]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[60]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[61]   = 6                                                                                                            
ntle[61]  = [3,3,2,2,1,1,1]                                                                                              
augm[61]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[61] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       5.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[61]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[61]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[62]   = 6                                                                                                            
ntle[62]  = [3,3,2,2,1,1,1]                                                                                              
augm[62]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[62] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       6.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[62]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[62]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[63]   = 6                                                                                                            
ntle[63]  = [3,3,2,2,1,1,1]                                                                                              
augm[63]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[63] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       7.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[63]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[63]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[64]   = 6                                                                                                            
ntle[64]  = [3,3,2,2,1,1,1]                                                                                              
augm[64]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[64] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,       7.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[64]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[64]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[65]   = 6                                                                                                            
ntle[65]  = [3,3,2,2,1,1,1]                                                                                              
augm[65]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[65] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       9.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[65]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[65]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[66]   = 6                                                                                                            
ntle[66]  = [3,3,2,2,1,1,1]                                                                                              
augm[66]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[66] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      10.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[66]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[66]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[67]   = 6                                                                                                            
ntle[67]  = [3,3,2,2,1,1,1]                                                                                              
augm[67]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[67] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      11.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[67]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[67]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[68]   = 6                                                                                                            
ntle[68]  = [3,3,2,2,1,1,1]                                                                                              
augm[68]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[68] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      12.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[68]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[68]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[69]   = 6                                                                                                            
ntle[69]  = [3,3,2,2,1,1,1]                                                                                              
augm[69]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[69] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      13.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[69]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[69]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[70]   = 6                                                                                                            
ntle[70]  = [3,3,2,2,1,1,1]                                                                                              
augm[70]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[70] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      14.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[70]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[70]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[71]   = 6                                                                                                            
ntle[71]  = [3,3,2,2,1,1,1]                                                                                              
augm[71]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
atocc[71] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,      14.0,  0.0,       0.0,       0.0,       0.0 ]
ptnl[71]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
idmd[71]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

lmb[72]   = 4
ntle[72]  = [3,3,2,3,1]
augm[72]  = ["LOC","APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[72] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,      2.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
ptnl[72]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
idmd[72]  = [  0,    0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

lmb[73]   = 4                                                                                           
ntle[73]  = [3,3,2,3,1]                                                                                 
augm[73]  = ["LOC","APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[73] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,      3.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
ptnl[73]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
idmd[73]  = [  0,    0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

lmb[74]   = 4                                                                                           
ntle[74]  = [3,3,2,3,1]                                                                                 
augm[74]  = ["LOC","APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[74] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,      4.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
ptnl[74]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
idmd[74]  = [  0,    0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

lmb[75]   = 4                                                                                           
ntle[75]  = [2,3,2,3,1]                                                                                
augm[75]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[75] = [ 2.0,  0.0,       6.0,  0.0,  0.0,      5.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
ptnl[75]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
idmd[75]  = [  0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

lmb[76]   = 4                                                                                     
ntle[76]  = [2,3,2,3,1]                                                                           
augm[76]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
atocc[76] = [ 2.0,  0.0,       6.0,  0.0,  0.0,      6.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
ptnl[76]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
idmd[76]  = [  0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

lmb[77]   = 4                                                                                     
ntle[77]  = [2,3,2,2,1]                                                                           
augm[77]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "APW","LOC",     "APW"]
atocc[77] = [ 2.0,  0.0,       6.0,  0.0,  0.0,      7.0,  0.0,       0.0,  0.0,       0.0 ]
ptnl[77]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       5.8,  6.8,       5.8 ]
idmd[77]  = [  0,    1,         0,    0,    1,        0,    1,         0,    1,         0  ]

lmb[78]   = 4                                                                               
ntle[78]  = [2,3,2,2,1]                                                                     
augm[78]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "APW","LOC",     "APW"]
atocc[78] = [ 1.0,  0.0,       6.0,  0.0,  0.0,      9.0,  0.0,       0.0,  0.0,       0.0 ]
ptnl[78]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       5.8,  6.8,       5.8 ]
idmd[78]  = [  0,    1,         0,    0,    1,        0,    1,         0,    1,         0  ]


lmb[79]   = 4                                                                               
ntle[79]  = [2,3,2,2,1]                                                                     
augm[79]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "APW","LOC",     "APW"]
atocc[79] = [ 1.0,  0.0,       6.0,  0.0,  0.0,     10.0,  0.0,       0.0,  0.0,       0.0 ]
ptnl[79]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       5.8,  6.8,       5.8 ]
idmd[79]  = [  0,    1,         0,    0,    1,        0,    1,         0,    1,         0  ]


lmb[80]   = 4                                                                               
ntle[80]  = [2,3,2,2,1]                                                                     
augm[80]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "APW","LOC",     "APW"]
atocc[80] = [ 2.0,  0.0,       6.0,  0.0,  0.0,     10.0,  0.0,       0.0,  0.0,       0.0 ]
ptnl[80]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       5.8,  6.8,       5.8 ]
idmd[80]  = [  0,    1,         0,    0,    1,        0,    1,         0,    1,         0  ]


lmb[81]   = 2
ntle[81]  = [2,2,3]
augm[81]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
atocc[81] = [ 2.0,  0.0,     1.0,  0.0,       10.0, 0.0,  0.0 ]
ptnl[81]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8 ]
idmd[81]  = [  0,    1,       0,    1,        0,    0,    1,  ]


lmb[82]   = 2
ntle[82]  = [2,2,3]
augm[82]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
atocc[82] = [ 2.0,  0.0,     2.0,  0.0,       10.0, 0.0,  0.0 ]
ptnl[82]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8 ]
idmd[82]  = [  0,    1,       0,    1,        0,    0,    1   ]


lmb[83]   = 2
ntle[83]  = [2,2,3]
augm[83]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
atocc[83] = [ 2.0,  0.0,     3.0,  0.0,       10.0, 0.0,  0.0 ]
ptnl[83]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8 ]
idmd[83]  = [  0,    1,       0,    1,        0,    0,    1   ]

lmb[84]   = 2
ntle[84]  = [2,2,3]
augm[84]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
atocc[84] = [ 2.0,  0.0,     4.0,  0.0,       10.0, 0.0,  0.0 ]
ptnl[84]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8 ]
idmd[84]  = [  0,    1,       0,    1,        0,    0,    1 ]

lmb[85]   = 2
ntle[85]  = [2,2,3]
augm[85]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
atocc[85] = [ 2.0,  0.0,     5.0,  0.0,       10.0, 0.0,  0.0]
ptnl[85]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8]
idmd[85]  = [  0,    1,       0,    1,        0,    0,    1]

lmb[86]   = 2
ntle[86]  = [2,2,3]
augm[86]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
atocc[86] = [ 2.0,  0.0,     6.0,  0.0,       10.0, 0.0,  0.0]
ptnl[86]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8]
idmd[86]  = [  0,    1,       0,    1,        0,    0,    1]


lmb[87]   = 2
ntle[87]  = [3,3,1]
augm[87]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW"]
atocc[87] = [ 2.0,  1.0,  0.0,      6.0,  0.0,  0.0,        0.0]
ptnl[87]  = [ 6.8,  7.8,  8.8,      6.8,  7.8,  8.8,        6.8]
idmd[87]  = [  0,    0,    1,        0,    0,    1,          0 ] 


lmb[88]   = 2                                                    
ntle[88]  = [3,3,1]                                              
augm[88]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW"]
atocc[88] = [ 2.0,  2.0,  0.0,      6.0,  0.0,  0.0,        0.0] 
ptnl[88]  = [ 6.8,  7.8,  8.8,      6.8,  7.8,  8.8,        6.8] 
idmd[88]  = [  0,    0,    1,        0,    0,    1,          0 ] 

lmb[89]   = 6
ntle[89]  = [3,3,2,2,1,1,1]
augm[89]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW","LOC",   "APW",    "APW",    "APW"]
atocc[89] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,      0.0,  0.0,     0.0,      0.0,      0.0 ]
ptnl[89]  = [ 6.8,  7.8,  8.8,       6.8,  7.8,  8.8,       6.8,  7.8,      5.8,  6.8,     5.8,      6.8,      7.8 ]
idmd[89]  = [  0,    0,    1,         0,    0,    1,         0,    1,        0,    1,       0,        0,        0  ]


lmb[90]   = 6
ntle[90]  = [3,3,2,2,1,1,1]
augm[90]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[90] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     2.0,  0.0,     0.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[90]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[90]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[91]   = 6                                                                                                   
ntle[91]  = [3,3,2,2,1,1,1]                                                                                     
augm[91]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[91] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     2.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[91]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[91]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[92]   = 6                                                                                                   
ntle[92]  = [3,3,2,2,1,1,1]                                                                                     
augm[92]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[92] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     3.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[92]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[92]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[93]   = 6                                                                                                   
ntle[93]  = [3,3,2,2,1,1,1]                                                                                     
augm[93]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[93] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     4.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[93]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[93]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[94]   = 6                                                                                                   
ntle[94]  = [3,3,2,2,1,1,1]                                                                                     
augm[94]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[94] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     5.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[94]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[94]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[95]   = 6                                                                                                   
ntle[95]  = [3,3,2,2,1,1,1]                                                                                     
augm[95]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[95] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     6.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[95]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[95]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[96]   = 6                                                                                                   
ntle[96]  = [3,3,2,2,1,1,1]                                                                                     
augm[96]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[96] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     7.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[96]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[96]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[97]   = 6                                                                                                   
ntle[97]  = [3,3,2,2,1,1,1]                                                                                     
augm[97]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[97] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     8.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[97]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[97]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[98]   = 6                                                                                                   
ntle[98]  = [3,3,2,2,1,1,1]                                                                                     
augm[98]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[98] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     9.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[98]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[98]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]

lmb[99]   = 6                                                                                                   
ntle[99]  = [3,3,2,2,1,1,1]                                                                                     
augm[99]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[99] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    10.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[99]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[99]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]

lmb[100]   = 6                                                                                                   
ntle[100]  = [3,3,2,2,1,1,1]                                                                                     
augm[100]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[100] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    11.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[100]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[100]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]

lmb[101]   = 6                                                                                                   
ntle[101]  = [3,3,2,2,1,1,1]                                                                                     
augm[101]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[101] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    12.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[101]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[101]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]

lmb[102]   = 6                                                                                                   
ntle[102]  = [3,3,2,2,1,1,1]                                                                                     
augm[102]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[102] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    13.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[102]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[102]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


lmb[103]   = 6                                                                                                   
ntle[103]  = [3,3,2,2,1,1,1]                                                                                     
augm[103]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
atocc[103] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    14.0,  0.0,      0.0,      0.0,      0.0 ]
ptnl[103]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
idmd[103]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]



for ii in range(1,104):
    # check lmb
    # if (len(ntle[ii]) != lmb[ii]+1):
    #     print(ii, 'lmb mismach\n')
    
    # # check the number of atomic orbitals
    # if (sum(ntle[ii]) !=len(augm[ii])):
    #     print(ii, 'augm length wrong')
    # if (sum(ntle[ii]) !=len(atocc[ii])):
    #     print(ii, 'atocc length wrong', sum(ntle[ii]), len(atocc[ii]), atocc[ii])
    # if (sum(ntle[ii]) !=len(ptnl[ii])):
    #     print(ii, 'ptnl length wrong')
    # if (sum(ntle[ii]) !=len(idmd[ii])):
    #     print(ii, 'idmd length wrong')            
    # # check total charge
    # ncore=0
    # nvalence=0
    # for jj in range(lmb[ii]+1):
    #     kstart=int(np.sum(ntle[ii][:jj]))
    #     # print(ii, jj, kstart)        
    #     for kk in range(kstart, kstart+ntle[ii][jj]):
    #         #core
            
    #         if (kk == kstart):
    #             ncore=ncore+(int(ptnl[ii][kk])-(jj+1))*(4*jj+2)
    # nvalence=np.sum(atocc[ii])
    # ntot=ncore+nvalence
    # if (ntot != ii):
    #     print(ii, ntot,ncore,nvalence,'total change wrong')                


    # #   check valence orbital
    # print(ii, end='  ')
    # for jj in range(lmb[ii]+1):
    #     kstart=int(np.sum(ntle[ii][:jj]))
    #     # print(ii, jj, kstart)        
    #     for kk in range(kstart, kstart+ntle[ii][jj]):    
    #         if((augm[ii][kk] == "APW") & (atocc[ii][kk]>0)):
    #             if (jj==0):
    #                 ll='s'
    #             if (jj==1):
    #                 ll='p'
    #             if (jj==2):
    #                 ll='d'
    #             if (jj==3):
    #                 ll='f'                    
    #             print(int(ptnl[ii][kk]), ll, atocc[ii][kk], end='         ')
    # print('\n')


    #   check APW: only 1 APW

    # for jj in range(lmb[ii]+1):
    #     kstart=int(np.sum(ntle[ii][:jj]))
    #     # print(ii, jj, kstart)
    #     cnt=0
    #     for kk in range(kstart, kstart+ntle[ii][jj]):    
    #         if(augm[ii][kk] == "APW"):
    #             cnt=cnt+1
    #     if (cnt !=1):
    #         print(ii, jj, augm[ii][kstart:(kstart+ntle[ii][jj])]=='APW', 'more than 1 APW')

    #  check idmd

    for jj in range(lmb[ii]+1):
        kstart=int(np.sum(ntle[ii][:jj]))
        for kk in range(kstart, kstart+ntle[ii][jj]):            
            # print(ii, jj, kstart)
            if (kk ==kstart):
                if (idmd[ii][kk]!=0):
                    print(ii, jj, kk, 'idmd wrong')
            else:
                if (augm[ii][kk-1] =='LOC'):
                    if (idmd[ii][kk]!=0):
                        print(ii, jj, kk, 'idmd wrong')
                elif (augm[ii][kk-1] =='APW'):                    
                    if (idmd[ii][kk]!=1):
                        print(ii, jj, kk, 'idmd wrong')
                else:
                    print(ii, jj, kk, 'idmd wrong')            

