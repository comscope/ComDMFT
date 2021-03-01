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

element_lmb[1]   = 1
element_ntle[1]  = [2,1]
element_augm[1]  = ["APW","LOC",    "APW"]
element_atocc[1] = [ 1.0,  0.0,      0.0]
element_ptnl[1]  = [ 1.8,  2.8,      2.8]
element_idmd[1]  = [  0,    1,        0]


element_lmb[2]   = 1
element_ntle[2]  = [2,1]
element_augm[2]  = ["APW","LOC",    "APW"]
element_atocc[2] = [ 2.0,  0.0,      0.0]
element_ptnl[2]  = [ 1.8,  2.8,      2.8]
element_idmd[2]  = [  0,    1,        0]


element_lmb[3]   = 2
element_ntle[3]  = [3,2,1]
element_augm[3]  = ["LOC","APW","LOC",    "APW","LOC",     "APW"]
element_atocc[3] = [ 2.0,  1.0,  0.0,      0.0,  0.0,       0.0]
element_ptnl[3]  = [ 1.8,  2.8,  3.8,      2.8,  3.8,       3.8]
element_idmd[3]  = [  0,    0,    1,        0,     1,        0]  


element_lmb[4]   = 2
element_ntle[4]  = [2,2,1]
element_augm[4]  = ["APW","LOC",    "APW","LOC",     "APW"]
element_atocc[4] = [ 2.0,  0.0,      0.0,  0.0,       0.0]
element_ptnl[4]  = [ 2.8,  3.8,      2.8,  3.8,       3.8]
element_idmd[4]  = [  0,    1,        0,     1,        0 ] 

element_lmb[5]   = 2                                       
element_ntle[5]  = [2,2,1]                                 
element_augm[5]  = ["APW","LOC",    "APW","LOC",     "APW"]
element_atocc[5] = [ 2.0,  0.0,      1.0,  0.0,       0.0] 
element_ptnl[5]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
element_idmd[5]  = [  0,    1,        0,     1,        0 ] 


element_lmb[6]   = 2                                       
element_ntle[6]  = [2,2,1]                                 
element_augm[6]  = ["APW","LOC",    "APW","LOC",     "APW"]
element_atocc[6] = [ 2.0,  0.0,      2.0,  0.0,       0.0] 
element_ptnl[6]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
element_idmd[6]  = [  0,    1,        0,     1,        0 ] 


element_lmb[7]   = 2                                       
element_ntle[7]  = [2,2,1]                                 
element_augm[7]  = ["APW","LOC",    "APW","LOC",     "APW"]
element_atocc[7] = [ 2.0,  0.0,      3.0,  0.0,       0.0] 
element_ptnl[7]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
element_idmd[7]  = [  0,    1,        0,     1,        0 ] 


element_lmb[8]   = 2                                       
element_ntle[8]  = [2,2,1]                                 
element_augm[8]  = ["APW","LOC",    "APW","LOC",     "APW"]
element_atocc[8] = [ 2.0,  0.0,      4.0,  0.0,       0.0] 
element_ptnl[8]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
element_idmd[8]  = [  0,    1,        0,     1,        0 ] 


element_lmb[9]   = 2                                       
element_ntle[9]  = [2,2,1]                                 
element_augm[9]  = ["APW","LOC",    "APW","LOC",     "APW"]
element_atocc[9] = [ 2.0,  0.0,      5.0,  0.0,       0.0] 
element_ptnl[9]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
element_idmd[9]  = [  0,    1,        0,     1,        0 ] 


element_lmb[10]   = 2                                       
element_ntle[10]  = [2,2,1]                                 
element_augm[10]  = ["APW","LOC",    "APW","LOC",     "APW"]
element_atocc[10] = [ 2.0,  0.0,      6.0,  0.0,       0.0] 
element_ptnl[10]  = [ 2.8,  3.8,      2.8,  3.8,       3.8] 
element_idmd[10]  = [  0,    1,        0,     1,        0 ] 

element_lmb[11]   = 2
element_ntle[11]  = [3,3,1]
element_augm[11]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[11] = [ 2.0,  1.0,  0.0,       6.0,  0.0,  0.0,       0.0]
element_ptnl[11]  = [ 2.8,  3.8,  4.8,       2.8,  3.8,  4.8,       3.8]
element_idmd[11]  = [  0,    0,    1,         0,    0,    1,         0] 


element_lmb[12]   = 2                                                    
element_ntle[12]  = [2,3,1]                                              
element_augm[12]  = ["APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[12] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       0.0] 
element_ptnl[12]  = [ 3.8,  4.8,       2.8,  3.8,  4.8,       3.8] 
element_idmd[12]  = [  0,    1,         0,    0,    1,         0 ] 


element_lmb[13]   = 2                                              
element_ntle[13]  = [2,2,1]                                        
element_augm[13]  = ["APW","LOC",     "APW","LOC",     "APW"]
element_atocc[13] = [ 2.0,  0.0,       1.0,  0.0,       0.0] 
element_ptnl[13]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
element_idmd[13]  = [  0,    1,         0,    1,         0 ] 


element_lmb[14]   = 2                                              
element_ntle[14]  = [2,2,1]                                        
element_augm[14]  = ["APW","LOC",     "APW","LOC",     "APW"]
element_atocc[14] = [ 2.0,  0.0,       2.0,  0.0,       0.0] 
element_ptnl[14]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
element_idmd[14]  = [  0,    1,         0,    1,         0] 


element_lmb[15]   = 2                                              
element_ntle[15]  = [2,2,1]                                        
element_augm[15]  = ["APW","LOC",     "APW","LOC",     "APW"]
element_atocc[15] = [ 2.0,  0.0,       3.0,  0.0,       0.0] 
element_ptnl[15]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
element_idmd[15]  = [  0,    1,         0,    1,         0 ] 

element_lmb[16]   = 2                                              
element_ntle[16]  = [2,2,1]                                        
element_augm[16]  = ["APW","LOC",     "APW","LOC",     "APW"]
element_atocc[16] = [ 2.0,  0.0,       4.0,  0.0,       0.0] 
element_ptnl[16]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
element_idmd[16]  = [  0,    1,         0,    1,         0] 


element_lmb[17]   = 2                                              
element_ntle[17]  = [2,2,1]                                        
element_augm[17]  = ["APW","LOC",     "APW","LOC",     "APW"]
element_atocc[17] = [ 2.0,  0.0,       5.0,  0.0,       0.0] 
element_ptnl[17]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
element_idmd[17]  = [  0,    1,         0,    1,         0] 


element_lmb[18]   = 2                                              
element_ntle[18]  = [2,2,1]                                        
element_augm[18]  = ["APW","LOC",     "APW","LOC",     "APW"]
element_atocc[18] = [ 2.0,  0.0,       6.0,  0.0,       0.0] 
element_ptnl[18]  = [ 3.8,  4.8,       3.8,  4.8,       3.8] 
element_idmd[18]  = [  0,    1,         0,    1,         0 ] 

element_lmb[19]   = 2
element_ntle[19]  = [3,3,1]
element_augm[19]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[19] = [ 2.0,  1.0,  0.0,       6.0,  0.0,  0.0,       0.0]
element_ptnl[19]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8]
element_idmd[19]  = [  0,    0,    1,         0,    0,    1,         0]

element_lmb[20]   = 2
element_ntle[20]  = [3,3,1]
element_augm[20]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[20] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0]
element_ptnl[20]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8]
element_idmd[20]  = [  0,    0,    1,         0,    0,    1,         0]

element_lmb[21]   = 4
element_ntle[21]  = [3,3,2,1,1]
element_augm[21]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[21] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,      0.0,      0.0]
element_ptnl[21]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8]
element_idmd[21]  = [  0,    0,    1,         0,    0,    1,         0,    1,        0,        0]  


element_lmb[22]   = 4                                                                              
element_ntle[22]  = [3,3,2,1,1]                                                                    
element_augm[22]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[22] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       2.0,  0.0,      0.0,      0.0] 
element_ptnl[22]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[22]  = [  0,    0,    1,         0,    0,    1,         0,    1,        0,        0]  


element_lmb[23]   = 4                                                                              
element_ntle[23]  = [3,3,2,1,1]                                                                    
element_augm[23]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[23] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       3.0,  0.0,      0.0,      0.0] 
element_ptnl[23]  = [ 3.8,  4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[23]  = [  0,    0,    1,         0,    0,    1,         0,    1,        0,        0]  

element_lmb[24]   = 4                                                                        
element_ntle[24]  = [2,3,2,1,1]                                                              
element_augm[24]  = ["APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[24] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       4.0,  0.0,      0.0,      0.0] 
element_ptnl[24]  = [ 4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[24]  = [  0,    1,         0,    0,    1,         0,    1,        0,        0]  


element_lmb[25]   = 4                                                                        
element_ntle[25]  = [2,3,2,1,1]                                                              
element_augm[25]  = ["APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[25] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       5.0,  0.0,      0.0,      0.0] 
element_ptnl[25]  = [ 4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[25]  = [  0,    1,         0,    0,    1,         0,    1,        0,        0]  

element_lmb[26]   = 4                                                                        
element_ntle[26]  = [2,3,2,1,1]                                                              
element_augm[26]  = ["APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[26] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       6.0,  0.0,      0.0,      0.0] 
element_ptnl[26]  = [ 4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[26]  = [  0,    1,         0,    0,    1,         0,    1,        0,        0]  


element_lmb[27]   = 4                                                                        
element_ntle[27]  = [2,3,2,1,1]                                                              
element_augm[27]  = ["APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[27] = [ 2.0,  0.0,       6.0,  0.0,  0.0,       7.0,  0.0,      0.0,      0.0] 
element_ptnl[27]  = [ 4.8,  5.8,       3.8,  4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[27]  = [  0,    1,         0,    0,    1,         0,    1,        0,        0]  


element_lmb[28]   = 4                                                                  
element_ntle[28]  = [2,2,2,1,1]                                                        
element_augm[28]  = ["APW","LOC",     "APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[28] = [ 2.0,  0.0,       0.0,  0.0,       8.0,  0.0,      0.0,      0.0] 
element_ptnl[28]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[28]  = [  0,    1,         0,    1,         0,    1,        0,        0]  


element_lmb[29]   = 4                                                                  
element_ntle[29]  = [2,2,2,1,1]                                                        
element_augm[29]  = ["APW","LOC",     "APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[29] = [ 2.0,  0.0,       0.0,  0.0,       9.0,  0.0,      0.0,      0.0] 
element_ptnl[29]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[29]  = [  0,    1,         0,    1,         0,    1,        0,        0]  

element_lmb[30]   = 4                                                                  
element_ntle[30]  = [2,2,2,1,1]                                                        
element_augm[30]  = ["APW","LOC",     "APW","LOC",     "APW","LOC",    "APW",    "APW"]
element_atocc[30] = [ 2.0,  0.0,       0.0,  0.0,      10.0,  0.0,      0.0,      0.0] 
element_ptnl[30]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8,      4.8,      5.8] 
element_idmd[30]  = [  0,    1,         0,    1,         0,    1,        0,        0]  


element_lmb[31]   = 2
element_ntle[31]  = [2,2,3]
element_augm[31]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[31] = [ 2.0,  0.0,       1.0,  0.0,       10.0, 0.0, 0.0]
element_ptnl[31]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]
element_idmd[31]  = [  0,    1,         0,    1,         0,    0,   1]

element_lmb[32]   = 2
element_ntle[32]  = [2,2,3]
element_augm[32]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[32] = [ 2.0,  0.0,       2.0,  0.0,       10.0, 0.0, 0.0]
element_ptnl[32]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]
element_idmd[32]  = [  0,    1,         0,    1,         0,    0,   1]


element_lmb[33]   = 2
element_ntle[33]  = [2,2,3]
element_augm[33]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[33] = [ 2.0,  0.0,       3.0,  0.0,       10.0, 0.0, 0.0]
element_ptnl[33]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]
element_idmd[33]  = [  0,    1,         0,    1,         0,    0,   1]


element_lmb[34]   = 2
element_ntle[34]  = [2,2,3]
element_augm[34]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[34] = [ 2.0,  0.0,       4.0,  0.0,       10.0, 0.0, 0.0]
element_ptnl[34]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]
element_idmd[34]  = [  0,    1,         0,    1,         0,    0,   1]   


element_lmb[35]   = 2                                                    
element_ntle[35]  = [2,2,3]                                              
element_augm[35]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[35] = [ 2.0,  0.0,       5.0,  0.0,       10.0, 0.0, 0.0]  
element_ptnl[35]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]  
element_idmd[35]  = [  0,    1,         0,    1,         0,    0,   1]   

element_lmb[36]   = 2                                                    
element_ntle[36]  = [2,2,3]                                              
element_augm[36]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[36] = [ 2.0,  0.0,       6.0,  0.0,       10.0, 0.0, 0.0]  
element_ptnl[36]  = [ 4.8,  5.8,       4.8,  5.8,       3.8,  4.8, 5.8]  
element_idmd[36]  = [  0,    1,         0,    1,         0,    0,   1]   

element_lmb[37]   = 2
element_ntle[37]  = [3,3,1]
element_augm[37]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[37] = [ 2.0,  1.0,  0.0,       6.0,  0.0,  0.0,       0.0 ]
element_ptnl[37]  = [ 4.8,  5.8,  6.8,       4.8,  5.8,  6.8,       4.8]
element_idmd[37]  = [  0,    0,    1,         0,    0,    1,         0]  

element_lmb[38]   = 2                                                    
element_ntle[38]  = [3,3,1]                                              
element_augm[38]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[38] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0 ]
element_ptnl[38]  = [ 4.8,  5.8,  6.8,       4.8,  5.8,  6.8,       4.8] 
element_idmd[38]  = [  0,    0,    1,         0,    0,    1,         0]  

element_lmb[39]   = 4
element_ntle[39]  = [3,3,2,1,1]
element_augm[39]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
element_atocc[39] = [ 2.0,  2.0,  0.0,      6.0,  0.0,  0.0,        1.0,  0.0,        0.0,        0.0]
element_ptnl[39]  = [ 4.8,  5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]
element_idmd[39]  = [  0,    0,    1,        0,    0,    1,          0,    1,          0,          0]   

element_lmb[40]   = 4                                                                                   
element_ntle[40]  = [3,3,2,1,1]                                                                         
element_augm[40]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
element_atocc[40] = [ 2.0,  2.0,  0.0,      6.0,  0.0,  0.0,        2.0,  0.0,        0.0,        0.0]  
element_ptnl[40]  = [ 4.8,  5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]  
element_idmd[40]  = [  0,    0,    1,        0,    0,    1,          0,    1,          0,          0]   

element_lmb[41]   = 4                                                                                   
element_ntle[41]  = [3,3,2,1,1]                                                                         
element_augm[41]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
element_atocc[41] = [ 2.0,  1.0,  0.0,      6.0,  0.0,  0.0,        4.0,  0.0,        0.0,        0.0]  
element_ptnl[41]  = [ 4.8,  5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]  
element_idmd[41]  = [  0,    0,    1,        0,    0,    1,          0,    1,          0,          0]   

element_lmb[42]   = 4                                                                                   
element_ntle[42]  = [3,3,2,1,1]                                                                         
element_augm[42]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
element_atocc[42] = [ 2.0,  1.0,  0.0,      6.0,  0.0,  0.0,        5.0,  0.0,        0.0,        0.0]  
element_ptnl[42]  = [ 4.8,  5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]  
element_idmd[42]  = [  0,    0,    1,        0,    0,    1,          0,    1,          0,          0]   

element_lmb[43]   = 4                                                                            
element_ntle[43]  = [2,3,2,1,1]                                                                   
element_augm[43]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
element_atocc[43] = [ 2.0,  0.0,      6.0,  0.0,  0.0,        5.0,  0.0,        0.0,        0.0]  
element_ptnl[43]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8]  
element_idmd[43]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]   

element_lmb[44]   = 4                                                                             
element_ntle[44]  = [2,3,2,1,1]                                                                   
element_augm[44]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"] 
element_atocc[44] = [ 1.0,  0.0,      6.0,  0.0,  0.0,        7.0,  0.0,        0.0,        0.0] 
element_ptnl[44]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
element_idmd[44]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

element_lmb[45]   = 4                                                                            
element_ntle[45]  = [2,3,2,1,1]                                                                  
element_augm[45]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
element_atocc[45] = [ 1.0,  0.0,      6.0,  0.0,  0.0,        8.0,  0.0,        0.0,        0.0] 
element_ptnl[45]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
element_idmd[45]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

element_lmb[46]   = 4                                                                            
element_ntle[46]  = [2,3,2,1,1]                                                                  
element_augm[46]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
element_atocc[46] = [ 0.0,  0.0,      6.0,  0.0,  0.0,       10.0,  0.0,        0.0,        0.0] 
element_ptnl[46]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
element_idmd[46]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

element_lmb[47]   = 4                                                                            
element_ntle[47]  = [2,3,2,1,1]                                                                  
element_augm[47]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
element_atocc[47] = [ 1.0,  0.0,      6.0,  0.0,  0.0,       10.0,  0.0,        0.0,        0.0] 
element_ptnl[47]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
element_idmd[47]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

element_lmb[48]   = 4                                                                            
element_ntle[48]  = [2,3,2,1,1]                                                                  
element_augm[48]  = ["APW","LOC",    "LOC","APW","LOC",      "APW","LOC",      "APW",      "APW"]
element_atocc[48] = [ 2.0,  0.0,      6.0,  0.0,  0.0,       10.0,  0.0,        0.0,        0.0] 
element_ptnl[48]  = [ 5.8,  6.8,      4.8,  5.8,  6.8,        4.8,  5.8,        4.8,        5.8] 
element_idmd[48]  = [  0,    1,        0,    0,    1,          0,    1,          0,          0]  

element_lmb[49]   = 2
element_ntle[49]  = [2,2,3]
element_augm[49]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[49] = [ 2.0,  0.0,       1.0,  0.0,      10.0,  0.0,  0.0 ]
element_ptnl[49]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
element_idmd[49]  = [  0,    1,         0,    1,         0,    0,    1  ]

element_lmb[50]   = 2                                                    
element_ntle[50]  = [2,2,3]                                              
element_augm[50]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[50] = [ 2.0,  0.0,       2.0,  0.0,      10.0,  0.0,  0.0 ]
element_ptnl[50]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
element_idmd[50]  = [  0,    1,         0,    1,         0,    0,    1  ]

element_lmb[51]   = 2                                                    
element_ntle[51]  = [2,2,3]                                              
element_augm[51]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[51] = [ 2.0,  0.0,       3.0,  0.0,      10.0,  0.0,  0.0 ]
element_ptnl[51]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
element_idmd[51]  = [  0,    1,         0,    1,         0,    0,    1  ]

element_lmb[52]   = 2                                                    
element_ntle[52]  = [2,2,3]                                              
element_augm[52]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[52] = [ 2.0,  0.0,       4.0,  0.0,      10.0,  0.0,  0.0 ]
element_ptnl[52]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
element_idmd[52]  = [  0,    1,         0,    1,         0,    0,    1  ]

element_lmb[53]   = 2                                                    
element_ntle[53]  = [2,2,3]                                              
element_augm[53]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[53] = [ 2.0,  0.0,       5.0,  0.0,      10.0,  0.0,  0.0 ]
element_ptnl[53]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
element_idmd[53]  = [  0,    1,         0,    1,         0,    0,    1  ]

element_lmb[54]   = 2                                                    
element_ntle[54]  = [2,2,3]                                              
element_augm[54]  = ["APW","LOC",     "APW","LOC",     "LOC","APW","LOC"]
element_atocc[54] = [ 2.0,  0.0,       6.0,  0.0,      10.0,  0.0,  0.0 ]
element_ptnl[54]  = [ 5.8,  6.8,       5.8,  6.8,       4.8,  5.8,  6.8 ]
element_idmd[54]  = [  0,    1,         0,    1,         0,    0,    1  ]

element_lmb[55]   = 2
element_ntle[55]  = [3,3,1]
element_augm[55]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[55] = [ 2.0,  1.0,  0.0,       6.0,  0.0,  0.0,       0.0 ]
element_ptnl[55]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8 ]
element_idmd[55]  = [  0,    0,    1,         0,    0,    1,         0  ]

element_lmb[56]   = 2                                                    
element_ntle[56]  = [3,3,1]                                              
element_augm[56]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[56] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0 ]
element_ptnl[56]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8 ]
element_idmd[56]  = [  0,    0,    1,         0,    0,    1,         0  ]


element_lmb[57]   = 6
element_ntle[57]  = [3,3,2,2,1,1,1]
element_augm[57]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",      "APW",     "APW",     "APW"]
element_atocc[57] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,       0.0,  0.0,        0.0,       0.0,       0.0 ]
element_ptnl[57]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,        5.8,       6.8,       7.8 ]
element_idmd[57]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,          0,         0,         0  ]

element_lmb[58]   = 6
element_ntle[58]  = [3,3,2,2,1,1,1]
element_augm[58]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[58] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,       1.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[58]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[58]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[59]   = 6                                                                                                            
element_ntle[59]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[59]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[59] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       3.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[59]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[59]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[60]   = 6                                                                                                            
element_ntle[60]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[60]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[60] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       4.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[60]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[60]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[61]   = 6                                                                                                            
element_ntle[61]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[61]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[61] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       5.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[61]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[61]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[62]   = 6                                                                                                            
element_ntle[62]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[62]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[62] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       6.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[62]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[62]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[63]   = 6                                                                                                            
element_ntle[63]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[63]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[63] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       7.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[63]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[63]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[64]   = 6                                                                                                            
element_ntle[64]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[64]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[64] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,       7.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[64]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[64]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[65]   = 6                                                                                                            
element_ntle[65]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[65]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[65] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,       9.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[65]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[65]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[66]   = 6                                                                                                            
element_ntle[66]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[66]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[66] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      10.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[66]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[66]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[67]   = 6                                                                                                            
element_ntle[67]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[67]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[67] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      11.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[67]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[67]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[68]   = 6                                                                                                            
element_ntle[68]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[68]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[68] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      12.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[68]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[68]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[69]   = 6                                                                                                            
element_ntle[69]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[69]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[69] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      13.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[69]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[69]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[70]   = 6                                                                                                            
element_ntle[70]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[70]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[70] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       0.0,  0.0,      14.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[70]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[70]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[71]   = 6                                                                                                            
element_ntle[71]  = [3,3,2,2,1,1,1]                                                                                              
element_augm[71]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",     "APW","LOC",     "APW",     "APW",     "APW"]
element_atocc[71] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,      14.0,  0.0,       0.0,       0.0,       0.0 ]
element_ptnl[71]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,       5.8,  6.8,       4.8,  5.8,       5.8,       6.8,       7.8 ]
element_idmd[71]  = [  0,    0,    1,         0,    0,    1,         0,    1,         0,    1,         0,         0,         0  ]

element_lmb[72]   = 4
element_ntle[72]  = [3,3,2,3,1]
element_augm[72]  = ["LOC","APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[72] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,      2.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
element_ptnl[72]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
element_idmd[72]  = [  0,    0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

element_lmb[73]   = 4                                                                                           
element_ntle[73]  = [3,3,2,3,1]                                                                                 
element_augm[73]  = ["LOC","APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[73] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,      3.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
element_ptnl[73]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
element_idmd[73]  = [  0,    0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

element_lmb[74]   = 4                                                                                           
element_ntle[74]  = [3,3,2,3,1]                                                                                 
element_augm[74]  = ["LOC","APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[74] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,      4.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
element_ptnl[74]  = [ 5.8,  6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
element_idmd[74]  = [  0,    0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

element_lmb[75]   = 4                                                                                           
element_ntle[75]  = [2,3,2,3,1]                                                                                
element_augm[75]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[75] = [ 2.0,  0.0,       6.0,  0.0,  0.0,      5.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
element_ptnl[75]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
element_idmd[75]  = [  0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

element_lmb[76]   = 4                                                                                     
element_ntle[76]  = [2,3,2,3,1]                                                                           
element_augm[76]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "LOC","APW","LOC",     "APW"]
element_atocc[76] = [ 2.0,  0.0,       6.0,  0.0,  0.0,      6.0,  0.0,      14.0,  0.0,  0.0,       0.0 ]
element_ptnl[76]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       4.8,  5.8,  6.8,       5.8 ]
element_idmd[76]  = [  0,    1,         0,    0,    1,        0,    1,         0,    0,    1,         0  ]

element_lmb[77]   = 4                                                                                     
element_ntle[77]  = [2,3,2,2,1]                                                                           
element_augm[77]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "APW","LOC",     "APW"]
element_atocc[77] = [ 2.0,  0.0,       6.0,  0.0,  0.0,      7.0,  0.0,       0.0,  0.0,       0.0 ]
element_ptnl[77]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       5.8,  6.8,       5.8 ]
element_idmd[77]  = [  0,    1,         0,    0,    1,        0,    1,         0,    1,         0  ]

element_lmb[78]   = 4                                                                               
element_ntle[78]  = [2,3,2,2,1]                                                                     
element_augm[78]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "APW","LOC",     "APW"]
element_atocc[78] = [ 1.0,  0.0,       6.0,  0.0,  0.0,      9.0,  0.0,       0.0,  0.0,       0.0 ]
element_ptnl[78]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       5.8,  6.8,       5.8 ]
element_idmd[78]  = [  0,    1,         0,    0,    1,        0,    1,         0,    1,         0  ]


element_lmb[79]   = 4                                                                               
element_ntle[79]  = [2,3,2,2,1]                                                                     
element_augm[79]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "APW","LOC",     "APW"]
element_atocc[79] = [ 1.0,  0.0,       6.0,  0.0,  0.0,     10.0,  0.0,       0.0,  0.0,       0.0 ]
element_ptnl[79]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       5.8,  6.8,       5.8 ]
element_idmd[79]  = [  0,    1,         0,    0,    1,        0,    1,         0,    1,         0  ]


element_lmb[80]   = 4                                                                               
element_ntle[80]  = [2,3,2,2,1]                                                                     
element_augm[80]  = ["APW","LOC",     "LOC","APW","LOC",    "APW","LOC",     "APW","LOC",     "APW"]
element_atocc[80] = [ 2.0,  0.0,       6.0,  0.0,  0.0,     10.0,  0.0,       0.0,  0.0,       0.0 ]
element_ptnl[80]  = [ 6.8,  7.8,       5.8,  6.8,  7.8,      5.8,  6.8,       5.8,  6.8,       5.8 ]
element_idmd[80]  = [  0,    1,         0,    0,    1,        0,    1,         0,    1,         0  ]


element_lmb[81]   = 2
element_ntle[81]  = [2,2,3]
element_augm[81]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
element_atocc[81] = [ 2.0,  0.0,     1.0,  0.0,       10.0, 0.0,  0.0 ]
element_ptnl[81]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8 ]
element_idmd[81]  = [  0,    1,       0,    1,        0,    0,    1,  ]


element_lmb[82]   = 2
element_ntle[82]  = [2,2,3]
element_augm[82]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
element_atocc[82] = [ 2.0,  0.0,     2.0,  0.0,       10.0, 0.0,  0.0 ]
element_ptnl[82]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8 ]
element_idmd[82]  = [  0,    1,       0,    1,        0,    0,    1   ]


element_lmb[83]   = 2
element_ntle[83]  = [2,2,3]
element_augm[83]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
element_atocc[83] = [ 2.0,  0.0,     3.0,  0.0,       10.0, 0.0,  0.0 ]
element_ptnl[83]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8 ]
element_idmd[83]  = [  0,    1,       0,    1,        0,    0,    1   ]

element_lmb[84]   = 2
element_ntle[84]  = [2,2,3]
element_augm[84]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
element_atocc[84] = [ 2.0,  0.0,     4.0,  0.0,       10.0, 0.0,  0.0 ]
element_ptnl[84]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8 ]
element_idmd[84]  = [  0,    1,       0,    1,        0,    0,    1 ]

element_lmb[85]   = 2
element_ntle[85]  = [2,2,3]
element_augm[85]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
element_atocc[85] = [ 2.0,  0.0,     5.0,  0.0,       10.0, 0.0,  0.0]
element_ptnl[85]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8]
element_idmd[85]  = [  0,    1,       0,    1,        0,    0,    1]

element_lmb[86]   = 2
element_ntle[86]  = [2,2,3]
element_augm[86]  = ["APW","LOC",   "APW","LOC",     "LOC","APW","LOC"]
element_atocc[86] = [ 2.0,  0.0,     6.0,  0.0,       10.0, 0.0,  0.0]
element_ptnl[86]  = [ 6.8,  7.8,     6.8,  7.8,       5.8,  6.8,  7.8]
element_idmd[86]  = [  0,    1,       0,    1,        0,    0,    1]


element_lmb[87]   = 2
element_ntle[87]  = [3,3,1]
element_augm[87]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW"]
element_atocc[87] = [ 2.0,  1.0,  0.0,      6.0,  0.0,  0.0,        0.0]
element_ptnl[87]  = [ 6.8,  7.8,  8.8,      6.8,  7.8,  8.8,        6.8]
element_idmd[87]  = [  0,    0,    1,        0,    0,    1,          0 ] 


element_lmb[88]   = 2                                                    
element_ntle[88]  = [3,3,1]                                              
element_augm[88]  = ["LOC","APW","LOC",    "LOC","APW","LOC",      "APW"]
element_atocc[88] = [ 2.0,  2.0,  0.0,      6.0,  0.0,  0.0,        0.0] 
element_ptnl[88]  = [ 6.8,  7.8,  8.8,      6.8,  7.8,  8.8,        6.8] 
element_idmd[88]  = [  0,    0,    1,        0,    0,    1,          0 ] 

element_lmb[89]   = 6
element_ntle[89]  = [3,3,2,2,1,1,1]
element_augm[89]  = ["LOC","APW","LOC",     "LOC","APW","LOC",     "APW","LOC",    "APW","LOC",   "APW",    "APW",    "APW"]
element_atocc[89] = [ 2.0,  2.0,  0.0,       6.0,  0.0,  0.0,       1.0,  0.0,      0.0,  0.0,     0.0,      0.0,      0.0 ]
element_ptnl[89]  = [ 6.8,  7.8,  8.8,       6.8,  7.8,  8.8,       6.8,  7.8,      5.8,  6.8,     5.8,      6.8,      7.8 ]
element_idmd[89]  = [  0,    0,    1,         0,    0,    1,         0,    1,        0,    1,       0,        0,        0  ]


element_lmb[90]   = 6
element_ntle[90]  = [3,3,2,2,1,1,1]
element_augm[90]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[90] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     2.0,  0.0,     0.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[90]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[90]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[91]   = 6                                                                                                   
element_ntle[91]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[91]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[91] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     2.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[91]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[91]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[92]   = 6                                                                                                   
element_ntle[92]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[92]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[92] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     3.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[92]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[92]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[93]   = 6                                                                                                   
element_ntle[93]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[93]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[93] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     4.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[93]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[93]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[94]   = 6                                                                                                   
element_ntle[94]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[94]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[94] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     5.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[94]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[94]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[95]   = 6                                                                                                   
element_ntle[95]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[95]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[95] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     6.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[95]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[95]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[96]   = 6                                                                                                   
element_ntle[96]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[96]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[96] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     7.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[96]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[96]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[97]   = 6                                                                                                   
element_ntle[97]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[97]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[97] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     8.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[97]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[97]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[98]   = 6                                                                                                   
element_ntle[98]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[98]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[98] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,     9.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[98]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[98]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]

element_lmb[99]   = 6                                                                                                   
element_ntle[99]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[99]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[99] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    10.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[99]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[99]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]

element_lmb[100]   = 6                                                                                                   
element_ntle[100]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[100]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[100] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    11.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[100]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[100]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]

element_lmb[101]   = 6                                                                                                   
element_ntle[101]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[101]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[101] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    12.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[101]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[101]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]

element_lmb[102]   = 6                                                                                                   
element_ntle[102]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[102]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[102] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    13.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[102]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[102]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]


element_lmb[103]   = 6                                                                                                   
element_ntle[103]  = [3,3,2,2,1,1,1]                                                                                     
element_augm[103]  = ["LOC","APW","LOC",   "LOC","APW","LOC",   "APW","LOC",   "APW","LOC",    "APW",    "APW",    "APW"]
element_atocc[103] = [ 2.0,  2.0,  0.0,     6.0,  0.0,  0.0,     1.0,  0.0,    14.0,  0.0,      0.0,      0.0,      0.0 ]
element_ptnl[103]  = [ 6.8,  7.8,  8.8,     6.8,  7.8,  8.8,     6.8,  7.8,     5.8,  6.8,      5.8,      6.8,      7.8 ]
element_idmd[103]  = [  0,    0,    1,       0,    0,    1,       0,    1,       0,    1,        0,        0,        0  ]
