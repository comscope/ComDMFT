"""
Functions providing required symmetry related information about the
material structure.
"""
#
# SymGen is a dictionary of generators of a spacegroup.
#
# The key into the dictionary is the Hall spacegroup symbol of which there
# 530 distinct ones.
# Each entry contains a string which describes the generators of the space
# group. These strings are parsed by parsgn.f in the GW code. 
#
# The elements on the string are separated by underscores "_". Some elements
# involve vectors, denoted by "v" below. Vectors may be represented in two 
# different ways:
# 1. A vector may be given by its direction represented by floating point 
#    numbers. Examples: (1.0,1.0,1.0), (2.0,2.0,0.0), etc.
# 2. Common vectors may be given by a single character:
#    X  = (1.0,0.0,0.0)
#    Y  = (0.0,1.0,0.0)
#    Z  = (0.0,0.0,1.0)
#    D  = (1.0,1.0,1.0)
#    Sn = (0.0,0.0,1.0)/n (where n implies a rotation over 360/n degrees)
# Vectors are parsed by parsvc.f in the GW code.
#
# The symmetry operators are given by:
# 1. Rnv: n-fold rotation axis around v.
# 2. Mv : mirror plane perpendicular to v.
# 3. I  : inversion
# 4. E  : identity (usually not given as it is implied as the first
#         operation anyway)
# These operations may be modified to produce glide planes, screw axes, etc.
# The translation component can be represented by adding ":Tv". For example
# R4X:TZ represents a 4-fold rotation axis around X where every rotation is 
# followed by a translation of one unit in the Z-direction.
#
# If we are given a symmetry operation as a 3x3 matrix and we would like to
# know which one of the 4 operators above applies then consider the following
# facts:
#
# 1. The identity E leaves all vectors the same
# 2. A rotation R leaves one vector the same (the axis of rotation)
# 3. A mirror plane M leaves two vectors that span the plane the same,
#    and changes the sign of the vector normal to the mirror plane.
# 4. An inversion center I changes the sign of all vectors.
#
# Candidate vectors to use for testing:
#
# 1.  X    = (1.0,0.0,0.0)
# 2.  Y    = (0.0,1.0,0.0)
# 3.  Z    = (0.0,0.0,1.0)
# 4.  XY   = (1.0,1.0,0.0)
# 5.  XZ   = (1.0,0.0,1.0)
# 6.  YZ   = (0.0,1.0,1.0)
# 7.  -XY  = (-1.0,1.0,0.0)
# 8.  -XZ  = (-1.0,0.0,1.0)
# 9.  -YZ  = (0.0,-1.0,1.0)
# 10. D    = (1.0,1.0,1.0)
# 11. -XYZ = (-1.0,1.0,1.0)
# 12. X-YZ = (1.0,-1.0,1.0)
# 13. XY-Z = (1.0,1.0,-1.0)
#
# Given the vectors above the procedure to determine the symmetry operator
# corresponding to a given 3x3 matrix is:
#
#   1. Multiply the 3x3 matrix with all vectors
#   2. Check which vector remained unchanged
#   3. Count the number of unchanged vectors
#   4. If all vector unchanged -> E
#      else if one vector unchanged -> R
#        V = unchanged vector
#        n = 360/(angle between RnV*v and v), n should be 2, 3, 4, or 6.
#      else if two vectors unchanged and one vector changes sign -> M
#        if vector v1 changed sign then
#          v1 is the normal to the mirror plane
#      else if no vectors unchanged
#        if all vectors changed sign -> I
#        else -> Unidentified (U)
#      endif
#   5. Characterize the operator by the type of operator (E,R,M,I,U),
#      the vector V for R or M, and the order n.
#   6. Identify all primitive operations (ones that cannot be written as 
#      a product of others).
#   7. Build a string of symmetry generators by concatenating all generators
#      except the identity.
#
# References:
#
#   S.C. Miller, W.F. Love, "Tables of irreducible representations of space
#   groups and co-representations of magnetic space groups", Pruett Press,
#   Boulder, Co, 1967. Library of Congress #67-30015.
#
#   J.K. Dewhurst, S. Sharma, L. Nordstrom, "The Spacegroup Manual:
#   Version 1.2.1", https://elk.sourceforge.net/spacegroup.pdf
#
symgen = {}
                                                                                                        # spcgrp#, symm ops
# spacegroups 1-10
symgen["P 1"]     = ""                                                                                   # 1, C1
symgen["-P 1"]    = "I_"                                                                                 # 2, Inversion
symgen["P 2"]     = "R2Z_"                                                                               # 3, C2
symgen["P 2x"]    = "R2X_"                                                                               # 3, C2x
symgen["P 2y"]    = "R2Y_"                                                                               # 3, C2y
symgen["P 2c"]    = "R2Z:T(0.0,0.0,0.5)_"                                                                # 4, C2(0,0,1/2)
symgen["P 2xa"]   = "R2X:T(0.5,0.0,0.0)_"                                                                # 4, C2x(1/2,0,0)
symgen["P 2yb"]   = "R2Y:T(0.0,0.5,0.0)_"                                                                # 4, C2y(0,1/2,0)
symgen["A 2"]     = "R2Z_"                                                                               # 5, C2
symgen["B 2"]     = "R2Z_"                                                                               # 5, C2
symgen["I 2"]     = "R2Z_"                                                                               # 5, C2
symgen["B 2x"]    = "R2X_"                                                                               # 5, C2
symgen["C 2x"]    = "R2X_"                                                                               # 5, C2
symgen["I 2x"]    = "R2X_"                                                                               # 5, C2
symgen["A 2y"]    = "R2Y_"                                                                               # 5, C2
symgen["C 2y"]    = "R2Y_"                                                                               # 5, C2
symgen["I 2y"]    = "R2Y_"                                                                               # 5, C2
symgen["P -2"]    = "MZ_"                                                                                # 6, SigH
symgen["P -2x"]   = "MX_"                                                                                # 6, SigH
symgen["P -2y"]   = "MY_"                                                                                # 6, SigH
symgen["P -2a"]   = "MZ:T(0.5,0.0,0.0)_"                                                                 # 7, SigH(1/2,0,0)
symgen["P -2ab"]  = "MZ:T(0.5,0.5,0.0)_"                                                                 # 7, SigH(1/2,1/2,0)
symgen["P -2b"]   = "MZ:T(0.0,0.5,0.0)_"                                                                 # 7, SigH(0,1/2,0)
symgen["P -2xb"]  = "MX:T(0.0,0.5,0.0)_"                                                                 # 7, SigH(0,1/2,0)
symgen["P -2xbc"] = "MX:T(0.0,0.5,0.5)_"                                                                 # 7, SigH(0,1/2,1/2)
symgen["P -2xc"]  = "MX:T(0.0,0.0,0.5)_"                                                                 # 7, SigH(0,0,1/2)
symgen["P -2ya"]  = "MY:T(0.5,0.0,0.0)_"                                                                 # 7, SigH(1/2,0,0)
symgen["P -2yac"] = "MY:T(0.5,0.0,0.5)_"                                                                 # 7, SigH(1/2,0,1/2)
symgen["P -2yc"]  = "MY:T(0.0,0.0,0.5)_"                                                                 # 7, SigH(0,0,1/2)
symgen["A -2"]    = "MZ_"                                                                                # 8, SigH
symgen["B -2"]    = "MZ_"                                                                                # 8, SigH
symgen["I -2"]    = "MZ_"                                                                                # 8, SigH
symgen["B -2x"]   = "MX_"                                                                                # 8, SigH
symgen["C -2x"]   = "MX_"                                                                                # 8, SigH
symgen["I -2x"]   = "MX_"                                                                                # 8, SigH
symgen["A -2y"]   = "MY_"                                                                                # 8, SigH
symgen["C -2y"]   = "MY_"                                                                                # 8, SigH
symgen["I -2y"]   = "MY_"                                                                                # 8, SigH
symgen["A -2a"]   = "MZ:T(0.5,0.0,0.0)_"                                                                 # 9, SigH(1/2,0,0)
symgen["A -2ac"]  = "MZ:T(0.5,0.0,0.5)_"                                                                 # 9, SigH(1/2,0,1/2)
symgen["B -2b"]   = "MZ:T(0.0,0.5,0.0)_"                                                                 # 9, SigH(0,1/2,0)
symgen["B -2bc"]  = "MZ:T(0.0,0.5,0.5)_"                                                                 # 9, SigH(0,1/2,1/2)
symgen["I -2a"]   = "MZ:T(0.5,0.0,0.0)_"                                                                 # 9, SigH(1/2,0,0)
symgen["I -2b"]   = "MZ:T(0.0,0.5,0.0)_"                                                                 # 9, SigH(0,1/2,0)
symgen["B -2xb"]  = "MX:T(0.0,0.5,0.0)_"                                                                 # 9, SigH(0,1/2,0)
symgen["B -2xbc"] = "MX:T(0.0,0.5,0.5)_"                                                                 # 9, SigH(0,1/2,1/2)
symgen["C -2xc"]  = "MX:T(0.0,0.0,0.5)_"                                                                 # 9, SigH(0,0,1/2)
symgen["C -2xbc"] = "MX:T(0.0,0.5,0.5)_"                                                                 # 9, SigH(0,1/2,1/2)
symgen["I -2xb"]  = "MX:T(0.0,0.5,0.0)_"                                                                 # 9, SigH(0,1/2,0)
symgen["I -2xc"]  = "MX:T(0.0,0.0,0.5)_"                                                                 # 9, SigH(0,0,1/2)
symgen["-P 2"]    = "R2Z_MZ_"                                                                            # 10, C2,SigH
symgen["-P 2x"]   = "R2X_MX_"                                                                            # 10, C2,SigH
symgen["-P 2y"]   = "R2Y_MY_"                                                                            # 10, C2,SigH
# spacegroups 11-20
symgen["-P 2c"]   = "R2Z:T(0.0,0.0,0.5)_MZ:T(0.0,0.0,0.5)_"                                              # 11, C2(0,0,1/2),SigH(0,0,1/2)
symgen["-P 2xa"]  = "R2X:T(0.5,0.0,0.0)_MX:T(0.5,0.0,0.0)_"                                              # 11, C2(1/2,0,0),SigH(1/2,0,0)
symgen["-P 2yb"]  = "R2Y:T(0.0,0.5,0.0)_MY:T(0.0,0.5,0.0)_"                                              # 11, C2(0,1/2,0),SigH(0,1/2,0)
symgen["-A 2"]    = "R2Z_MZ_"                                                                            # 12, C2,SigH
symgen["-B 2"]    = "R2Z_MZ_"                                                                            # 12, C2,SigH
symgen["-I 2"]    = "R2Z_MZ_"                                                                            # 12, C2,SigH
symgen["-B 2x"]   = "R2X_MX_"                                                                            # 12, C2,SigH
symgen["-C 2x"]   = "R2X_MX_"                                                                            # 12, C2,SigH
symgen["-I 2x"]   = "R2X_MX_"                                                                            # 12, C2,SigH
symgen["-A 2y"]   = "R2Y_MY_"                                                                            # 12, C2,SigH
symgen["-C 2y"]   = "R2Y_MY_"                                                                            # 12, C2,SigH
symgen["-I 2y"]   = "R2Y_MY_"                                                                            # 12, C2,SigH
symgen["-P 2a"]   = "R2Z:T(0.5,0.0,0.0)_MZ:T(0.5,0.0,0.0)_"                                              # 13, C2(1/2,0,0),SigH(1/2,0,0)
symgen["-P 2ab"]  = "R2Z:T(0.5,0.5,0.0)_MZ:T(0.5,0.5,0.0)_"                                              # 13, C2(1/2,1/2,0),SigH(1/2,1/2,0)
symgen["-P 2b"]   = "R2Z:T(0.0,0.5,0.0)_MZ:T(0.0,0.5,0.0)_"                                              # 13, C2(0,1/2,0),SigH(0,1/2,0)
symgen["-P 2xb"]  = "R2X:T(0.0,0.5,0.0)_MX:T(0.0,0.5,0.0)_"                                              # 13, C2(0,1/2,0),SigH(0,1/2,0)
symgen["-P 2xbc"] = "R2X:T(0.0,0.5,0.5)_MX:T(0.0,0.5,0.5)_"                                              # 13, C2(0,1/2,1/2),SigH(0,1/2,1/2)
symgen["-P 2xc"]  = "R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                                              # 13, C2(0,0,1/2),SigH(0,0,1/2)
symgen["-P 2ya"]  = "R2Y:T(0.5,0.0,0.0)_MY:T(0.5,0.0,0.0)_"                                              # 13, C2(1/2,0,0),SigH(1/2,0,0)
symgen["-P 2yac"] = "R2Y:T(0.5,0.0,0.5)_MY:T(0.5,0.0,0.5)_"                                              # 13, C2(1/2,0,1/2),SigH(1/2,0,1/2)
symgen["-P 2yc"]  = "R2Y:T(0.0,0.0,0.5)_MY:T(0.0,0.0,0.5)_"                                              # 13, C2(0,0,1/2),SigH(0,0,1/2)
symgen["-P 2ac"]  = "R2Z:T(0.5,0.0,0.5)_MZ:T(0.5,0.0,0.5)_"                                              # 14, C2(1/2,0,1/2),SigH(1/2,0,1/2)
symgen["-P 2bc"]  = "R2Z:T(0.0,0.5,0.5)_MZ:T(0.0,0.5,0.5)_"                                              # 14, C2(0,1/2,1/2),SigH(0,1/2,1/2)
symgen["-P 2n"]   = "R2Z:T(0.5,0.5,0.5)_MZ:T(0.5,0.5,0.5)_"                                              # 14, C2(1/2,1/2,1/2),SigH(1/2,1/2,1/2)
symgen["-P 2xab"] = "R2X:T(0.5,0.5,0.0)_MX:T(0.5,0.5,0.0)_"                                              # 14, C2(1/2,1/2,0),SigH(1/2,1/2,0)
symgen["-P 2xac"] = "R2X:T(0.5,0.0,0.5)_MX:T(0.5,0.0,0.5)_"                                              # 14, C2(1/2,0,1/2),SigH(1/2,0,1/2)
symgen["-P 2xn"]  = "R2X:T(0.5,0.5,0.5)_MX:T(0.5,0.5,0.5)_"                                              # 14, C2(1/2,1/2,1/2),SigH(1/2,1/2,1/2)
symgen["-P 2yab"] = "R2Y:T(0.5,0.5,0.0)_MY:T(0.5,0.5,0.0)_"                                              # 14, C2(1/2,1/2,0),SigH(1/2,1/2,0)
symgen["-P 2ybc"] = "R2Y:T(0.0,0.5,0.5)_MY:T(0.0,0.5,0.5)_"                                              # 14, C2(0,1/2,1/2),SigH(0,1/2,1/2)
symgen["-P 2yn"]  = "R2Y:T(0.5,0.5,0.5)_MY:T(0.5,0.5,0.5)_"                                              # 14, C2(1/2,1/2,1/2),SigH(1/2,1/2,1/2)
symgen["-A 2a"]   = "R2Z:T(0.5,0.0,0.0)_MZ:T(0.5,0.0,0.0)_"                                              # 15, C2(1/2,0,0),SigH(1/2,0,0)
symgen["-A 2ac"]  = "R2Z:T(0.5,0.0,0.5)_MZ:T(0.5,0.0,0.5)_"                                              # 15, C2(1/2,0,1/2),SigH(1/2,0,1/2)
symgen["-B 2b"]   = "R2Z:T(0.0,0.5,0.0)_MZ:T(0.0,0.5,0.0)_"                                              # 15, C2(0,1/2,0),SigH(0,1/2,0)
symgen["-B 2bc"]  = "R2Z:T(0.0,0.5,0.5)_MZ:T(0.0,0.5,0.5)_"                                              # 15, C2(0,1/2,1/2),SigH(0,1/2,1/2)
symgen["-I 2a"]   = "R2Z:T(0.5,0.0,0.0)_MZ:T(0.5,0.0,0.0)_"                                              # 15, C2(1/2,0,0),SigH(1/2,0,0)
symgen["-I 2b"]   = "R2Z:T(0.0,0.5,0.0)_MZ:T(0.0,0.5,0.0)_"                                              # 15, C2(0,1/2,0),SigH(0,1/2,0)
symgen["-B 2xb"]  = "R2X:T(0.0,0.5,0.0)_MX:T(0.0,0.5,0.0)_"                                              # 15, C2(0,1/2,0),SigH(0,1/2,0)
symgen["-B 2xbc"] = "R2X:T(0.0,0.5,0.5)_MX:T(0.0,0.5,0.5)_"                                              # 15, C2(0,1/2,1/2),SigH(0,1/2,1/2)
symgen["-C 2xbc"] = "R2X:T(0.0,0.5,0.5)_MX:T(0.0,0.5,0.5)_"                                              # 15, C2(0,1/2,1/2),SigH(0,1/2,1/2)
symgen["-C 2xc"]  = "R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                                              # 15, C2(0,0,1/2),SigH(0,0,1/2)
symgen["-I 2xb"]  = "R2X:T(0.0,0.5,0.0)_MX:T(0.0,0.5,0.0)_"                                              # 15, C2(0,1/2,0),SigH(0,1/2,0)
symgen["-I 2xc"]  = "R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                                              # 15, C2(0,0,1/2),SigH(0,0,1/2)
symgen["-A 2ya"]  = "R2Y:T(0.5,0.0,0.0)_MY:T(0.5,0.0,0.0)_"                                              # 15, C2(1/2,0,0),SigH(1/2,0,0)
symgen["-A 2yac"] = "R2Y:T(0.5,0.0,0.5)_MY:T(0.5,0.0,0.5)_"                                              # 15, C2(1/2,0,1/2),SigH(1/2,0,1/2)
symgen["-C 2ybc"] = "R2Y:T(0.0,0.5,0.5)_MY:T(0.0,0.5,0.5)_"                                              # 15, C2(0,1/2,1/2),SigH(0,1/2,1/2)
symgen["-C 2yc"]  = "R2Y:T(0.0,0.0,0.5)_MY:T(0.0,0.0,0.5)_"                                              # 15, C2(0,0,1/2),SigH(0,0,1/2)
symgen["-I 2ya"]  = "R2Y:T(0.5,0.0,0.0)_MY:T(0.5,0.0,0.0)_"                                              # 15, C2(1/2,0,0),SigH(1/2,0,0)
symgen["-I 2yc"]  = "R2Y:T(0.0,0.0,0.5)_MY:T(0.0,0.0,0.5)_"                                              # 15, C2(0,0,1/2),SigH(0,0,1/2)
symgen["P 2 2"]   = "R2Z_R2X_"                                                                           # 16, C2,C2x
symgen[17]  = "R2Z:T(0.0,0.0,0.5)_R2X_"                                                            # C2(0,0,1/2),C2x
symgen[18]  = "R2Z_R2X:T(0.5,0.5,0.0)_"                                                            # C2,C2x(1/2,1/2,0)
symgen[19]  = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_"                                             # C2(1/2,0,1/2),C2x(1/2,1/2,0)
symgen[20]  = "R2Z:T(0.0,0.0,0.5)_R2X_"                                                            # C2(0,0,1/2),C2x
symgen[21]  = "R2Z_R2X_"                                                                           # C2,C2x
symgen[22]  = "R2Z_R2X_"                                                                           # C2,C2x
symgen[23]  = "R2Z_R2X_"                                                                           # C2,C2x
symgen[24]  = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_"                                             # C2(1/2,0,1/2),C2x(1/2,1/2,0)
symgen[25]  = "R2Z_MX_"                                                                            # C2,SigV
symgen[26]  = "R2Z:T(0.0,0.0,0.5)_MX_"                                                             # C2(0,0,1/2),SigV
symgen[27]  = "R2Z_MX:T(0.0,0.0,0.5)_"                                                             # C2,SigV(0,0,1/2)
symgen[28]  = "R2Z_MX:T(0.5,0.0,0.0)_"                                                             # C2,SigV(1/2,0,0)
symgen[29]  = "R2Z:T(0.0,0.0,0.5)_MX:T(0.5,0.0,0.5)_"                                              # C2(0,0,1/2),SigV(1/2,0,1/2)
symgen[30]  = "R2Z_MX:T(0.0,0.5,0.5)_"                                                             # C2,SigV(0,1/2,1/2)
symgen[31]  = "R2Z:T(0.5,0.0,0.5)_MX_"                                                             # C2(1/2,0,1/2),SigV
symgen[32]  = "R2Z_MX:T(0.5,0.5,0.0)_"                                                             # C2,SigV(1/2,1/2,0)
symgen[33]  = "R2Z:T(0.0,0.0,0.5)_MX:T(0.5,0.5,0.5)_"                                              # C2(0,0,1/2),SigV(1/2,1/2,1/2)
symgen[34]  = "R2Z_MX:T(0.5,0.5,0.5)_"                                                             # C2,SigV(1/2,1/2,1/2)
symgen[35]  = "R2Z_MX_"                                                                            # C2,SigV
symgen[36]  = "R2Z:T(0.0,0.0,0.5)_MX_"                                                             # C2(0,0,1/2),SigV
symgen[37]  = "R2Z_MX:T(0.0,0.0,0.5)_"                                                             # C2,SigV(0,0,1/2)
symgen[38]  = "R2X_MZ_"                                                                            # C2X,SigH
symgen[39]  = "R2X_MZ:T(0.0,0.5,0.0)_"                                                             # C2X,SigH(0,1/2,0)
symgen[40]  = "R2X_MZ:T(0.0,0.0,0.5)_"                                                             # C2X,SigH(0,0,1/2)
symgen[41]  = "R2X_MZ:T(0.0,0.5,0.5)_"                                                             # C2X,SigH(0,1/2,1/2)
symgen[42]  = "R2Z_MX_"                                                                            # C2,SigV
symgen[43]  = "R2Z_MX:T(0.25,0.25,0.25)_"                                                          # C2,SigV(1/4,1/4,1/4)
symgen[44]  = "R2Z_MX_"                                                                            # C2,SigV
symgen[45]  = "R2Z_MX:T(0.0,0.0,0.5)_"                                                             # C2,SigV(0,0,1/2)
symgen[46]  = "R2Z_MX:T(0.5,0.0,0.0)_"                                                             # C2,SigV(1/2,0,0)
symgen[47]  = "R2Z_R2X_MX_"                                                                        # C2,C2X,SigV
symgen[48]  = "R2Z_R2X_MX:T(0.5,0.5,0.5)_"                                                         # C2,C2X,SigV(1/2,1/2,1/2)
symgen[49]  = "R2Z_R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                                          # C2,C2X(0,0,1/2),SigV(0,0,1/2)
symgen[50]  = "R2Z_R2X_MX:T(0.5,0.5,0.0)_"                                                         # C2,C2X,SigV(1/2,1/2,0)
symgen[51]  = "R2Z:T(0.5,0.0,0.0)_R2X:T(0.5,0.0,0.0)_MX:T(0.5,0.0,0.0)_"                           # C2(1/2,0,0),C2X(1/2,0,0),SigV(1/2,0,0)
symgen[52]  = "R2Z:T(0.5,0.0,0.0)_R2X:T(0.0,0.5,0.5)_MX:T(0.0,0.5,0.5)_"                           # C2(1/2,0,0),C2X(0,1/2,1/2),SigV(0,1/2,1/2)
symgen[53]  = "R2Z:T(0.5,0.0,0.5)_R2X_MX_"                                                         # C2(1/2,0,1/2),C2X,SigV
symgen[54]  = "R2Z:T(0.5,0.0,0.0)_R2X:T(0.5,0.0,0.5)_MX:T(0.5,0.0,0.5)_"                           # C2(1/2,0,0),C2X(1/2,0,1/2),SigV(1/2,0,1/2)
symgen[55]  = "R2Z_R2X:T(0.5,0.5,0.0)_MX:T(0.5,0.5,0.0)_"                                          # C2,C2X(1/2,1/2,0),SigV(1/2,1/2,0)
symgen[56]  = "R2Z:T(0.5,0.5,0.0)_R2X:T(0.5,0.0,0.5)_MX:T(0.5,0.0,0.5)_"                           # C2(1/2,1/2,0),C2X(1/2,0,1/2),SigV(1/2,0,1/2)
symgen[57]  = "R2Z:T(0.0,0.0,0.5)_R2X:T(0.0,0.5,0.0)_MX:T(0.0,0.5,0.0)_"                           # C2(0,0,1/2),C2X(0,1/2,0),SigV(0,1/2,0)
symgen[58]  = "R2Z_R2X:T(0.5,0.5,0.5)_MX:T(0.5,0.5,0.5)_"                                          # C2,C2X(1/2,1/2,1/2),SigV(1/2,1/2,1/2)
symgen[59]  = "R2Z_R2X:T(0.5,0.5,0.0)_MX_"                                                         # C2,C2X(1/2,1/2,0),SigV
symgen[60]  = "R2Z:T(0.5,0.5,0.5)_R2X:T(0.5,0.5,0.0)_MX:T(0.5,0.5,0.0)_"                           # C2(1/2,1/2,1/2),C2X(1/2,1/2,0),SigV(1/2,1/2,0)
symgen[61]  = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_MX:T(0.5,0.5,0.0)_"                           # C2(1/2,0,1/2),C2X(1/2,1/2,0),SigV(1/2,1/2,0)
symgen[62]  = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.5)_MX:T(0.5,0.5,0.5)_"                           # C2(1/2,0,1/2),C2X(1/2,1/2,1/2),SigV(1/2,1/2,1/2)
symgen[63]  = "R2Z:T(0.0,0.0,0.5)_R2X_MX_"                                                         # C2(0,0,1/2),C2X,SigV
symgen[64]  = "R2Z:T(0.0,0.5,0.5)_R2X_MX_"                                                         # C2(0,1/2,1/2),C2X,SigV
symgen[65]  = "R2Z_R2X_MX_"                                                                        # C2,C2X,SigV
symgen[66]  = "R2Z_R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                                          # C2,C2X(0,0,1/2),SigV(0,0,1/2)
symgen[67]  = "R2Z:T(0.5,0.0,0.0)_R2X_MX_"                                                         # C2(1/2,0,0),C2X,SigV
symgen[68]  = "R2Z_R2X_MX:T(0.0,0.5,0.5)_"                                                         # C2,C2X,SigV(0,1/2,1/2)
symgen[69]  = "R2Z_R2X_MX_"                                                                        # C2,C2X,SigV
symgen[70]  = "R2Z_R2X_MX:T(0.25,0.25,0.25)_"                                                      # C2,C2X,SigV(1/4,1/4,1/4)
symgen[71]  = "R2Z_R2X_MX_"                                                                        # C2,C2X,SigV
symgen[72]  = "R2Z_R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                                          # C2,C2X(0,0,1/2),SigV(0,0,1/2)
symgen[73]  = "R2Z:T(0.0,0.5,0.0)_R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                           # C2(0,1/2,0),C2X(0,0,1/2),SigV(0,0,1/2)
symgen[74]  = "R2Z:T(0.0,0.5,0.0)_R2X_MX_"                                                         # C2(0,1/2,0),C2X,SigV
symgen[75]  = "R2Z_R4Z_"                                                                           # C2,C4
symgen[76]  = "R2Z:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.25)_"                                            # C2(0,0,1/2),C4(0,0,1/4)
symgen[77]  = "R2Z_R4Z:T(0.0,0.0,0.5)_"                                                            # C2,C4(0,0,1/4)
symgen[78]  = "R2Z:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.75)_"                                            # C2(0,0,1/2),C4(0,0,3/4)
symgen[79]  = "R2Z_R4Z_"                                                                           # C2,C4
symgen[80]  = "R2Z_R4Z:T(0.0,0.5,0.25)_"                                                           # C2,C4(0,1/2,1/4)
symgen[81]  = "R2Z_R4Z:T(0.0,0.0,0.25)_"                                                           # C2,S4
symgen[82]  = "R2Z_R4Z:T(0.0,0.0,0.25)_"                                                           # C2,S4
symgen[83]  = "R2Z_MZ_R4Z_"                                                                        # C2,SigH,C4
symgen[84]  = "R2Z_MZ_R4Z:T(0.0,0.0,0.5)_"                                                         # C2,SigH,C4(0,0,1/2)
symgen[85]  = "R2Z_MZ:T(0.5,0.5,0.0)_R4Z:T(0.5,0.5,0.0)_"                                          # C2,SigH(1/2,1/2,0),C4(1/2,1/2,0)
symgen[86]  = "R2Z_MZ:T(0.5,0.5,0.5)_R4Z:T(0.5,0.5,0.5)_"                                          # C2,SigH(1/2,1/2,1/2),C4(1/2,1/2,1/2)
symgen[87]  = "R2Z_MZ_R4Z_"                                                                        # C2,SigH,C4
symgen[88]  = "R2Z_MZ:T(0.0,0.5,0.25)_R4Z:T(0.0,0.5,0.25)_"                                        # C2,SigH(0,1/2,1/4),C4(0,1/2,1/4)
symgen[89]  = "R2Z_R2X_R4Z_"                                                                       # C2,C2X,C4
symgen[90]  = "R2Z_R2X:T(0.5,0.5,0.0)_R4Z:T(0.5,0.5,0.0)_"                                         # C2,C2X(1/2,1/2,0),C4(1/2,1/2,0)
symgen[91]  = "R2Z:T(0.0,0.0,0.5)_R2X:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.25)_"                         # C2(0,0,1/2),C2X(0,0,1/2),C4(0,0,1/4)
symgen[92]  = "R2Z:T(0.0,0.0,0.5)_R2X:T(0.5,0.5,0.75)_R4Z:T(0.5,0.5,0.25)_"                        # C2(0,0,1/2),C2X(1/2,1/2,3/4),C4(1/2,1/2,1/4)
symgen[93]  = "R2Z_R2X_R4Z:T(0.0,0.0,0.5)_"                                                        # C2,C2X,C4(0,0,1/2)
symgen[94]  = "R2Z_R2X:T(0.5,0.5,0.5)_R4Z:T(0.5,0.5,0.5)_"                                         # C2,C2X(1/2,1/2,1/2),C4(1/2,1/2,1/2)
symgen[95]  = "R2Z:T(0.0,0.0,0.5)_R2X:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.75)_"                         # C2(0,0,1/2),C2X(0,0,1/2),C4(0,0,3/4)
symgen[96]  = "R2Z:T(0.0,0.0,0.5)_R2X:T(0.5,0.5,0.25)_R4Z:T(0.5,0.5,0.75)_"                        # C2(0,0,1/2),C2X(1/2,1/2,1/4),C4(1/2,1/2,3/4)
symgen[97]  = "R2Z_R2X_R4Z_"                                                                       # C2,C2X,C4
symgen[98]  = "R2Z_R2X:T(0.0,0.5,0.25)_R4Z:T(0.0,0.5,0.25)_"                                       # C2,C2X(0,1/2,1/4),C4(0,1/2,1/4)
symgen[99]  = "R2Z_MX_R4Z_"                                                                        # C2,SigV,C4
symgen[100] = "R2Z_MX:T(0.5,0.5,0.0)_R4Z_"                                                         # C2,SigV(1/2,1/2,0),C4
symgen[101] = "R2Z_MX:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.5)_"                                          # C2,SigV(0,0,1/2),C4(0,0,1/2)
symgen[102] = "R2Z_MX:T(0.5,0.5,0.5)_R4Z:T(0.5,0.5,0.5)_"                                          # C2,SigV(1/2,1/2,1/2),C4(1/2,1/2,1/2)
symgen[103] = "R2Z_MX:T(0.0,0.0,0.5)_R4Z_"                                                         # C2,SigV(0,0,1/2),C4
symgen[104] = "R2Z_MX:T(0.5,0.5,0.5)_R4Z_"                                                         # C2,SigV(1/2,1/2,1/2),C4
symgen[105] = "R2Z_MX_R4Z:T(0.0,0.0,0.5)_"                                                         # C2,SigV,C4(0,0,1/2)
symgen[106] = "R2Z_MX:T(0.5,0.5,0.0)_R4Z:T(0.0,0.0,0.5)_"                                          # C2,SigV(1/2,1/2,0),C4(0,0,1/2)
symgen[107] = "R2Z_MX_R4Z_"                                                                        # C2,SigV,C4
symgen[108] = "R2Z_MX:T(0.5,0.5,0.0)_R4Z_"                                                         # C2,SigV(1/2,1/2,0),C4
symgen[109] = "R2Z_MX_R4Z:T(0.0,0.5,0.25)_"                                                        # C2,SigV,C4(0,1/2,1/4)
symgen[110] = "R2Z_MX:T(0.0,0.0,0.5)_R4Z:T(0.0,0.5,0.25)_"                                         # C2,SigV(0,0,1/2),C4(0,1/2,1/4)
symgen[111] = "R2Z_R2X_R4Z:T(0.0,0.0,0.25)_"                                                       # C2,C2X,S4
symgen[112] = "R2Z_R2X:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.25)_"                                        # C2,C2X(0,0,1/2),S4
symgen[113] = "R2Z_R2X:T(0.5,0.5,0.0)_R4Z:T(0.0,0.0,0.25)_"                                        # C2,C2X(1/2,1/2,0),S4
symgen[114] = "R2Z_R2X:T(0.5,0.5,0.5)_R4Z:T(0.0,0.0,0.25)_"                                        # C2,C2X(1/2,1/2,1/2),S4
symgen[115] = "R2Z_MX_R4Z:T(0.0,0.0,0.25)_"                                                        # C2,SigV,S4
symgen[116] = "R2Z_MX:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.25)_"                                         # C2,SigV(0,0,1/2),S4
symgen[117] = "R2Z_MX:T(0.5,0.5,0.0)_R4Z:T(0.0,0.0,0.25)_"                                         # C2,SigV(1/2,1/2,0),S4
symgen[118] = "R2Z_MX:T(0.5,0.5,0.5)_R4Z:T(0.0,0.0,0.25)_"                                         # C2,SigV(1/2,1/2,1/2),S4
symgen[119] = "R2Z_MX_R4Z:T(0.0,0.0,0.25)_"                                                        # C2,SigV,S4
symgen[120] = "R2Z_MX:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.25)_"                                         # C2,SigV(0,0,1/2),S4
symgen[121] = "R2Z_R2X_R4Z:T(0.0,0.0,0.25)_"                                                       # C2,C2X,S4
symgen[122] = "R2Z_R2X:T(0.0,0.5,0.25)_R4Z:T(0.0,0.0,0.25)_"                                       # C2,C2X(0,1/2,1/4),S4
symgen[123] = "R2Z_R2X_MX_R4Z_"                                                                    # C2,C2X,SigV,C4
symgen[124] = "R2Z_R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_R4Z_"                                      # C2,C2X(0,0,1/2),SigV(0,0,1/2),C4
symgen[125] = "R2Z_R2X_MX:T(0.5,0.5,0.0)_R4Z_"                                                     # C2,C2X,SigV(1/2,1/2,0),C4
symgen[126] = "R2Z_R2X_MX:T(0.5,0.5,0.5)_R4Z_"                                                     # C2,C2X,SigV(1/2,1/2,1/2),C4
symgen[127] = "R2Z_R2X:T(0.5,0.5,0.0)_MX:T(0.5,0.5,0.0)_R4Z_"                                      # C2,C2X(1/2,1/2,0),SigV(1/2,1/2,0),C4
symgen[128] = "R2Z_R2X:T(0.5,0.5,0.5)_MX:T(0.5,0.5,0.5)_R4Z_"                                      # C2,C2X(1/2,1/2,1/2),SigV(1/2,1/2,1/2),C4
symgen[129] = "R2Z_R2X:T(0.5,0.5,0.0)_MX_R4Z:T(0.5,0.5,0.0)_"                                      # C2,C2X(1/2,1/2,0),SigV,C4(1/2,1/2,0)
symgen[130] = "R2Z_R2X:T(0.5,0.5,0.5)_MX:T(0.0,0.0,0.5)_R4Z:T(0.5,0.5,0.0)_"                       # C2,C2X(1/2,1/2,1/2),SigV(0,0,1/2),C4(1/2,1/2,0)
symgen[131] = "R2Z_R2X_MX_R4Z:T(0.0,0.0,0.5)_"                                                     # C2,C2X,SigV,C4(0,0,1/2)
symgen[132] = "R2Z_R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_R4Z:T(0.0,0.0,0.5)_"                       # C2,C2X(0,0,1/2),SigV(0,0,1/2),C4(0,0,1/2)
symgen[133] = "R2Z_R2X:T(0.0,0.0,0.5)_MX:T(0.5,0.5,0.0)_R4Z:T(0.5,0.5,0.5)_"                       # C2,C2X(0,0,1/2),SigV(1/2,1/2,0),C4(1/2,1/2,1/2)
symgen[134] = "R2Z_R2X_MX:T(0.5,0.5,0.5)_R4Z:T(0.5,0.5,0.5)_"                                      # C2,C2X,SigV(1/2,1/2,1/2),C4(1/2,1/2,1/2)
symgen[135] = "R2Z_R2X:T(0.5,0.5,0.0)_MX:T(0.5,0.5,0.0)_R4Z:T(0.0,0.0,0.5)_"                       # C2,C2X(1/2,1/2,0),SigV(1/2,1/2,0),C4(0,0,1/2)
symgen[136] = "R2Z_R2X:T(0.5,0.5,0.5)_MX:T(0.5,0.5,0.5)_R4Z:T(0.5,0.5,0.5)_"                       # C2,C2X(1/2,1/2,1/2),SigV(1/2,1/2,1/2),C4(1/2,1/2,1/2)
symgen[137] = "R2Z_R2X:T(0.5,0.5,0.5)_MX_R4Z:T(0.5,0.5,0.5)_"                                      # C2,C2X(1/2,1/2,1/2),SigV,C4(1/2,1/2,1/2)
symgen[138] = "R2Z_R2X:T(0.5,0.5,0.0)_MX:T(0.0,0.0,0.5)_R4Z:T(0.5,0.5,0.5)_"                       # C2,C2X(1/2,1/2,0),SigV(0,0,1/2),C4(1/2,1/2,1/2)
symgen[139] = "R2Z_R2X_MX_R4Z_"                                                                    # C2,C2X,SigV,C4
symgen[140] = "R2Z_R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_R4Z_"                                      # C2,C2X(0,0,1/2),SigV(0,0,1/2),C4
symgen[141] = "R2Z_R2X:T(0.0,0.5,0.25)_MX_R4Z:T(0.0,0.5,0.25)_"                                    # C2,C2X(0,1/2,1/4),SigV,C4(0,1/2,1/4)
symgen[142] = "R2Z_R2X:T(0.0,0.5,0.75)_MX:T(0.0,0.0,0.5)_R4Z:T(0.0,0.5,0.25)_"                     # C2,C2X(0,1/2,3/4),SigV(0,0,1/2),C4(0,1/2,1/4)
symgen[143] = "R3Z_"                                                                               # C3
symgen[144] = "R3Z:T(0.0,0.0,0.3333333333333333)_"                                                 # C3(0,0,1/3)
symgen[145] = "R3Z:T(0.0,0.0,0.6666666666666667)_"                                                 # C3(0,0,2/3)
symgen[146] = "R3Z_"                                                                               # C3
symgen[147] = "R3Z_I_"                                                                             # C3,I
symgen[148] = "R3Z_I_"                                                                             # C3,I
symgen[149] = "R3Z_R2D_"                                                                           # C3,C2D
symgen[150] = "R3Z_R2X_"                                                                           # C3,C2X
symgen[151] = "R3Z:T(0.0,0.0,0.3333333333333333)_R2D_"                                             # C3(0,0,1/3),C2D
symgen[152] = "R3Z:T(0.0,0.0,0.3333333333333333)_R2X:T(0.0,0.0,0.6666666666666667)_"               # C3(0,0,1/3),C2X(0,0,2/3)
symgen[153] = "R3Z:T(0.0,0.0,0.6666666666666667)_R2D_"                                             # C3(0,0,2/3),C2D
symgen[154] = "R3Z:T(0.0,0.0,0.6666666666666667)_R2X_"                                             # C3(0,0,2/3),C2X
symgen[155] = "R3Z_R2X_"                                                                           # C3,C2X
symgen[156] = "R3Z_M(-1.0,1.0,0.0)_"                                                               # C3,SigD
symgen[157] = "R3Z_MX_"                                                                            # C3,SigV
symgen[158] = "R3Z_M(-1.0,1.0,0.0):T(0.0,0.0,0.5)_"                                                # C3,SigD(0,0,1/2)
symgen[159] = "R3Z_MX:T(0.0,0.0,0.5)_"                                                             # C3,SigV(0,0,1/2)     
symgen[160] = "R3Z_M(-1.0,1.0,0.0)_"                                                               # C3,SigD
symgen[161] = "R3Z_M(-1.0,1.0,0.0):T(0.0,0.0,0.5)_"                                                # C3,SigD(0,0,1/2)
symgen[162] = "R3Z_I_MX_"                                                                          # C3,I,SigV
symgen[163] = "R3Z_I_MX:T(0.0,0.0,0.5)_"                                                           # C3,I,SigV(0,0,1/2)
symgen[164] = "R3Z_I_M(-1.0,1.0,0.0)_"                                                             # C3,I,SigD
symgen[165] = "R3Z_I_M(-1.0,1.0,0.0):T(0.0,0.0,0.5)_"                                              # C3,I,SigD(0,0,1/2)
symgen[166] = "R3Z_I_M(-1.0,1.0,0.0)_"                                                             # C3,I,SigD
symgen[167] = "R3Z_I_M(-1.0,1.0,0.0):T(0.0,0.0,0.5)_"                                              # C3,I,SigD(0,0,1/2)
symgen["P 6"]              = "R3Z_R2Z_"                                                            # C3,C2 (actually C6)
symgen["P 61"]             = "R3Z:T(0.0,0.0,0.3333333333333333)_R2Z:T(0.0,0.0,0.5)_"               # C3(0,0,1/3),C2(0,0,1/2)
symgen["P 65"]             = "R3Z:T(0.0,0.0,0.6666666666666667)_R2Z:T(0.0,0.0,0.5)_"               # C3(0,0,2/3),C2(0,0,1/2)
symgen["P 62"]             = "R3Z:T(0.0,0.0,0.6666666666666667)_R2Z_"                              # C3(0,0,2/3),C2
symgen["P 64"]             = "R3Z:T(0.0,0.0,0.3333333333333333)_R2Z_"                              # C3(0,0,1/3),C2
symgen["P 6c"]             = "R3Z_R2Z:T(0.0,0.0,0.5)_"                                             # C3,C2(0,0,1/2)
symgen["P -6"]             = "R3Z_MZ_"                                                             # C3,SigH
symgen["-P 6"]             = "R3Z_I_MZ_"                                                           # C3,I,SigH
symgen["-P 6c"]            = "R3Z_I_MZ:T(0.0,0.0,0.5)_"                                            # C3,I,SigH(0,0,1/2)
symgen["P 6 2"]            = "R3Z_R2Z_R2X_"                                                        # C3,C2,C2X 
symgen["P 61 2 (0 0 -1)"]  = "R3Z:T(0.0,0.0,0.3333333333333333)_R2Z:T(0.0,0.0,0.5)_R2X_"           # C3(0,0,1/3),C2(0,0,1/2),C2X 
symgen["P 65 2 (0 0 1)"]   = "R3Z:T(0.0,0.0,0.6666666666666667)_R2Z:T(0.0,0.0,0.5)_R2X_"           # C3(0,0,2/3),C2(0,0,1/2),C2X 
symgen["P 62 2c (0 0 1)"]  = "R3Z:T(0.0,0.0,0.6666666666666667)_R2Z_R2X_"                          # C3(0,0,2/3),C2,C2X 
symgen["P 64 2c (0 0 -1)"] = "R3Z:T(0.0,0.0,0.3333333333333333)_R2Z_R2X_"                          # C3(0,0,1/3),C2,C2X 
symgen["P 6c 2c"]          = "R3Z_R2Z:T(0.0,0.0,0.5)_R2X_"                                         # C3,C2(0,0,1/2),C2X 
symgen["P 6 -2"]           = "R3Z_R2Z_MX_"                                                         # C3,C2,SigV 
symgen["P 6 -2c"]          = "R3Z_R2Z_MX:T(0.0,0.0,0.5)_"                                          # C3,C2,SigV(0,0,1/2)
symgen["P 6c -2"]          = "R3Z_R2Z:T(0.0,0.0,0.5)_MX_"                                          # C3,C2(0,0,1/2),SigV
symgen["P 6c -2c"]         = "R3Z_R2Z:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                           # C3,C2(0,0,1/2),SigV(0,0,SigV)
symgen["P -6 2"]           = "R3Z_MZ_R2(-1.0,1.0,0.0)_"                                            # C3,SigH,C2D
symgen["P -6c 2"]          = "R3Z_MZ(0,0,0.5)_R2(-1.0,1.0,0.0)_"                                   # C3,SigH(0,0,1/2),C2D
symgen["P -6 -2"]          = "R3Z_MZ_R2X_"                                                         # C3,SigH,C2X
symgen["P -6c -2c"]        = "R3Z_MZ(0,0,0.5)_R2X_"                                                # C3,SigH(0,0,1/2),C2X
symgen["-P 6 2"]           = "R3Z_R2Z_R2X_MX_"                                                     # C3,C2,C2X,SigV
symgen["-P 6 2c"]          = "R3Z_R2Z_R2X:T(0.0,0.0,0.5)_MX:T(0.0,0.0,0.5)_"                       # C3,C2,C2X(0,0,1/2),SigV(0,0,1/2)
symgen["-P 6c 2"]          = "R3Z_R2Z:T(0.0,0.0,0.5)_R2X:T(0.0,0.0,0.5)_MX_"                       # C3,C2(0,0,1/2),C2X(0,0,1/2),SigV
symgen["-P 6c 2c"]         = "R3Z_R2Z:T(0.0,0.0,0.5)_M(1.0,-1.0,0.0):T(0.0,0.0,0.5)_"              # C3,C2(0,0,1/2),SigV(0,0,1/2)
symgen[195] = "R2Z_R2X_R3D_"                                                                       # C2,C2X,C3XYZ
symgen[196] = "R2Z_R2X_R3D_"                                                                       # C2,C2X,C3XYZ
symgen[197] = "R2Z_R2X_R3D_"                                                                       # C2,C2X,C3XYZ
symgen[198] = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_R3D_"                                         # C2(1/2,0,1/2),C2X(1/2,1/2,0),C3XYZ
symgen[199] = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_R3D_"                                         # C2(1/2,0,1/2),C2X(1/2,1/2,0),C3XYZ
symgen[200] = "R2Z_R2X_MX_R3D_"                                                                    # C2,C2X,SigV,C3XYZ
symgen[201] = "R2Z_R2X_MX:T(0.5,0.5,0.5)_R3D_"                                                     # C2,C2X,SigV(1/2,1/2,1/2),C3XYZ
symgen[202] = "R2Z_R2X_MX_R3D_"                                                                    # C2,C2X,SigV,C3XYZ
symgen[203] = "R2Z_R2X_MX:T(0.25,0.25,0.25)_R3D_"                                                  # C2,C2X,SigV(1/4,1/4,1/4),C3XYZ
symgen[204] = "R2Z_R2X_R3D_MZ_"                                                                    # C2,C2X,C3XYZ,SigH
symgen[205] = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_MX:T(0.5,0.5,0.0)_R3D_"                       # C2(1/2,0,1/2),C2X(1/2,1/2,0),SigV(1/2,1/2,0),C3XYZ
symgen[206] = "R2Z:T(0.0,0.5,0.0)_R2X:T(0.0,0.0,0.5)_R3D_MZ:T(0.0,0.5,0.0)_"                       # C2(0,1/2,0),C2X(0,0,1/2),C3XYZ,SigH(0,1/2,0)
symgen[207] = "R2Z_R2X_R3D_R4Z_"                                                                   # C2,C2X,C3XYZ,C4
symgen[208] = "R2Z_R2X_R3D_R4Z:T(0.5,0.5,0.5)_"                                                    # C2,C2X,C3XYZ,C4(1/2,1/2,1/2)
symgen[209] = "R2Z_R2X_R3D_R4Z_"                                                                   # C2,C2X,C3XYZ,C4
symgen[210] = "R2Z_R2X_R3D_R4Z:T(0.25,0.25,0.25)_"                                                 # C2,C2X,C3XYZ,C4(1/4,1/4,1/4)
symgen[211] = "R2Z_R2X_R3D_R4Z_"                                                                   # C2,C2X,C3XYZ,C4
symgen[212] = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_R3D_R4Z:T(0.75,0.25,0.75)_"                   # C2(1/2,0,1/2),C2X(1/2,1/2,0),C3XYZ,C4(3/4,1/4,3/4)
symgen[213] = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_R3D_R4Z:T(0.25,0.75,0.25)_"                   # C2(1/2,0,1/2),C2X(1/2,1/2,0),C3XYZ,C4(1/4,3/4,1/4)
symgen[214] = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_R3D_R4Z:T(0.75,0.25,0.75)_"                   # C2(1/2,0,1/2),C2X(1/2,1/2,0),C3XYZ,C4(3/4,1/4,3/4)
symgen[215] = "R2Z_R2X_R3D_M(-1.0,1.0,0.0)_"                                                       # C2,C2X,C3XYZ,SigD
symgen[216] = "R2Z_R2X_R3D_M(-1.0,1.0,0.0)_"                                                       # C2,C2X,C3XYZ,SigD
symgen[217] = "R2Z_R2X_R3D_M(-1.0,1.0,0.0)_"                                                       # C2,C2X,C3XYZ,SigD
symgen[218] = "R2Z_R2X_R3D_M(-1.0,1.0,0.0):T(0.5,0.5,0.5)_"                                        # C2,C2X,C3XYZ,SigD(1/2,1/2,1/2)
symgen[219] = "R2Z_R2X_R3D_M(-1.0,1.0,0.0):T(0.0,0.0,0.5)_"                                        # C2,C2X,C3XYZ,SigD(0,0,1/2)
symgen[220] = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_R3D_M(-1.0,1.0,0.0):T(0.25,0.25,0.25)_"       # C2(1/2,0,1/2),C2X(1/2,1/2,0),C3XYZ,SigD(1/4,1/4,1/4)
symgen[221] = "R2Z_R2X_MX_R3D_M(-1.0,1.0,0.0)_"                                                    # C2,C2X,SigV,C3XYZ,SigD
symgen[222] = "R2Z_R2X_MX:T(0.5,0.5,0.5)_R3D_M(-1.0,1.0,0.0):T(0.5,0.5,0.5)_"                      # C2,C2X,SigV(1/2,1/2,1/2),C3XYZ,SigD(1/2,1/2,1/2)
symgen[223] = "R2Z_R2X_MX_R3D_M(-1.0,1.0,0.0):T(0.5,0.5,0.5)_"                                     # C2,C2X,SigV,C3XYZ,SigD(1/2,1/2,1/2)
symgen[224] = "R2Z_R2X_R3D_M(-1.0,1.0,0.0)_I:T(0.5,0.5,0.5)_"                                      # C2,C2X,C3XYZ,SigD,I(1/2,1/2,1/2)
symgen[225] = "R2Z_R2X_MX_R3D_M(-1.0,1.0,0.0)_"                                                    # C2,C2X,SigV,C3XYZ,SigD
symgen[226] = "R2Z_R2X_MX_R3D_M(-1.0,1.0,0.0):T(0.0,0.0,0.5)_"                                     # C2,C2X,SigV,C3XYZ,SigD(0,0,1/2)
symgen[227] = "R2Z:T(0.25,0.25,0.0)_R2X:T(0.0,0.25,0.25)_MX:T(0.0,0.75,0.75)_R3D_M(-1.0,1.0,0.0)_" # C2(1/4,1/4,0),C2X(0,1/4,1/4),SigV(0,3/4,3/4),C3XYZ,SigD
symgen[228] = "R2Z:T(0.25,0.25,0.0)_R2X:T(0.0,0.25,0.25)_R3D_M(-1.0,1.0,0.0):T(0.0,0.0,0.5)_I_"    # C2(1/4,1/4,0),C2X(0,1/4,1/4),C3XYZ,SigD(0,0,1/2),I
symgen["-I 4 2 3"] = "R2Z_R2X_MX_R3D_M(-1.0,1.0,0.0)_"                                                    # C2,C2X,SigV,C3XYZ,SigD
symgen[230] = "R2Z:T(0.5,0.0,0.5)_R2X:T(0.5,0.5,0.0)_R3D_M(-1.0,1.0,0.0):T(0.25,0.25,0.25)_I_"     # C2(1/2,0,1/2),C2X(1/2,1/2,0),C3XYZ,SigD(1/4,1/4,1/4),I

def retr_symmetry_generators(struct,ini):
    """
    Retrieve the symmetry generators.

    The Hall symbol is extracted from the structure instance. With this 
    symbol the corresponding string of symmetry generators is looked up in 
    the symgen dictionary. This string is added to the dictionary passed in 
    and the updated dictionary is returned.
    """
    #hall = struct.spacegroup_hall()
    ini["symgen"] = struct.get_symmetry_generators()
    return ini

def retr_symmetry_operations(struct,ini):
    """
    Retrieve the symmetry operations.

    The symmetry operations are extracted from the structure instance. 
    A list of rotation matrices and translation vectors is returned.
    """
    ini["symgen"] = struct.get_symmetry_operations()
    return ini

def retr_spacegroup_number(struct,ini):
    """
    Retrieve the space group number.

    When choosing real space grids some space groups imply restrictions
    on the number of points along the different lattice vectors.
    To impose the correct restrictions we need to know the space group
    number of the current system. This function add this number to
    the dictionary.
    """
    ini["spacegroup"] = struct.get_space_group_number()
    return ini
