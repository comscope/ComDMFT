# Copyright (C) 2006 Huub van Dam, Science and Technology Facilities Council,
# Daresbury Laboratory.
# All rights reserved.
#
# Developed by:        Huub van Dam
#                      Science and Technology Facilities Council
#                      Daresbury Laboratory
#                      Computational Science and Engineering Department
#                      Computational Chemistry Group
#                      http://www.cse.clrc.ac.uk/ccg
#
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"),
# to deal with the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions: 
#
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimers. 
# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution. 
# Neither the names of the Science and Technology Facilities Council,
# Daresbury Laboratory, the Computational Science and Engineering Department,
# the Computational Chemistry Group, nor the names of its contributors may be
# used to endorse or promote products derived from this Software without
# specific prior written permission. 
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS WITH THE SOFTWARE. 

import string
import toldiff_lcs
import toldiff_tokens

def lcs_to_change(lcs,ref,Nb,Ne,dat,Mb,Me,change,feps,ieps,itol_scale,ftol_scale,ctol_scale):
  """This routine takes the snake list, the reference file, the data file,
     and the change section from the tolerance file. Based on this information
     it merges any new tolerances into the tolerance data structures.
     The new change tolerance data structures are returned.
     """
  nsnake = len(lcs)
  isnake = 0
  while (isnake < nsnake):
    (xi2,yj2,xi3,yj3,itype) = lcs[isnake]
    #
    # deal with changed lines...
    #
    if itype == 3:
      #
      # The lines in the snake do not match exactly and not after considering
      # the tolerances. They are paired however so as to maximize similarity.
      # So we need to update the tolerances to make them a type 2 match.
      #
      i = 0
      while (yj2+i <= yj3):
        tmp_ref = ref.token[yj2+i]
        tmp_dat = dat.token[xi2+i]
        tmp_tol = ""
        if change.has_key(yj2+i):
          tmp_tol = change[yj2+i]
        tmp_tol = toldiff_tokens.tolerance(tmp_ref,tmp_dat,tmp_tol,feps,ieps,itol_scale,ftol_scale,ctol_scale)
        if tmp_tol != "":
          change[yj2+i] = tmp_tol
        i = i + 1
    isnake = isnake + 1
  return change

def lcs_to_add(lcs,ref,Nb,Ne,dat,Mb,Me,add):
  """This routine takes the snake list, the reference file, the data file,
     and the add section from the tolerance file. Based on this information it
     merges any new tolerances into the add tolerance data structure.
     The new add tolerance data structure are returned.
     """
  nsnake = len(lcs)
  isnake = 0
  xi1 = Mb-1
  yj1 = Nb-1
  if nsnake == 0:
    #
    # deal with added lines
    #
    if yj1 == 0:
      lnt = 0
      tnt = 0
    else:
      (tr,dr,lnt,tnt) = ref.token[yj1]
    Nadd = Me-Mb+1
    if Nadd > 0:
      if add.has_key(yj1):
        (ntok,lnt,tnt) = add[yj1]
        if ntok < Nadd:
          add[yj1] = (Nadd,lnt,tnt)
      else:
        add[yj1] = (Nadd,lnt,tnt)
  else:
    while (isnake < nsnake):
      (xi2,yj2,xi3,yj3,itype) = lcs[isnake]
      #
      # deal with added lines
      #
      if yj1 == 0:
        lnt = 0
        tnt = 0
      else:
        (tr,dr,lnt,tnt) = ref.token[yj1]
      Nadd = xi2-(xi1+1)
      if Nadd > 0:
        if add.has_key(yj1):
          (ntok,lnt,tnt) = add[yj1]
          if ntok < Nadd:
            add[yj1] = (Nadd,lnt,tnt)
        else:
          add[yj1] = (Nadd,lnt,tnt)
      isnake = isnake + 1
      xi1 = xi3
      yj1 = yj3
    if yj1 == 0:
      lnt = 0
      tnt = 0
    else:
      (tr,dr,lnt,tnt) = ref.token[yj1]
    Nadd = Me-xi1
    if Nadd > 0:
      if add.has_key(yj1):
        (ntok,lnt,tnt) = add[yj1]
        if ntok < Nadd:
          add[yj1] = (Nadd,lnt,tnt)
      else:
        add[yj1] = (Nadd,lnt,tnt)
  return add

def lcs_to_delete(lcs,ref,Nb,Ne,dat,Mb,Me,delete):
  """This routine takes the snake list, the reference file, the data file,
     and the delete section from the tolerance file. Based on this information
     it merges any new tolerances into the delete tolerance data structure.
     The new delete tolerance data structure is returned.
     """
  nsnake = len(lcs)
  isnake = 0
  xi1 = Mb-1
  yj1 = Nb-1
  if nsnake == 0:
    #
    # deal with deleted lines
    #
    i = Nb
    while (i <= Ne):
      (tr,dr,lnt,tnt) = ref.token[i]
      delete[i] = (lnt,tnt)
      i = i + 1
  else:
    while (isnake < nsnake):
      (xi2,yj2,xi3,yj3,itype) = lcs[isnake]
      #
      # deal with deleted lines
      #
      i = yj1+1
      while (i < yj2):
        (tr,dr,lnt,tnt) = ref.token[i]
        delete[i] = (lnt,tnt)
        i = i + 1
      isnake = isnake + 1
      xi1 = xi3
      yj1 = yj3
    i = yj1+1
    while (i <= Ne):
      (tr,dr,lnt,tnt) = ref.token[i]
      delete[i] = (lnt,tnt)
      i = i + 1
  return delete


