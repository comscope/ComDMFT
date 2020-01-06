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

def transfer_tol(lcs,Nb,Ne,Mb,Me,chg_txt,add_txt,del_txt):
  """Take the snake list and transfer all tolerances from the type 1 and 2
     snakes to a new tolerance file.
     """
  chg_new = {}
  add_new = {}
  del_new = {}
  nsnake  = len(lcs)
  isnake  = 0
  while (isnake < nsnake):
    (xi2,yj2,xi3,yj3,itype) = lcs[isnake]
    i = 0
    while (yj2+i <= yj3):
      if chg_txt.has_key(yj2+i):
        chg_new[xi2+i] = chg_txt[yj2+i]
      if add_txt.has_key(yj2+i):
        add_new[xi2+i] = add_txt[yj2+i]
      if del_txt.has_key(yj2+i):
        del_new[xi2+i] = del_txt[yj2+i]
      i = i + 1
    isnake = isnake + 1
  #
  # In case an external diff program is used there are some potential funnies
  # in that a potential leading snake and a potential final snake are missing.
  # The tolerances with these snakes need to be transfered nevertheless.
  #
  # Do the change tolerances for these implied snakes first
  #
  offset = 0
  if nsnake == 0:
    if Nb == 0 and Mb == 0 and Ne == -1 and Me == -1:
      #
      # The two files are actually identical...
      # 
      end_ref_1st    = -1
      begin_ref_last = -1
      end_dat_1st    = -1
      begin_dat_last = -1
      #
    else:
      if Nb <= Ne:
        #
        # The first implied snake terminates at Nb-1 and
        # the last implied snake starts at Ne+1
        #
        end_ref_1st    = Nb-1
        begin_ref_last = Ne+1
        #
      else:
        #
        # The first implied snake terminates at Nb and
        # the last implied snake starts at Nb+1
        #
        offset         = +1
        end_ref_1st    = Nb
        begin_ref_last = Nb+1
        #
      if Mb <= Me:
        #
        # The first implied snake terminates at Mb-1 and
        # the last implied snake starts at Me+1
        #
        end_dat_1st    = Mb-1
        begin_dat_last = Me+1
        #
      else:
        #
        # The first implied snake terminates at Mb and
        # the last implied snake starts at Mb+1
        #
        offset         = -1
        end_dat_1st    = Mb
        begin_dat_last = Mb+1
        #
  else:
    (xi2,yj2,xi3,yj3,itype) = lcs[0]
    if Nb < yj2:
      #
      # The first implied snake terminates at Nb-1
      #
      end_ref_1st = Nb-1
      #
    else:
      #
      # The first implied snake terminates at Nb
      #
      end_ref_1st = Nb
      #
    if Mb < xi2:
      #
      # The first implied snake terminates at Mb-1
      #
      end_dat_1st = Mb-1
      #
    else:
      #
      # The first implied snake terminates at Mb
      #
      end_dat_1st = Mb
      #
    (xi2,yj2,xi3,yj3,itype) = lcs[nsnake-1]
    if yj3 < Ne:
      #
      # The last implied snake starts at Ne+1
      #
      begin_ref_last = Ne+1
      #
    else:
      #
      # The last implied snake starts at yj3+1
      #
      offset         = +1
      begin_ref_last = yj3+1
      #
    if xi3 < Me:
      #
      # The last implied snake starts at Me+1
      #
      begin_dat_last = Me+1
      #
    else:
      #
      # The last implied snake starts at xi3+1
      #
      offset         = -1
      begin_dat_last = xi3+1
      #
  #
  # Now transfer the change tolerances
  #
  Tkeys = chg_txt.keys()
  while (len(Tkeys) > 0):
    key = Tkeys.pop()
    if key <= end_ref_1st:
      chg_new[key] = chg_txt[key]
    if begin_ref_last <= key:
      chg_new[key-begin_ref_last+begin_dat_last+offset] = chg_txt[key]
  #
  # Now transfer the add tolerances
  #
  Tkeys = add_txt.keys()
  while (len(Tkeys) > 0):
    key = Tkeys.pop()
    if key <= end_ref_1st:
      add_new[key] = add_txt[key]
    if begin_ref_last <= key:
      add_new[key-begin_ref_last+begin_dat_last+offset] = add_txt[key]
  #
  # Now transfer the del tolerances
  #
  Tkeys = del_txt.keys()
  while (len(Tkeys) > 0):
    key = Tkeys.pop()
    if key <= end_ref_1st:
      del_new[key] = del_txt[key]
    if begin_ref_last <= key:
      del_new[key-begin_ref_last+begin_dat_last+offset] = del_txt[key]
  return (chg_new,add_new,del_new)

