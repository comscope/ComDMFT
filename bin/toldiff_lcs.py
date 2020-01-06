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

import sys
import copy
import string
import toldiff_tokens

true = (0 == 0)
false = not true

class search_path:
  """The search_path class is used to record the progress of a search
     path through the two files. There are three components that should be 
     tracked
     1. A list of snakes each snake is recorded by its starting and end point.
     2. As the last point is not recorded if there is no match the last point
        explored in the search path needs to be recorded separately.
     3. For the purpose of eliminating failing paths it is useful to keep 
        track of the length longest common subsequence sofar.
     """
  def __init__(self):
    """Initialize the search path.
       1. Create an empty snake list.
       2. Create an initial last point (although this will require resetting
          for backward searches).
       3. Set the length sofar to zero.
       """
    self.snakelist = []
    self.lastpoint = (0,0)
    self.lcslength = 0

  def set_lastpoint(self,i,j):
    """Set the last point searched to (i,j)"""
    self.lastpoint = (i,j)

  def increase_length(self,n=1):
    """Increase the LCS length sofar by n"""
    self.lcslength = self.lcslength + n

  def add_snake(self,i1,j1,i2,j2,itype):
    """Add a new snake defined by points (i1,j1) and (i2,j2)
       also store the type of the snake."""
    self.snakelist.append((i1,j1,i2,j2,itype))

class search_path_linked:
  """The search_path class is used to record the progress of a search
     path through the two files. However the search path class as implemented
     has as a severe disadvantage that if a search splits into two different
     directions a new search has to be created that is a copy of the one it 
     split off from. This copy operation is extremely expensive. 

     So in the search_path_linked implementation the snake list is stored as
     a linked list instead. Every snake links back to the previous snake in the
     search path. Thus avoiding the whole snakelist to be copied if a search
     splits. 

     To generate the final snake list that is returned from the search the
     linked snake list will need to be transformed.

     The linked search path has three components that should be tracked
     1. A linked list of snakes each snake is recorded by its starting and end
        point.
     2. As the last point is not recorded if there is no match the last point
        explored in the search path needs to be recorded separately.
     3. For the purpose of eliminating failing paths it is useful to keep 
        track of the length longest common subsequence sofar.
     """
  def __init__(self):
    """Initialize the search path.
       1. Create an empty linked snake list.
       2. Create an initial last point (although this will require resetting
          for backward searches).
       3. Set the length sofar to zero.
       """
    self.snakelist = None
    self.lastpoint = (0,0)
    self.lcslength = 0

  def set_lastpoint(self,i,j):
    """Set the last point searched to (i,j)"""
    self.lastpoint = (i,j)

  def increase_length(self,n=1):
    """Increase the LCS length sofar by n"""
    self.lcslength = self.lcslength + n

  def add_snake(self,i1,j1,i2,j2,itype,previous_snake=None):
    """Add a new snake defined by points (i1,j1) and (i2,j2)
       also store the type of the snake."""
    if (previous_snake == None):
      self.snakelist = (i1,j1,i2,j2,itype,self.snakelist)
    else:
      self.snakelist = (i1,j1,i2,j2,itype,previous_snake)

  def get_previous_snake(self):
    """Return the one but last snake. If there is no previous snake return
       zero for the coordinates and type."""
    if self.snakelist == None:
      snake = (0,0,0,0,0)
    else:
      (i1,j1,i2,j2,itype,some) = self.snakelist
      snake = (i1,j1,i2,j2,itype)
    return snake


def min(a,b):
  """Return the minimum of two values"""
  if (a<b):
    return a
  else:
    return b

def max(a,b):
  """Return the maximum of two values"""
  if (a>b):
    return a
  else:
    return b

def transform_snake_list(linked_snake):
  """Transforms a linked snake list to a canonical snake list, i.e.
     from something of the form:

       (x1,y1,x2,y2,itype,previous_snake)

     to something of the form:

       [ ..., (x1,y1,x2,y2,ityp), ... ]

     This is needed to go from the linked list snake search to a snake list
     that is more suitable to other parts of the program."""

  lcs = []
  while (linked_snake):
    (x1,y1,x2,y2,itype,linked_snake) = linked_snake
    lcs.insert(0,(x1,y1,x2,y2,itype))
  return lcs

def trim_snakes(lcs,ref,Nb,Ne,dat,Mb,Me):
  """Previously found matches can cause problems if they are not optimal.
     In such a case sticking with the matches as found prevents subsequent
     more advanced diff routines from recovering from an early sub-optimal
     choice. To counter this some snakes and pseudo-snakes are trimmed down,
     or equivalently the mismatches are expanded to whole lines.
     Note that it is clear that if a mismatch involves the same number of tokens
     of the reference file as of the data file then expanding it does not bring
     any benefits. Only when different numbers of tokens and hence adds or
     deletes are implied does this trick help.
     The process is:
     1. Merge subsequent snakes to build a list in which each pair of
        snakes is separated by a non-empty section of mismatching tokens.
     2. Replace the snake list by a list of mismatching sections.
     3. Go over the mismatching sections and expand each one in which the
        number of reference tokens is not the same as the number of data tokens.
        Eliminate mismatch sections that are sub-sections of larger ones along
        the way.
     4. Go back from the list of mismatching sections to a snake list.
     The routine returns the revised snake list.
     """
  true  = (0 == 0)
  false = (0 != 0)
  #
  # Collapse the snake list by merging adjacent snakes.
  #
  #DEBUG
  #print "--- trim_snakes in ---"
  #print_lcs(sys.stdout,lcs)
  #DEBUG
  nsnake = len(lcs)
  isnake = 0
  if nsnake > 0:
    lcs_tmp = []
    (xi1,yj1,xi2,yj2,itype) = lcs[isnake]
    isnake = isnake + 1
    while (isnake < nsnake):
      (xi3,yj3,xi4,yj4,itype) = lcs[isnake]
      isnake = isnake + 1
      if (xi2+1 == xi3 and yj2+1 == yj3):
        #
        # This snake continues from the previous one so merge the two.
        #
        xi2 = xi4
        yj2 = yj4
        #
      else:
        #
        # This snake is separated from the previous one so store the 
        # previous one and restart the merge procedure.
        #
        lcs_tmp.append((xi1,yj1,xi2,yj2,itype))
        xi1 = xi3
        yj1 = yj3
        xi2 = xi4
        yj2 = yj4
    #
    # Store the last snake.
    #
    lcs_tmp.append((xi1,yj1,xi2,yj2,itype))
    lcs = lcs_tmp
  #
  # Replace the snake list by a list of mismatch sections
  #
  nsnake = len(lcs)
  isnake = 0
  lcs_tmp = []
  txi = Mb
  tyj = Nb
  if nsnake > 0:
    while (isnake < nsnake):
      (xi1,yj1,xi2,yj2,itype) = lcs[isnake]
      if (xi1 > txi or yj1 > tyj):
        lcs_tmp.append((txi,tyj,xi1-1,yj1-1,4))
      txi = xi2+1
      tyj = yj2+1
      isnake = isnake + 1
    if (txi <= Me or tyj <= Ne):
      lcs_tmp.append((txi,tyj,Me,Ne,4))
  else:
    lcs_tmp.append((Mb,Nb,Me,Ne,4))
  lcs = lcs_tmp
  #
  # Expand certain differences
  # 
  nsnake = len(lcs)
  isnake = 0
  tsnake = 0
  lcs_tmp = []
  while (isnake < nsnake):
    (xi3,yj3,xi4,yj4,itype) = lcs[isnake]
    isnake = isnake + 1
    rsnake = tsnake
    if xi4-xi3 != yj4-yj3:
      #
      # Expand the difference out to a set of whole lines
      #
      # Ensure that the end point is at least on the same line as the 
      # start point otherwise the expanded difference may not make any sense.
      # There are at least two tricky bits:
      # 1. Assume extra tokens have been added at the beginning of a line.
      #    Yj4 will then be yj3-1 and point to the line before yj3, whereas
      #    xi3 and xi4 will point to the same line. If nothing is done then
      #    after expansion yj3 and yj4 will still be the same but xi3 and xi4
      #    will have been expanded out. The additional differences are now
      #    unrecoverable.
      # 2. Assume that xi3 and xi4 cover the addition of a full set of lines
      #    and yj4 == yj3-1. Simply adjusting yj4 leads to an extra reference
      #    token being added to the difference but no data token. Again an
      #    unintended difference would occur.
      # To address both of the above one has to increase both end points 
      # simultaneously.
      #
      if xi3 > xi4:
        yj4 = yj4 + (xi3 - xi4)
        xi4 = xi3
      elif yj3 > yj4:
        xi4 = xi4 + (yj3 - yj4)
        yj4 = yj3
      #
      lxi3 = toldiff_tokens.tokenno2lineno(dat,xi3)
      xi3  = toldiff_tokens.lineno2tokenno(dat,lxi3)
      lyj3 = toldiff_tokens.tokenno2lineno(ref,yj3)
      yj3  = toldiff_tokens.lineno2tokenno(ref,lyj3)
      lxi4 = toldiff_tokens.tokenno2lineno(dat,xi4)+1
      xi4  = toldiff_tokens.lineno2tokenno(dat,lxi4)-1
      lyj4 = toldiff_tokens.tokenno2lineno(ref,yj4)+1
      yj4  = toldiff_tokens.lineno2tokenno(ref,lyj4)-1
      tsnake = tsnake + 1
      lcs_tmp.append((xi3,yj3,xi4,yj4,4))
    else:
      #
      # Keep the difference the way it is
      #
      tsnake = tsnake + 1
      lcs_tmp.append((xi3,yj3,xi4,yj4,4))
    #
    # Check if the new mismatch section overlaps with another one
    #
    done = false
    while (rsnake > 0 and not done):
      (xi1,yj1,xi2,yj2,itype) = lcs_tmp[rsnake-1]
      if xi2 < xi3 and yj2 < yj3:
        done = true
      elif xi3 <= xi1 and yj3 <= yj1 and xi2 <= xi4 and yj2 <= yj4:
        #
        # rsnake is a sub set of the latest set so remove it
        #
        rsnake = rsnake - 1
        tsnake = tsnake - 1
        tt = lcs_tmp.pop(rsnake)
      elif xi1 <= xi3 and yj1 <= yj3 and xi4 <= xi2 and yj4 <= yj2:
        #
        # tsnake is a sub set of rsnake so remove it
        #
        tsnake = tsnake - 1
        tt = lcs_tmp.pop(tsnake)
        done = true
      elif xi1 <= xi3 and xi2 >= xi3 and xi2 <= xi4:
        #
        # overlapping mismatches merge both
        #
        lcs_tmp[rsnake-1] = (xi1,min(yj1,yj3),xi4,max(yj2,yj4),itype)
        tt = lcs_tmp.pop(tsnake-1)
        rsnake = rsnake - 1
        tsnake = tsnake - 1
        xi3 = xi1
        yj3 = min(yj1,yj3)
        #
      elif yj1 <= yj3 and yj2 >= yj3 and yj2 <= yj4:
        #
        # overlapping mismatches merge both
        #
        lcs_tmp[rsnake-1] = (min(xi1,xi3),yj1,max(xi2,xi4),yj4,itype)
        tt = lcs_tmp.pop(tsnake-1)
        rsnake = rsnake - 1
        tsnake = tsnake - 1
        xi3 = min(xi1,xi3)
        yj3 = yj1
        #
      else:
        print "backtrack: this is odd..."
        print "(yj1,yj2)<->(xi1,xi2) : ("+str(yj1)+","+str(yj2)+")<->("+str(xi1)+","+str(xi2)+")"
        print "(yj3,yj4)<->(xi3,xi4) : ("+str(yj3)+","+str(yj4)+")<->("+str(xi3)+","+str(xi4)+")"
        sys.exit(10)
  lcs = lcs_tmp
  #DEBUG
  #print "--- trim_snakes different ---"
  #print_lcs(sys.stdout,lcs)
  #DEBUG
  #
  # Convert from mismatches back to a snake list
  #
  nsnake = len(lcs)
  isnake = 0
  lcs_tmp = []
  if nsnake == 0:
    lcs_tmp.append((Mb,Nb,Me,Ne,1))
  else:
    xi0 = Mb
    yj0 = Nb
    while (isnake < nsnake):
      (xi1,yj1,xi2,yj2,itype) = lcs[isnake]
      isnake = isnake + 1
      lcs_tmp.append((xi0,yj0,xi1-1,yj1-1,1))
      xi0 = xi2+1
      yj0 = yj2+1
    lcs_tmp.append((xi0,yj0,Me,Ne,1))
  lcs = lcs_tmp
  #DEBUG
  #print "--- trim_snakes out ---"
  #print_lcs(sys.stdout,lcs)
  #DEBUG
  return lcs
     
def trim_snakes_old(lcs,ref,Nb,Ne,dat,Mb,Me):
  """Previously found matches can cause problems if they are not optimal.
     In such a case sticking with the matches as found prevents subsequent
     more advanced diff routines from recovering from an early sub-optimal
     choice. To counter this all snakes and pseudo-snakes are trimmed down
     such that they involve whole lines only.
     The process is:
     1. Merge subsequent snakes to build a list in which each pair of
        snakes is separated by a non-empty section of mismatching tokens.
     2. Trim each snake by increasing the starting point to the first token
        on the next line, and decreasing the end point to the last token on
        the previous line. If as a result the begin token exceeds the end
        token then eliminate the snake.
     The routine returns the revised snake list.
     """
  #
  # Collapse the snake list by merging adjacent snakes.
  #
  nsnake = len(lcs)
  isnake = 0
  if nsnake > 0:
    lcs_tmp = []
    (xi1,yj1,xi2,yj2,itype) = lcs[isnake]
    isnake = isnake + 1
    while (isnake < nsnake):
      (xi3,yj3,xi4,yj4,itype) = lcs[isnake]
      isnake = isnake + 1
      if (xi2+1 == xi3 and yj2+1 == yj3):
        #
        # This snake continues from the previous one so merge the two.
        #
        xi2 = xi4
        yj2 = yj4
        #
      else:
        #
        # This snake is separated from the previous one so store the 
        # previous one and restart the merge procedure.
        #
        lcs_tmp.append((xi1,yj1,xi2,yj2,itype))
        xi1 = xi3
        yj1 = yj3
        xi2 = xi4
        yj2 = yj4
    #
    # Store the last snake.
    #
    lcs_tmp.append((xi1,yj1,xi2,yj2,itype))
    lcs = lcs_tmp
  #
  # Trim the snakes to precisely matching lines.
  #
  nsnake = len(lcs)
  isnake = 0
  lcs_tmp = []
  txi = 0
  tyj = 0
  while (isnake < nsnake):
    (xi1,yj1,xi2,yj2,itype) = lcs[isnake]
    isnake = isnake + 1
    #
    # Move the starting point to the first token on the next line unless
    # the token is the first token on the current line.
    #
    lxi1 = toldiff_tokens.tokenno2lineno(dat,xi1)
    txi1 = toldiff_tokens.lineno2tokenno(dat,lxi1)
    lyj1 = toldiff_tokens.tokenno2lineno(ref,yj1)
    tyj1 = toldiff_tokens.lineno2tokenno(ref,lyj1)
    if txi1 != xi1 or tyj1 != yj1:
      xi1 = toldiff_tokens.lineno2tokenno(dat,lxi1+1)
      yj1 = toldiff_tokens.lineno2tokenno(ref,lyj1+1)
    #
    # Move the end point to the last token on the previous line unless
    # the token is the last token on the current line.
    #
    lxi2 = toldiff_tokens.tokenno2lineno(dat,xi2)
    txi2 = toldiff_tokens.lineno2tokenno(dat,lxi2+1)-1
    lyj2 = toldiff_tokens.tokenno2lineno(ref,yj2)
    tyj2 = toldiff_tokens.lineno2tokenno(ref,lyj2+1)-1
    if txi2 != xi2 or tyj2 != yj2:
      xi2 = toldiff_tokens.lineno2tokenno(dat,lxi2)-1
      yj2 = toldiff_tokens.lineno2tokenno(ref,lyj2)-1
    if xi1-1 <= xi2 and yj1-1 <= yj2 and (xi1 > txi or yj1 > tyj):
      #
      # There is a non-empty snake remaining so store it.
      #
      lcs_tmp.append((xi1,yj1,xi2,yj2,itype))
      txi = max(xi1,xi2)
      tyj = max(yj1,yj2)
      #
  lcs = lcs_tmp
  return lcs
     
def find_lcs1(ref,Nb,Ne,dat,Mb,Me):
  """Compares the data stored in 'dat' against the data in 'ref', 
     and returns the longest common subsequence (LCS) in 'lcs'. The LCS
     is stored as a list of snakes. A snake is a sequence of line pairs
     (Xi,Yj) to (Xi+p,Yj+p) where the lines X and Y in every pair match.
     Whatever happens between two snakes in a path is irrelevant.
     As this routine looks for exact matches it produces type 1 snakes.

     The algorithm used here is inspired by:

       E. Myers, 'An O(ND) Difference Algorithm and Its Variations'
       Algorithmica 1, 2 (1986), 251-266
       http://www.cs.arizona.edu/people/gene/PAPERS/diff.ps

     however I cannot guarantee that understood it well enough to reproduce
     the actual published algorithm.
     
     Huub van Dam, SciTech Daresbury Laboratory, June 2006.
     """

  lcs = { }
  # FP - Forward Pij
  # Records the maximum number of diagonal lines of all candidates paths that
  # passed through node (i,j). P is a dictionary with tuples (i,j) as keys and
  # the maximum number as data.
  FP = { }

  # FV - Forward search path vector
  # Stores the forwards search paths. 
  FV = { }
  # NF - counter for generating forward search path keys
  #NF = 1
  # 
  s = search_path_linked()
  s.set_lastpoint(Mb-1,Nb-1)
  FV[(Mb-1,Nb-1)] = s

  # flij - forward last i+j
  # foij - forward old i+j
  # lij is the smallest i+j of the end point of any search path in the current
  # pass.
  # oij is the value of lij of the previous pass. These values will be used
  # to eliminate any entries in P that are no longer nessecary.
  flij = Me+Ne
  foij = Mb+Nb

  # the length of the longest LCS sofar
  max_len_lcs = 0

  finished = (0 != 0)

  # D - the number of non-diagonal steps
  D = 0
  while (D < (Me+Ne) and not(finished)):
    #
    # Work on the forward searches
    #
    D = D + 1
    flij = Me+Ne
    Fkeys = FV.keys()
    Fkeys.sort()
    while (len(Fkeys) > 0):
      key = Fkeys.pop()
      s = FV[key]
      del FV[key]
      (xi,yj) = s.lastpoint
      opt_len = s.lcslength + min(Me-xi+1,Ne-yj+1)
      if (opt_len > max_len_lcs):
        #
        # There (still) is hope that this search can beat the best search sofar
        #
        #
        # First try whether we are onto a snake
        #
        xi1 = xi+1
        yj1 = yj+1
        if yj1 <= Ne and xi1 <= Me:
          (type,token,lineno,tokenno) = ref.token[yj1]
          reftok = (type,token)
          (type,token,lineno,tokenno) = dat.token[xi1]
          dattok = (type,token)
          if reftok == dattok:
            # yes, we are onto a snake
            xi2 = xi1 + 1
            yj2 = yj1 + 1
            while (yj2 <=Ne and xi2 <= Me):
              (type,token,lineno,tokenno) = ref.token[yj2]
              reftok = (type,token)
              (type,token,lineno,tokenno) = dat.token[xi2]
              dattok = (type,token)
              if reftok == dattok:
                xi2 = xi2 + 1
                yj2 = yj2 + 1
              else:
                break
            xi2 = xi2 - 1
            yj2 = yj2 - 1
            s.add_snake(xi1,yj1,xi2,yj2,1)
            s.increase_length(xi2-xi1+1)
            s.set_lastpoint(xi2,yj2)
            xi = xi2
            yj = yj2
            finished = (yj2 == Ne and xi2 == Me)
            if finished:
              lcs = transform_snake_list(s.snakelist)
              break
        #
        # update the maximum LCS length
        #
        max_len_lcs = max(max_len_lcs,s.lcslength)
        #
        # now explore the way forward, horizontal first
        #
        keep_horizontal = false
        xih = xi+1
        yjh = yj
        if xih <= Me:
          if FP.has_key((xih,yjh)):
            if FP[(xih,yjh)] < s.lcslength:
              keep_horizontal = true
          else:
            keep_horizontal = true
          if xih+yjh < flij:
            flij = xih+yjh
          finished = (yjh == Ne and xih == Me)
          if finished:
            lcs = transform_snake_list(s.snakelist)
            break
        #
        # now explore the vertical direction
        #
        keep_vertical = false
        xiv = xi
        yjv = yj+1
        if yjv <= Ne:
          if FP.has_key((xiv,yjv)):
            if FP[(xiv,yjv)] < s.lcslength:
              keep_vertical = true
          else:
            keep_vertical = true
          if xiv+yjv < flij:
            flij = xiv+yjv
          finished = (yjv == Ne and xiv == Me)
          if finished:
            lcs = transform_snake_list(s.snakelist)
            break
        if keep_vertical:
          if keep_horizontal:
            # Keeping both horizontal and vertical search direction
            # So generate a new search path
            sa = copy.copy(s)
            sa.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = sa
            FP[(xiv,yjv)] = sa.lcslength
            #
            s.set_lastpoint(xih,yjh)
            FV[(xih,yjh)] = s
            FP[(xih,yjh)] = s.lcslength
          else:
            # Keeping only the vertical search direction
            # So simply update the current search path
            s.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = s
            FP[(xiv,yjv)] = s.lcslength
        else:
          if keep_horizontal:
            # Keeping only the horizontal search direction
            # So simply update the current search path
            s.set_lastpoint(xih,yjh)
            FV[(xih,yjh)] = s
            FP[(xih,yjh)] = s.lcslength
          else:
            # Keeping neither the horizontal or vertical search direction
            # So remove the current path from the search list
            pass
      else:
        pass
    #
    # now tidy up FP
    #
    ij = foij - (Mb+Nb)
    while (ij <= flij - (Mb+Nb)):
      xi = Mb + ij
      yj = Nb
      while (xi >= Mb):
        if FP.has_key((xi,yj)):
          del FP[(xi,yj)]
        xi = xi - 1
        yj = yj + 1
      ij = ij + 1
    foij = flij
  return lcs

def tol_decode(tol_line,txt_line):
  """Decode a single line of the tolerance change section and mark the
     relevant characters in the text line. Then return the modified 
     text line.
     """
  clist = string.split(tol_line,',')
  c1 = 0
  c2 = len(clist)
  while (c1 < c2):
    p = clist[c1]
    (s1,s2) = string.split(p,':')
    s1 = int(s1)
    s2 = int(s2)
    s3 = len(txt_line)
    while (s3 <= s2):
      txt_line = txt_line + " "
      s3 = s3 + 1
    while (s1 <= s2):
      txt_line = txt_line[:s1] + "#" + txt_line[s1+1:]
      s1 = s1 + 1
    c1 = c1 + 1
  return txt_line

def tol_encode(txt_line):
  """We assume that the text line contains hash character in all places where
     differences should be tolerated. This function will find the hash 
     characters and encode their locations in the format of a tolerance 
     change section line.
     """
  #
  # s1,s2 - beginning and end of a tolerance
  # intol - whether we are in a tolerance part
  # ntol  - the number of tolerance parts found
  #
  tol_line = ""
  true  = (0==0)
  false = not true
  ntol  = 0
  i     = 0
  n     = len(txt_line)
  intol = false
  while (i < n):
    if txt_line[i] == "#":
      if intol:
        s2 = i
      else:
        s1 = i
        s2 = i
        intol = true
        ntol  = ntol + 1
    else:
      if intol:
        if ntol > 1:
          tol_line = tol_line + ","
        tol_line = tol_line + str(s1)+":"+str(s2)
        intol = false
      else:
        pass
    i = i + 1
  if intol:
    if ntol > 1:
      tol_line = tol_line + ","
    tol_line = tol_line + str(s1)+":"+str(s2)
  return tol_line


def tol_compare_token(change,ref,yj,dat,xi,feps,ieps):
  """Compares two tokens taking into account tolerable differences
     stored in 'tol' if any.
     """
  reftok = ref.token[yj]
  dattok = dat.token[xi]
  if change.has_key(yj):
    toltok = change[yj]
  else:
    toltok = ""
  result = toldiff_tokens.tokens_match(reftok,dattok,toltok,feps,ieps)
  return result


def find_lcs2(tol,ref,Nb,Ne,dat,Mb,Me,feps,ieps):
  """Compares the data stored in 'dat' against the data in 'ref', 
     and returns the longest common subsequence (LCS) in 'lcs'. The LCS
     is stored as a list of snakes. A snake is a sequence of line pairs
     (Xi,Yj) to (Xi+p,Yj+p) where the lines X and Y in every pair match.
     Whatever happens between two snakes in a path is irrelevant.
     In this particular routine the string comparison is modified based
     on the information held 'tol'. For every relevant line the tol
     dictionary holds a string of splices of characters where differences
     should be tolerated. As this routine uses a tolerant comparison it
     generates type 2 snakes.

     The algorithm used here is inspired by:

       E. Myers, 'An O(ND) Difference Algorithm and Its Variations'
       Algorithmica 1, 2 (1986), 251-266
       http://www.cs.arizona.edu/people/gene/PAPERS/diff.ps

     however I cannot guarantee that understood it well enough to reproduce
     the actual published algorithm.
     
     Huub van Dam, SciTech Daresbury Laboratory, June 2006.
     """

  lcs = { }
  # FP - Forward Pij
  # Records the maximum number of diagonal lines of all candidates paths that
  # passed through node (i,j). P is a dictionary with tuples (i,j) as keys and
  # the maximum number as data.
  FP = { }

  # FV - Forward search path vector
  # Stores the forwards search paths. 
  FV = { }
  # NF - counter for generating forward search path keys
  # 
  s = search_path_linked()
  s.set_lastpoint(Mb-1,Nb-1)
  FV[(Mb-1,Nb-1)] = s

  # flij - forward last i+j
  # foij - forward old i+j
  # lij is the smallest i+j of the end point of any search path in the current
  # pass.
  # oij is the value of lij of the previous pass. These values will be used
  # to eliminate any entries in P that are no longer nessecary.
  flij = Me+Ne
  foij = Mb+Nb

  # the length of the longest LCS sofar
  max_len_lcs = 0

  finished = (0 != 0)

  # D - the number of non-diagonal steps
  D = 0
  while (D < (Me+Ne) and not(finished)):
    #
    # Work on the forward searches
    #
    D = D + 1
    flij = Me+Ne
    Fkeys = FV.keys()
    Fkeys.sort()
    while (len(Fkeys) > 0):
      key = Fkeys.pop()
      s = FV[key]
      del FV[key]
      (xi,yj) = s.lastpoint
      opt_len = s.lcslength + min(Me-xi+1,Ne-yj+1)
      if (opt_len > max_len_lcs):
        #
        # There (still) is hope that this search can beat the best search sofar
        #
        # First try whether we are onto a snake
        #
        xi1 = xi+1
        yj1 = yj+1
        if yj1 <= Ne and xi1 <= Me:
          # line comparison 1
          if tol_compare_token(tol,ref,yj1,dat,xi1,feps,ieps):
            # yes, we are onto a snake
            xi2 = xi1 + 1
            yj2 = yj1 + 1
            while (yj2 <=Ne and xi2 <= Me):
              # line comparison 2
              if tol_compare_token(tol,ref,yj2,dat,xi2,feps,ieps):
                xi2 = xi2 + 1
                yj2 = yj2 + 1
              else:
                break
            xi2 = xi2 - 1
            yj2 = yj2 - 1
            s.add_snake(xi1,yj1,xi2,yj2,2)
            s.increase_length(xi2-xi1+1)
            s.set_lastpoint(xi2,yj2)
            xi = xi2
            yj = yj2
            finished = (yj2 == Ne and xi2 == Me)
            if finished:
              lcs = transform_snake_list(s.snakelist)
              break
        #
        # update the maximum LCS length
        #
        max_len_lcs = max(max_len_lcs,s.lcslength)
        #
        # now explore the way forward, horizontal first
        #
        keep_horizontal = false
        xih = xi+1
        yjh = yj
        if xi1 <= Me:
          if FP.has_key((xih,yjh)):
            if FP[(xih,yjh)] < s.lcslength:
              keep_horizontal = true
          else:
            keep_horizontal = true
          if xih+yjh < flij:
            flij = xih+yjh
          finished = (yjh == Ne and xih == Me)
          if finished:
            lcs = transform_snake_list(s.snakelist)
            break
        #
        # now explore the vertical direction
        #
        keep_vertical = false
        xiv = xi
        yjv = yj+1
        if yjv <= Ne:
          if FP.has_key((xiv,yjv)):
            if FP[(xiv,yjv)] < s.lcslength:
              keep_vertical = true
          else:
            keep_vertical = true
          if xiv+yjv < flij:
            flij = xiv+yjv
          finished = (yjv == Ne and xiv == Me)
          if finished:
            lcs = transform_snake_list(s.snakelist)
            break
        if keep_vertical:
          if keep_horizontal:
            # Keeping both horizontal and vertical search directions
            # So generate a new search path
            sa = copy.copy(s)
            sa.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = sa
            FP[(xiv,yjv)] = sa.lcslength
            #
            s.set_lastpoint(xih,yjh)
            FV[(xih,yjh)] = s
            FP[(xih,yjh)] = s.lcslength
          else:
            # Keeping only the vertical search direction
            # So simply update the current search path
            s.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = s
            FP[(xiv,yjv)] = s.lcslength
        else:
          if keep_horizontal:
            # Keeping only the horizontal search direction
            # So simply update the current search path
            s.set_lastpoint(xih,yjh)
            FV[(xih,yjh)] = s
            FP[(xih,yjh)] = s.lcslength
          else:
            # Keeping neither the horizontal or vertical search direction
            # So simply let the current search path die
            pass
    #
    # now tidy up FP
    #
    ij = foij - (Mb+Nb)
    while (ij <= flij - (Mb+Nb)):
      xi = Mb + ij
      yj = Nb
      while (xi >= Mb):
        if FP.has_key((xi,yj)):
          del FP[(xi,yj)]
        xi = xi - 1
        yj = yj + 1
      ij = ij + 1
    foij = flij
  return lcs


def num_same_token(tol,ref,yj,dat,xi,feps,ieps):
  """Compares two lines of text taking into account tolerable differences
     stored in 'tol' if any and computing the number of characters that 
     are the same in the two lines.
     """
  if tol.has_key(yj):
    toltok = tol[yj]
  else:
    toltok = ""
  reftok = ref.token[yj]
  dattok = dat.token[xi]
  result = toldiff_tokens.compare_tokens(reftok,dattok,toltok,feps,ieps)
  return result

def find_lcs3(tol,ref,Nb,Ne,dat,Mb,Me,feps,ieps):
  """Compares the data stored in 'dat' against the data in 'ref', 
     and returns a common subsequence where lines are matched up so as
     to maximise the number of equal characters. The common subsequence
     is stored as a list of snakes. As the matching is based on a very
     much relaxed line comparison it generates type 3 snakes.

     The algorithm used here is inspired by:

       E. Myers, 'An O(ND) Difference Algorithm and Its Variations'
       Algorithmica 1, 2 (1986), 251-266
       http://www.cs.arizona.edu/people/gene/PAPERS/diff.ps

     however as the algorithm has to find a maximum number of matching 
     characters the algorithm has to be substantially modified. In particular
     the greedy approach does not work, neither does the first across the line
     criterion.

     The algorithm we use is to progress in three directions from every point.
     1. Horizontally,
     2. Vertically
     3. Diagonally
     As diagonal moves are the only ones that can up the similarity score it
     is expected that they are most likely to survive. Therefore we use the
     deepcopy for the horizontal and vertical moves. Finally, as there are no
     truely matching lines all type 3 snakes are only of length 1.
     
     Huub van Dam, SciTech Daresbury Laboratory, June 2006.
     """

  lcs_t = None
  # FP - Forward Pij
  # Records the maximum number of diagonal lines of all candidates paths that
  # passed through node (i,j). P is a dictionary with tuples (i,j) as keys and
  # the maximum number as data.
  FP = { }

  # FV - Forward search path vector
  # Stores the forwards search paths. 
  FV = { }
  # 
  s = search_path_linked()
  s.set_lastpoint(Mb-1,Nb-1)
  FV[(Mb-1,Nb-1)] = s

  # flij - forward last i+j
  # foij - forward old i+j
  # lij is the smallest i+j of the end point of any search path in the current
  # pass.
  # oij is the value of lij of the previous pass. These values will be used
  # to eliminate any entries in P that are no longer nessecary.
  flij = Me+Ne
  foij = Mb+Nb

  finished = (0 != 0)

  # D - the number of non-diagonal steps
  D = 0
  while (D < (Me+Ne) and not(finished)):
    #
    # Work on the forward searches
    #
    D = D + 1
    flij = Me+Ne
    Fkeys = FV.keys()
    Fkeys.sort()
    while (len(Fkeys) > 0):
      key = Fkeys.pop()
      s = FV[key]
      del FV[key]
      (xi,yj) = s.lastpoint
      #
      # now explore the way forward, horizontal first
      #
      keep_horizontal = false
      xih = xi+1
      yjh = yj
      if xih <= Me:
        if FP.has_key((xih,yjh)):
          if FP[(xih,yjh)] < s.lcslength:
            FP[(xih,yjh)] = s.lcslength
            if (yjh == Ne) and (xih == Me):
              #
              # This search path has reached the end point so remove it from
              # the active search paths
              #
              lcs_t = s.snakelist
            else:
              #
              # This search path has not hit the end yet and is still viable
              # so keep it
              #
              keep_horizontal = true
        else:
          FP[(xih,yjh)] = s.lcslength
          if (yjh == Ne) and (xih == Me):
            #
            # This search path has reached the end point so remove it from
            # the active search paths
            #
            lcs_t = s.snakelist
          else:
            #
            # This search path has not hit the end yet and is still viable
            # so keep it
            #
            keep_horizontal = true
        if xih+yjh < flij:
          flij = xih+yjh
      #
      # now explore the vertical direction
      #
      keep_vertical = false
      xiv = xi
      yjv = yj+1
      if yjv <= Ne:
        if FP.has_key((xiv,yjv)):
          if FP[(xiv,yjv)] < s.lcslength:
            FP[(xiv,yjv)] = s.lcslength
            if (yjv == Ne) and (xiv == Me):
              #
              # This search path has reached the end point so remove it from
              # the active search paths
              #
              lcs_t = s.snakelist
            else:
              #
              # This search path has not hit the end yet and is still viable
              # so keep it
              #
              keep_vertical = true
        else:
          FP[(xiv,yjv)] = s.lcslength
          if (yjv == Ne) and (xiv == Me):
            #
            # This search path has reached the end point so remove it from
            # the active search paths
            #
            lcs_t = s.snakelist
          else:
            #
            # This search path has not hit the end yet and is still viable
            # so keep it
            #
            keep_vertical = true
        if xiv+yjv < flij:
          flij = xiv+yjv
      #
      # finally explore the diagonal direction
      #
      keep_diagonal = false
      xid = xi+1
      yjd = yj+1
      if yjd <= Ne and xid <= Me:
        length_diag = s.lcslength \
                    + num_same_token(tol,ref,yjd,dat,xid,feps,ieps)
        if FP.has_key((xid,yjd)):
          if FP[(xid,yjd)] < length_diag:
            keep_diagonal = true
            FP[(xid,yjd)] = length_diag
        else:
          keep_diagonal = true
          FP[(xid,yjd)] = length_diag
        if xid+yjd < flij:
          flij = xid+yjd
      if keep_diagonal:
        if keep_horizontal:
          # Keeping at least both the horizontal and diagonal search directions
          # So create a new search path for the horizontal direction
          sa = copy.copy(s)
          sa.set_lastpoint(xih,yjh)
          FV[(xih,yjh)] = sa
        if keep_vertical:
          # Keeping at least both the vertical and diagonal search directions
          # So create a new search path for the vertical direction
          sa = copy.copy(s)
          sa.set_lastpoint(xiv,yjv)
          FV[(xiv,yjv)] = sa
        # Keeping the diagonal search direction as well
        # So update the search path
        s.set_lastpoint(xid,yjd)
        s.add_snake(xid,yjd,xid,yjd,3)
        s.lcslength = length_diag
        if (yjd == Ne) and (xid == Me):
          #
          # This search path has reached the end point so remove it from
          # the active search paths. The diagonal case is special because
          # the snake needed to be added first.
          #
          lcs_t = s.snakelist
        else:
          FV[(xid,yjd)] = s
      else:
        if keep_horizontal:
          if keep_vertical:
            # Keeping both the horizontal and vertical search directions
            # So create a new search path for the vertical direction
            sa = copy.copy(s)
            sa.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = sa
          # Keeping the horizontal direction in any case
          # So update the search path for the horizontal direction
          s.set_lastpoint(xih,yjh)
          FV[(xih,yjh)] = s
        else:
          if keep_vertical:
            # Keeping only the vertical search direction
            # So update the search path for the vertical direction
            s.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = s
          else:
            # Not keeping anything so let the search path die
            pass
    #
    # now tidy up FP
    #
    ij = foij - (Mb+Nb)
    while (ij <= flij - (Mb+Nb)):
      xi = Mb + ij
      yj = Nb
      while (xi >= Mb):
        if FP.has_key((xi,yj)):
          del FP[(xi,yj)]
        xi = xi - 1
        yj = yj + 1
      ij = ij + 1
    foij = flij
    finished = (len(FV) == 0)
  lcs = transform_snake_list(lcs_t)
  return lcs

def filter_lcs(lcs,Nb,Ne,Mb,Me,add,delete):
  """Take the list of snakes and the add and delete dictionaries and remove
     all tolerable adds and/or deletes from the snake list."""
  lcsout = [ ]
  nsnake = len(lcs)
  if nsnake == 0:
    #
    # Check added lines
    #
    ok_add = false
    if Me-Mb+1 <= 0:
      ok_add = true
      xi1 = Mb
      xi2 = Me
    elif add.has_key(Nb-1):
      (ntoken,lineno,tokenno) = add[Nb-1]
      if ntoken >= Me-Mb+1:
        ok_add = true
        xi1 = Mb
        xi2 = Me
      else:
        xi1 = Mb
        xi2 = Mb
    else:
      xi1 = Mb
      xi2 = Mb
    #
    # Check deleted lines
    #
    ok_del = true
    itoken = Nb
    while (itoken <= Ne):
      ok_del = ok_del and delete.has_key(itoken)
      itoken = itoken + 1
    if ok_del:
      yj1 = Nb
      yj2 = Ne
    else:
      yj1 = Nb
      yj2 = Nb
    if ok_add and ok_del:
      if xi2 < xi1 and yj2 < yj1:
        pass
      else:
        lcsout.append((xi1,yj1,xi2,yj2,2))
  else:
    (xi1,yj1,xi2,yj2,ityp) = lcs.pop(0)
    #
    # Check added lines
    #
    ok_add = false
    if xi1-Mb <= 0:
      xi0 = Mb
      ok_add = true
    elif add.has_key(Nb-1):
      (ntoken,lineno,tokenno) = add[Nb-1]
      if ntoken >= xi1-Mb:
        xi0 = Mb
        ok_add = true
      else:
        xi0 = xi1-1
    else:
      xi0 = xi1-1
    #
    # Check deleted lines
    #
    itoken = Nb
    ok_del = true
    while (itoken < yj1):
      ok_del = ok_del and delete.has_key(itoken)
      itoken = itoken + 1
    if ok_del:
      yj0 = Nb
    else:
      yj0 = yj1-1
    if ok_add and ok_del:
      if xi1 == Mb and yj1 == Nb:
        pass
      else:
        lcsout.append((xi0,yj0,xi1-1,yj1-1,2))
    lcsout.append((xi1,yj1,xi2,yj2,ityp))
    xi3 = xi2
    yj3 = yj2
    while (len(lcs) > 0):
      (xi1,yj1,xi2,yj2,ityp) = lcs.pop(0)
      #
      # Check added lines
      #
      ok_add = false
      if xi1 - xi3 - 1 <= 0:
        xi0 = xi3+1
        ok_add = true
      elif add.has_key(yj3):
        (ntoken,lineno,tokenno) = add[yj3]
        if ntoken >= xi1 - xi3 - 1:
          xi0 = xi3+1
          ok_add = true
        else:
          xi0 = xi1-1
      else:
        xi0 = xi1-1
      #
      # Check deleted lines
      #
      itoken = yj3+1
      ok_del = true
      while (itoken < yj1):
        ok_del = ok_del and delete.has_key(itoken)
        itoken = itoken + 1
      if ok_del:
        yj0 = yj3+1
      else:
        yj0 = yj1-1
      if ok_add and ok_del:
        lcsout.append((xi0,yj0,xi1-1,yj1-1,2))
      lcsout.append((xi1,yj1,xi2,yj2,ityp))
      xi3 = xi2
      yj3 = yj2
    #
    # final bit
    #
    # Check added lines
    #
    ok_add = false
    if Me - xi3 <= 0:
      ok_add = true
      xi4 = Me
    elif add.has_key(yj3):
      (ntoken,lineno,tokenno) = add[yj3]
      if ntoken >= Me - xi3:
        ok_add = true
        xi4 = Me
      else:
        xi4 = xi3+1
    else:
      xi4 = xi3+1
    #
    # Check added lines
    #
    itoken = yj3+1
    ok_del = true
    while (itoken <= Ne):
      ok_del = ok_del and delete.has_key(itoken)
      itoken = itoken + 1
    if ok_del:
      yj4 = Ne
    else:
      yj4 = yj3+1
    if ok_add and ok_del:
      if xi3 == Me and yj3 == Ne:
        pass
      else:
        lcsout.append((xi3+1,yj3+1,xi4,yj4,2))
  return lcsout

def lcs_tokens2lines(lcs,ref,Nrefb,Nrefe,dat,Ndatb,Ndate,nguides):
  """Converts the snake list given in terms of tokens to a snake list in
     terms of lines. This is a two stage process. First the token snake list
     is converted into a list of (pseudo) snakes. For this purpose a new type
     of snake is introduced with type = 4. This type of pseudo snake indicates
     that the lines do not match. Subsequently, the list is parsed to compress
     the snakes so that every line appears only once in a snake or in a gap
     between snakes. In this phase the type 4 pseudo snakes are eliminated
     again, their role was solely to preserve information in the initial
     token to line conversion."""
  #
  # First work out the starting and final line numbers of the both files
  # once and for all.
  #
  list = ref.line2token.keys()
  list.sort()
  nlist = len(list)
  lNrefb = list[0]
  lNrefe = list[nlist-1]
  list = dat.line2token.keys()
  list.sort()
  nlist = len(list)
  if nlist > 0:
    lNdatb = list[0]
    lNdate = list[nlist-1]
  else:
    lNdatb = 0
    lNdate = 0
  #
  # Now convert the snake list
  #
  nsnake = len(lcs)
  line_lcs = []
  if nsnake != 0:
    #
    # First construct a list of pseudo snakes of mis-matching tokens.
    #
    line_tmp_lcs = []
    yj1 = Nrefb
    xi1 = Ndatb
    isnake = 0
    while (isnake < nsnake):
      (xi2,yj2,xi3,yj3,itype) = lcs[isnake]
      xit3 = xi3+1
      yjt3 = yj3+1
      if xi1 < xi2 or yj1 < yj2:
        #
        # There is a mismatch before the next snake
        #
        yjt2  = yj2-1
        xit2  = xi2-1
        line_tmp_lcs.append((xi1,yj1,xit2,yjt2,4))
      xi1 = xi3+1
      yj1 = yj3+1
      isnake = isnake + 1
    #
    # Handle any potential trailing differences
    #
    yj2 = Nrefe
    xi2 = Ndate
    if xi1 < xi2 or yj1 < yj2:
      #
      # There is a mismatch before the end token
      #
      line_tmp_lcs.append((xi1,yj1,xi2,yj2,4))
    line_lcs = line_tmp_lcs
    #
    # Now convert the snake list of mis-matching pseudo snakes from tokens
    # to lines. If a mis-match does not affect an integer number of lines
    # the mis-match cannot be an add or a delete but must be a change.
    #
    line_tmp_lcs = []
    nsnake = len(line_lcs)
    isnake = 0
    while (isnake < nsnake):
      (xi2,yj2,xi3,yj3,itype) = line_lcs[isnake]
      lyj2 = toldiff_tokens.tokenno2lineno(ref,yj2)
      lxi2 = toldiff_tokens.tokenno2lineno(dat,xi2)
      lyj3 = toldiff_tokens.tokenno2lineno(ref,yj3)
      lxi3 = toldiff_tokens.tokenno2lineno(dat,xi3)
      xit2  = xi2-(nguides+1)
      yjt2  = yj2-(nguides+1)
      lyjt2 = toldiff_tokens.tokenno2lineno(ref,yjt2)
      lxit2 = toldiff_tokens.tokenno2lineno(dat,xit2)
      xit3  = xi3+(nguides+1)
      yjt3  = yj3+(nguides+1)
      lyjt3 = toldiff_tokens.tokenno2lineno(ref,yjt3)
      lxit3 = toldiff_tokens.tokenno2lineno(dat,xit3)
      #
      # If the mismatch does not involve exactly a whole set of lines
      # then it cannot be an add or delete but must be a change difference.
      #
      if not(lxit2 < lxi2 and lyjt2 < lyj2 and lxi3 < lxit3 and lyj3 < lyjt3):
        #DEBUG
        #print "lyj2:lyj3<->lxi2:lxi3 = "+str(lyj2)+":"+str(lyj3)+"<->"+str(lxi2)+":"+str(lxi3)
        #DEBUG
        #
        # A correction is needed here to ensure that token adds or deletes
        # just on a line boundary do not result in displaying lines that are
        # actually the same.
        #
        #FIX
        #if lxi2 > lxi3:
        #  lxi3 = max(1,max(lxi2-1,lxi3))
        #if lyj2 > lyj3:
        #  lyj3 = max(1,max(lyj2-1,lyj3))
        #FIX
        # The above is commented out for now. More testing is needed because
        # there still are some spurious differences reported. Also I wanted
        # to test the effect of trimming the snakes. For these tests the rest
        # of the machinery should not change. Nevertheless the local issue
        # needs to be revisited soon.
        lxi3 = max(lxi2,lxi3)
        lyj3 = max(lyj2,lyj3)
      #
      # If the mismatch does not involve anything then that is very strange...
      #
      if (lxi3 < lxi2 and lyj3 < lyj2):
        print "strange: ("+str(lxi2)+","+str(lyj2)+")->("+str(lxi3)+","+str(lyj3)+")"
      #
      # If the mismatch is a change then it might be situated next to or 
      # embedded in white space lines. In that case it should be considered 
      # that the change might actually involve some of the white space lines.
      #
      if lxi2 <= lxi3 and lyj2 <= lyj3:
        lybefore = lyj2-lyjt2-1
        lyafter  = lyjt3-lyj3-1
        lxbefore = lxi2-lxit2-1
        lxafter  = lxit3-lxi3-1
        if lxbefore > lybefore:
          #
          # The data file has more white space lines before the change than
          # the reference file. So some of the data file's white space lines
          # have changed.
          #
          lxi2 = lxi2 - lxbefore + lybefore
          #
        elif lxbefore < lybefore:
          #
          # The reference file has more white space lines before the change than
          # the data file. So some of the reference file's white space lines
          # have changed.
          #
          lyj2 = lyj2 - lybefore + lxbefore
          #
        else:
          #  lxbefore == lybefore
          #
          # The amount of white space in front is the same so no adjustment 
          # is needed.
          #
          pass
          #
        if lxafter > lyafter:
          #
          # The data file has more white space lines after the change than
          # the reference file. So some of the data file's white space lines
          # have changed.
          #
          lxi3 = lxi3 - lyafter + lxafter
          #
        elif lxafter < lyafter:
          #
          # The reference file has more white space lines after the change than
          # the data file. So some of the reference file's white space lines
          # have changed.
          #
          lyj3 = lyj3 - lxafter + lyafter
          #
        else:
          #  lxafter == lyafter
          #
          # The amount of white space in front is the same so no adjustment 
          # is needed.
          #
          pass
          #
      #
      # If the mismatch is a delete then it might be situated next to or
      # embedded in white space lines. In that case an educated guess is needed
      # to pin-point the location in the data file.
      #
      if lxi3 < lxi2:
        lybefore = lyj2-lyjt2-1
        lyafter  = lyjt3-lyj3-1
        lxbefore = lxi2-lxit2-1
        lxafter  = lxit3-lxi2
        if lybefore+lyafter == lxbefore+lxafter:
          #
          # The amount of white space is the same in both the reference and
          # data file, so just reposition the point in the data file to match
          # that in the reference file.
          #
          lxi2 = lxi2 - lxbefore + lybefore
          lxi3 = lxi3 - lxbefore + lybefore
          #
        elif lybefore+lyafter > lxbefore+lxafter:
          #
          # In addition to the tokens some white space lines have been deleted
          # as well. So extend the section of deleted lines as well as 
          # reposition the point in the data file.
          # The GNU diff will first consider white space after the non-white
          # space to have been deleted and then white space before.
          #
          ldel = (lybefore+lyafter)-(lxbefore+lxafter)
          ldelafter  = min(lyafter,ldel)
          ldelbefore = ldel - ldelafter
          lyj2 = lyj2 - ldelbefore
          lyj3 = lyj3 + ldelafter
          lybefore = lybefore - ldelbefore
          lyafter  = lyafter  - ldelafter
          lxi2 = lxi2 - lxbefore + lybefore
          lxi3 = lxi3 - lxbefore + lybefore
          #
        else:
          #  lybefore+lyafter < lxbefore+lxafter
          #
          # Some tokens have been deleted but some white space lines have been
          # added. Together that makes the difference a change difference
          # rather than a delete.
          # The GNU diff will consider any additional white space lines to
          # appear after the changed tokens.
          #
          ladd = (lxbefore+lxafter)-(lybefore+lyafter)
          lxi3 = lxi2 + ladd - 1
          lxi2 = lxi2 - lxbefore + lybefore
          lxi3 = lxi3 - lxbefore + lybefore
      #
      # If the mismatch is an add then it might be situated next to or
      # embedded in white space lines. In that case an educated guess is needed
      # to pin-point the location in the reference file.
      #
      if lyj3 < lyj2:
        lybefore = lyj2-lyjt2-1
        lyafter  = lyjt3-lyj2
        lxbefore = lxi2-lxit2-1
        lxafter  = lxit3-lxi3-1
        if lybefore+lyafter == lxbefore+lxafter:
          #
          # The amount of white space is the same in both the reference and
          # data file, so just reposition the point in the reference file to
          # match that in the data file.
          #
          lyj2 = lyj2 - lybefore + lxbefore
          lyj3 = lyj3 - lybefore + lxbefore
          #
        elif lybefore+lyafter < lxbefore+lxafter:
          #
          # In addition to the tokens some white space lines have been added
          # as well. So extend the section of added lines as well as 
          # reposition the point in the reference file.
          # The GNU diff will first consider white space after the non-white
          # space to have been added and then white space before.
          #
          ladd = (lxbefore+lxafter)-(lybefore+lyafter)
          laddafter  = min(lxafter,ladd)
          laddbefore = ladd - laddafter
          lxi2 = lxi2 - laddbefore
          lxi3 = lxi3 + laddafter
          lxbefore = lxbefore - laddbefore
          lxafter  = lxafter  - laddafter
          lyj2 = lyj2 - lybefore + lxbefore
          lyj3 = lyj3 - lybefore + lxbefore
          #
        else:
          #  lybefore+lyafter > lxbefore+lxafter
          #
          # Some tokens have been added but some white space lines have been
          # deleted. Together that makes the difference a change difference
          # rather than an add.
          # The GNU diff will consider any deleted white space lines to
          # appear after the changed tokens.
          #
          ldel = (lybefore+lyafter)-(lxbefore+lxafter)
          lyj3 = lyj2 + ldel - 1
          lyj2 = lyj2 - lybefore + lxbefore
          lyj3 = lyj3 - lybefore + lxbefore
      #
      # In any case the end point can at worst be one step before the start
      # point (in the case of adds or deletes).
      #
      lxi3 = max(lxi2-1,lxi3)
      lyj3 = max(lyj2-1,lyj3)
      line_tmp_lcs.append((lxi2,lyj2,lxi3,lyj3,itype))
      isnake = isnake + 1
    line_lcs = line_tmp_lcs
    #
    # Finally the snake list of matching lines is constructed from the gaps
    # between the mis-matching lines.
    #
    line_tmp_lcs = []
    nsnake = len(line_lcs)
    isnake = 0
    lxi1 = lNdatb
    lyj1 = lNrefb
    while (isnake < nsnake):
      (lxi2,lyj2,lxi3,lyj3,itype4) = line_lcs[isnake]
      if lxi1 < lxi2 or lyj1 < lyj2:
        line_tmp_lcs.append((lxi1,lyj1,lxi2-1,lyj2-1,1))
      lxi1 = lxi3 + 1
      lyj1 = lyj3 + 1
      isnake = isnake + 1
    lxi2 = toldiff_tokens.tokenno2lineno(dat,Ndate)+1
    lyj2 = toldiff_tokens.tokenno2lineno(ref,Nrefe)+1
    lxi2 = lNdate+1
    lyj2 = lNrefe+1
    if lxi1 < lxi2 or lyj1 < lyj2:
      line_tmp_lcs.append((lxi1,lyj1,lxi2-1,lyj2-1,1))
    #
    # Done, now just return the result
    #
    line_lcs = line_tmp_lcs
  #lNdatb = toldiff_tokens.tokenno2lineno(dat,Ndatb)
  #lNdatb = 1
  #lNdate = toldiff_tokens.tokenno2lineno(dat,Ndate)
  return (line_lcs,lNrefb,lNrefe,lNdatb,lNdate)

def print_lcs(fp,lcs):
  """Print the snake list.
     """
  fp.write("=== printing the snake list ===\n")
  i = 0
  nsnake = len(lcs)
  while (i < nsnake):
    (xi1,yj1,xi2,yj2,ityp) = lcs[i]
    fp.write("(%d:%d)<->(%d:%d) : %d\n" % (yj1,yj2,xi1,xi2,ityp))
    i = i + 1
  fp.write("=== end of the snake list =====\n")
