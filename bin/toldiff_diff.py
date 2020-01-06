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
import sys
import toldiff_tokens

def lcs_analysis(Nb,Ne,Mb,Me,lcs,identical,equivalent,different):
  """This routine is derived from lcs_to_diff. Instead of writing out the
     diff based on the snake list it analyses the snake list and establishes
     whether the reference file and the data file are identical, equivalent,
     or different. 

     In this case the three results mean:
     - identical: there are no differences between the two files at all.
     - equivalent: there are only tolerated differences.
     - different : there are some non-tolerated differences.

     The values for the results are taken from the argument list.
     """
  analysis = identical
  xi1 = Mb-1
  yj1 = Nb-1
  mxtype = 0
  Nsnake = len(lcs)
  Isnake = 0
  itype  = 0
  if Nsnake == 0:
    if Nb <= Ne:
      analysis = different
    else:
      if Mb <= Me:
        analysis = different
      else:
        pass
  else:
    while (Isnake < Nsnake):
      (xi2,yj2,xi3,yj3,itype) = lcs[Isnake]
      Isnake = Isnake + 1
      if itype > mxtype:
        mxtype = itype
    if mxtype == 1:
      # there are only exact matches so identical still is the best possible
      pass
    elif mxtype == 2:
      # there are tolerated differences so equivalent is the best possible
      analysis = equivalent
    elif mxtype == 3:
      # there are non-tolerated differences so different is the best possible
      analysis = different
      return analysis
    Isnake = -1
    while (Isnake < Nsnake):
      Isnake = Isnake + 1
      if (Isnake < Nsnake):
        (xi2,yj2,xi3,yj3,itype) = lcs[Isnake]
      else:
        xi2 = Me+1
        yj2 = Ne+1
        xi3 = Me+1
        yj3 = Ne+1
      if xi1+1 <= xi2 and yj1+1 <= yj2:
        if xi1+1 == xi2:
          if yj1+1 == yj2:
            #
            # This is a continuation of the previous snake (of a different type)
            #
            pass
          else:
            analysis = different
        else:
          analysis = different
      xi1 = xi3
      yj1 = yj3
  return analysis


def reconstruct_line(dat,lineno,nguides):
  """This routine takes data representing a file (reference or data file)
     and a line number. It then reconstructs the requested line from the 
     token data and returns both the reconstructed line and the number of
     the next line."""
  tokenno = toldiff_tokens.lineno2tokenno(dat,lineno)
  linet = toldiff_tokens.tokenno2lineno(dat,tokenno)
  if linet == lineno:
    #
    # This is a non-empty line, so reconstruct it.
    #
    (line,tokenno) = toldiff_tokens.reconstruct_line(dat,tokenno,nguides)
    #
  else:
    #
    # This is an empty line, so just return that.
    #
    line = ""
    #
  return (line,lineno+1)

def lcs_to_diff(ref,Nb,Ne,dat,Mb,Me,lcs,fout,ferr,nguides):
  """This routine takes the reference file and the data file as well as the
     longest common subsequence stored as a snake list. Based on this
     information it writes the diff file.

     The exact form of the diff file depends on the types of snake present
     in the snake list as well. If there are type 3 snakes then lines have 
     been matched up as best as one possibly can. In this case only adds,
     deletes and changes with equal numbers of reference and data lines are
     possible.

     Otherwise the snake list is interpreted in the classical diff way.

     Although the snake list is presented in terms of lines the reference
     and data files are given in terms of tokens. So whenever a line of 
     either file has to be printed this line first has to be reconstructed
     from the corresponding tokens. The routine reconstruct_line
     performs this function given the reference or data file information as well
     as a line number (this is not to be confused with the routine
     toldiff_tokens.reconstruct_line which takes a token number).
     """
  try:
    xi1 = Mb-1
    yj1 = Nb-1
    mxtype = 0
    Nsnake = len(lcs)
    Isnake = 0
    itype  = 0
    if Nsnake == 0:
      if Nb == Ne:
        if Mb == Me:
          fout.write("%dc%d\n" % (Nb,Mb))
          (line,iline) = reconstruct_line(ref,Nb,nguides)
          fout.write("< "+line+"\n")
          fout.write("---\n")
          (line,iline) = reconstruct_line(dat,Mb,nguides)
          fout.write("> "+line+"\n")
        elif Mb < Me:
          fout.write("%dc%d,%d\n" % (Nb,Mb,Me))
          (line,iline) = reconstruct_line(ref,Nb,nguides)
          fout.write("< "+line+"\n")
          fout.write("---\n")
          iline = Mb
          while (iline <= Me):
            (line,iline) = reconstruct_line(dat,iline,nguides)
            fout.write("> "+line+"\n")
        else:
          fout.write("%dd%d\n" % (Nb,Mb-1))
          (line,iline) = reconstruct_line(ref,Nb,nguides)
          fout.write("< "+line+"\n")
      elif Nb < Ne:
        if Mb == Me:
          fout.write("%d,%dc%d\n" % (Nb,Ne,Mb))
          iline = Nb
          while (iline <= Ne):
            (line,iline) = reconstruct_line(ref,iline,nguides)
            fout.write("< "+line+"\n")
          fout.write("---\n")
          (line,iline) = reconstruct_line(dat,Mb,nguides)
          fout.write("> "+line+"\n")
        elif Mb < Me:
          fout.write("%d,%dc%d,%d\n" % (Nb,Ne,Mb,Me))
          iline = Nb
          while (iline <= Ne):
            (line,iline) = reconstruct_line(ref,iline,nguides)
            fout.write("< "+line+"\n")
          fout.write("---\n")
          iline = Mb
          while (iline <= Me):
            (line,iline) = reconstruct_line(dat,iline,nguides)
            fout.write("> "+line+"\n")
        else:
          fout.write("%d,%dd%d\n" % (Nb,Ne,Mb-1))
          iline = Nb
          while (iline <= Ne):
            (line,iline) = reconstruct_line(ref,iline,nguides)
            fout.write("< "+line+"\n")
      else:
        if Mb == Me:
          fout.write("%da%d\n" % (Nb-1,Mb))
          (line,iline) = reconstruct_line(dat,Mb,nguides)
          fout.write("> "+line+"\n")
        elif Mb < Me:
          fout.write("%da%d,%d\n" % (Nb-1,Mb,Me))
          iline = Mb
          while (iline <= Me):
            (line,iline) = reconstruct_line(dat,iline,nguides)
            fout.write("> "+line+"\n")
        else:
          pass
    else:
      while (Isnake < Nsnake):
        (xi2,yj2,xi3,yj3,itype) = lcs[Isnake]
        Isnake = Isnake + 1
        if itype > mxtype:
          mxtype = itype
      Isnake = -1
      while (Isnake < Nsnake):
        Isnake = Isnake + 1
        if (Isnake < Nsnake):
          (xi2,yj2,xi3,yj3,itype) = lcs[Isnake]
        else:
          xi2 = Me+1
          yj2 = Ne+1
          xi3 = Me+1
          yj3 = Ne+1
        if xi1+1 <= xi2 and yj1+1 <= yj2:
          if xi1+1 == xi2:
            if yj1+1 == yj2:
              #
              # This is a continuation of the previous snake (of a different type)
              #
              pass
            else:
              if yj1+1 == yj2-1:
                fout.write("%dd%d\n" % (yj1+1,xi1))
              else:
                fout.write("%d,%dd%d\n" % (yj1+1,yj2-1,xi1))
              iline = yj1+1
              while (iline <= yj2-1):
                (line,iline) = reconstruct_line(ref,iline,nguides)
                fout.write("< "+line+"\n")
          else:
            if yj1+1 == yj2:
              if xi1+1 == xi2-1:
                fout.write("%da%d\n" % (yj1,xi1+1))
              else:
                fout.write("%da%d,%d\n" % (yj1,xi1+1,xi2-1))
              iline = xi1+1
              while (iline <= xi2-1):
                (line,iline) = reconstruct_line(dat,iline,nguides)
                fout.write("> "+line+"\n")
            else:
              if xi1+1 == xi2-1:
                if yj1+1 == yj2-1:
                  if mxtype < 3:
                    # use the classical interpretation
                    fout.write("%dc%d\n" % (yj1+1,xi1+1))
                    (line,iline) = reconstruct_line(ref,yj1+1,nguides)
                    fout.write("< "+line+"\n")
                    fout.write("---\n")
                    (line,iline) = reconstruct_line(dat,xi1+1,nguides)
                    fout.write("> "+line+"\n")
                  else:
                    # these lines were not matched with any others in a best fit
                    fout.write("%dd%d\n" % (yj1+1,xi1+1))
                    (line,iline) = reconstruct_line(ref,yj1+1,nguides)
                    fout.write("< "+line+"\n")
                    fout.write("%da%d\n" % (yj1+1,xi1+1))
                    (line,iline) = reconstruct_line(dat,xi1+1,nguides)
                    fout.write("> "+line+"\n")
                else:
                  if mxtype < 3:
                    # use the classical interpretation
                    fout.write("%d,%dc%d\n" % (yj1+1,yj2-1,xi1+1))
                    iline = yj1+1
                    while (iline <= yj2-1):
                      (line,iline) = reconstruct_line(ref,iline,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("---\n")
                    iline = xi1+1
                    while (iline <= xi2-1):
                      (line,iline) = reconstruct_line(dat,iline,nguides)
                      fout.write("> "+line+"\n")
                  else:
                    # this line was not matched with any other in a best fit
                    fout.write("%d,%dd%d\n" % (yj1+1,yj2-1,xi1+1))
                    iline = yj1+1
                    while (iline <= yj2-1):
                      (line,iline) = reconstruct_line(ref,iline,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("%da%d,%d\n" % (yj1+1,xi1+1,xi2-1))
                    iline = xi1+1
                    while (iline <= xi2-1):
                      (line,iline) = reconstruct_line(dat,iline,nguides)
                      fout.write("> "+line+"\n")
              else:
                if yj1+1 == yj2-1:
                  if mxtype < 3:
                    # use the classical interpretation
                    fout.write("%dc%d,%d\n" % (yj1+1,xi1+1,xi2-1))
                    iline = yj1+1
                    while (iline <= yj2-1):
                      (line,iline) = reconstruct_line(ref,iline,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("---\n")
                    iline = xi1+1
                    while (iline <= xi2-1):
                      (line,iline) = reconstruct_line(dat,iline,nguides)
                      fout.write("> "+line+"\n")
                  else:
                    # this line was not matched with any other in a best fit
                    fout.write("%d,%dd%d\n" % (yj1+1,yj2-1,xi1+1))
                    iline = yj1+1
                    while (iline <= yj2-1):
                      (line,iline) = reconstruct_line(ref,iline,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("%da%d,%d\n" % (yj1+1,xi1+1,xi2-1))
                    iline = xi1+1
                    while (iline <= xi2-1):
                      (line,iline) = reconstruct_line(dat,iline,nguides)
                      fout.write("> "+line+"\n")
                else:
                  if mxtype < 3:
                    # use the classical interpretation
                    fout.write("%d,%dc%d,%d\n" % (yj1+1,yj2-1,xi1+1,xi2-1))
                    iline = yj1+1
                    while (iline <= yj2-1):
                      (line,iline) = reconstruct_line(ref,iline,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("---\n")
                    iline = xi1+1
                    while (iline <= xi2-1):
                      (line,iline) = reconstruct_line(dat,iline,nguides)
                      fout.write("> "+line+"\n")
                  else:
                    # these lines were not matched with any others in a best fit
                    fout.write("%d,%dd%d\n" % (yj1+1,yj2-1,xi1+1))
                    iline = yj1+1
                    while (iline <= yj2-1):
                      (line,iline) = reconstruct_line(ref,iline,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("%da%d,%d\n" % (yj1+1,xi1+1,xi2-1))
                    iline = xi1+1
                    while (iline <= xi2-1):
                      (line,iline) = reconstruct_line(dat,iline,nguides)
                      fout.write("> "+line+"\n")
        #
        if itype == 3:
          #
          # Write out a type 3 snake (best fit, but no matches)
          #
          if xi2 == xi3:
            fout.write("%dc%d\n" % (yj2,xi2))
          else:
            fout.write("%d,%dc%d,%d\n" % (yj2,yj3,xi2,xi3))
          iline = yj2
          while (iline <= yj3):
            (line,iline) = reconstruct_line(ref,iline,nguides)
            fout.write("< "+line+"\n")
          fout.write("---\n")
          iline = xi2
          while (iline <= xi3):
            (line,iline) = reconstruct_line(dat,iline,nguides)
            fout.write("> "+line+"\n")
        xi1 = xi3
        yj1 = yj3
  except IOError, e:
    (errno,errmsg) = e
    try:
      ferr.write("toldiff: error writing diff output\n")
      ferr.write("toldiff: error message: ")
      ferr.write(errmsg)
      ferr.write("\n")
    except IOError, e:
      pass
    sys.exit(5)

def lcs_to_diff_tokens_broken(ref,Nb,Ne,dat,Mb,Me,lcs,fout,ferr):
  # This routine is broken as the conversion from tokens to lines is more
  # complicated than simply translating the numbers. Some additional filtering
  # is needed. The proper routine is lcs_to_diff which still drives of a snake
  # list in terms of line numbers. The token to line conversion is performed
  # by a separate routine.
  """This routine takes the reference file and the data file as well as the
     longest common subsequence stored as a snake list. Based on this
     information it writes the diff file.

     This procedure is slightly complicated by the fact that the snake list
     is given in terms of tokens whereas diff presents its output in terms
     of lines. There it is required to convert everything back to lines.

     The exact form of the diff file depends on the types of snake present
     in the snake list as well. If there are type 3 snakes then lines have 
     been matched up as best as one possibly can. In this case only adds,
     deletes and changes with equal numbers of reference and data lines
     will be generated.

     Otherwise the snake list is interpreted in the classical diff way.
     """
  try:
    #
    # Comment on notation the variable names starting with "l" refer to 
    # line numbers, other variables refer to tokens.
    #
    xi1 = Mb-1
    yj1 = Nb-1
    lxi1 = toldiff_tokens.tokenno2lineno(dat,Mb-1)
    lyj1 = toldiff_tokens.tokenno2lineno(ref,Nb-1)
    mxtype = 0
    Nsnake = len(lcs)
    Isnake = 0
    itype  = 0
    if Nsnake == 0:
      lNb = toldiff_tokens.tokenno2lineno(ref,Nb)
      lNe = toldiff_tokens.tokenno2lineno(ref,Ne)
      lMb = toldiff_tokens.tokenno2lineno(dat,Mb)
      lMe = toldiff_tokens.tokenno2lineno(dat,Me)
      if lNb == lNe:
        if lMb == lMe:
          fout.write("%dc%d\n" % (lNb,lMb))
          (line,itoken) = toldiff_tokens.reconstruct_line(ref,Nb,nguides)
          fout.write("< "+line+"\n")
          fout.write("---\n")
          (line,itoken) = toldiff_tokens.reconstruct_line(dat,Mb,nguides)
          fout.write("> "+line+"\n")
        elif lMb < lMe:
          fout.write("%dc%d,%d\n" % (lNb,lMb,lMe))
          (line,itoken) = toldiff_tokens.reconstruct_line(ref,Nb,nguides)
          fout.write("< "+line+"\n")
          fout.write("---\n")
          itoken = Mb
          while (itoken <= Me):
            (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
            fout.write("> "+line+"\n")
        else:
          fout.write("%dd%d\n" % (lNb,lMb-1))
          (line,itoken) = toldiff_tokens.reconstruct_line(ref,Nb,nguides)
          fout.write("< "+line+"\n")
      elif lNb < lNe:
        if lMb == lMe:
          fout.write("%d,%dc%d\n" % (lNb,lNe,lMb))
          itoken = Nb
          while (itoken <= Ne):
            (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
            fout.write("< "+line+"\n")
          fout.write("---\n")
          (line,itoken) = toldiff_tokens.reconstruct_line(dat,Mb,nguides)
          fout.write("> "+line+"\n")
        elif lMb < lMe:
          fout.write("%d,%dc%d,%d\n" % (lNb,lNe,lMb,lMe))
          itoken = Nb
          while (itoken <= Ne):
            (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
            fout.write("< "+line+"\n")
          fout.write("---\n")
          itoken = Mb
          while (itoken <= Me):
            (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
            fout.write("> "+line+"\n")
        else:
          fout.write("%d,%dd%d\n" % (lNb,lNe,lMb-1))
          itoken = Nb
          while (itoken <= Ne):
            (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
            fout.write("< "+line+"\n")
      else:
        if lMb == lMe:
          fout.write("%da%d\n" % (lNb-1,lMb))
          (line,itoken) = toldiff_tokens.reconstruct_line(dat,Mb,nguides)
          fout.write("> "+line+"\n")
        elif lMb < lMe:
          fout.write("%da%d,%d\n" % (lNb-1,lMb,lMe))
          itoken = Mb
          while (itoken <= Me):
            (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
            fout.write("> "+line+"\n")
        else:
          pass
    else:
      while (Isnake < Nsnake):
        (xi2,yj2,xi3,yj3,itype) = lcs[Isnake]
        Isnake = Isnake + 1
        if itype > mxtype:
          mxtype = itype
      Isnake = -1
      while (Isnake < Nsnake):
        Isnake = Isnake + 1
        if (Isnake < Nsnake):
          (xi2,yj2,xi3,yj3,itype) = lcs[Isnake]
        else:
          xi2 = Me+1
          yj2 = Ne+1
          xi3 = Me+1
          yj3 = Ne+1
        lxi2 = toldiff_tokens.tokenno2lineno(dat,xi2)
        lxi3 = toldiff_tokens.tokenno2lineno(dat,xi3)
        lyj2 = toldiff_tokens.tokenno2lineno(ref,yj2)
        lyj3 = toldiff_tokens.tokenno2lineno(ref,yj3)
        if lxi1+1 <= lxi2 and lyj1+1 <= lyj2:
          if lxi1+1 == lxi2:
            if lyj1+1 == lyj2:
              #
              # This is a continuation of the previous snake (of a different type)
              #
              pass
            else:
              if lyj1+1 == lyj2-1:
                fout.write("%dd%d\n" % (lyj1+1,lxi1))
              else:
                fout.write("%d,%dd%d\n" % (lyj1+1,lyj2-1,lxi1))
              itoken = yj1+1
              while (itoken <= yj2-1):
                (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
                fout.write("< "+line+"\n")
          else:
            if lyj1+1 == lyj2:
              if lxi1+1 == lxi2-1:
                fout.write("%da%d\n" % (lyj1,lxi1+1))
              else:
                fout.write("%da%d,%d\n" % (lyj1,lxi1+1,lxi2-1))
              itoken = xi1+1
              while (itoken <= xi2-1):
                (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
                fout.write("> "+line+"\n")
            else:
              if lxi1+1 == lxi2-1:
                if lyj1+1 == lyj2-1:
                  if mxtype < 3:
                    # use the classical interpretation
                    fout.write("%dc%d\n" % (lyj1+1,lxi1+1))
                    (line,itoken) = toldiff_tokens.reconstruct_line(ref,yj1+1,nguides)
                    fout.write("< "+line+"\n")
                    fout.write("---\n")
                    (line,itoken) = toldiff_tokens.reconstruct_line(dat,xi1+1,nguides)
                    fout.write("> "+line+"\n")
                  else:
                    # these lines were not matched with any others in a best fit
                    fout.write("%dd%d\n" % (lyj1+1,lxi1+1))
                    (line,itoken) = toldiff_tokens.reconstruct_line(ref,yj1+1,nguides)
                    fout.write("< "+line+"\n")
                    fout.write("%da%d\n" % (lyj1+1,lxi1+1))
                    (line,itoken) = toldiff_tokens.reconstruct_line(dat,xi1+1,nguides)
                    fout.write("> "+line+"\n")
                else:
                  if mxtype < 3:
                    # use the classical interpretation
                    fout.write("%d,%dc%d\n" % (lyj1+1,lyj2-1,lxi1+1))
                    itoken = yj1+1
                    while (itoken <= yj2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("---\n")
                    itoken = xi1+1
                    while (itoken <= xi2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
                      fout.write("> "+line+"\n")
                  else:
                    # this line was not matched with any other in a best fit
                    fout.write("%d,%dd%d\n" % (lyj1+1,lyj2-1,lxi1+1))
                    itoken = yj1+1
                    while (itoken <= yj2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("%da%d,%d\n" % (lyj1+1,lxi1+1,lxi2-1))
                    itoken = xi1+1
                    while (itoken <= xi2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
                      fout.write("> "+line+"\n")
              else:
                if lyj1+1 == lyj2-1:
                  if mxtype < 3:
                    # use the classical interpretation
                    fout.write("%dc%d,%d\n" % (lyj1+1,lxi1+1,xi2-1))
                    itoken = yj1+1
                    while (itoken <= yj2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("---\n")
                    itoken = xi1+1
                    while (itoken <= xi2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
                      fout.write("> "+line+"\n")
                  else:
                    # this line was not matched with any other in a best fit
                    fout.write("%d,%dd%d\n" % (lyj1+1,lyj2-1,lxi1+1))
                    itoken = yj1+1
                    while (itoken <= yj2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("%da%d,%d\n" % (lyj1+1,lxi1+1,lxi2-1))
                    itoken = xi1+1
                    while (itoken <= xi2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
                      fout.write("> "+line+"\n")
                else:
                  if mxtype < 3:
                    # use the classical interpretation
                    fout.write("%d,%dc%d,%d\n" % (lyj1+1,lyj2-1,lxi1+1,lxi2-1))
                    itoken = yj1+1
                    while (itoken <= yj2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("---\n")
                    itoken = xi1+1
                    while (itoken <= xi2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
                      fout.write("> "+line+"\n")
                  else:
                    # these lines were not matched with any others in a best fit
                    fout.write("%d,%dd%d\n" % (lyj1+1,lyj2-1,lxi1+1))
                    itoken = yj1+1
                    while (itoken <= yj2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
                      fout.write("< "+line+"\n")
                    fout.write("%da%d,%d\n" % (lyj1+1,lxi1+1,lxi2-1))
                    itoken = xi1+1
                    while (itoken <= xi2-1):
                      (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
                      fout.write("> "+line+"\n")
        #
        if itype == 3:
          #
          # Write out a type 3 snake (best fit, but no matches)
          #
          if lxi2 == lxi3:
            fout.write("%dc%d\n" % (lyj2,lxi2))
          else:
            fout.write("%d,%dc%d,%d\n" % (lyj2,lyj3,lxi2,lxi3))
          itoken = yj2
          while (itoken <= yj3):
            (line,itoken) = toldiff_tokens.reconstruct_line(ref,itoken,nguides)
            fout.write("< "+line+"\n")
          fout.write("---\n")
          itoken = xi2
          while (itoken <= xi3):
            (line,itoken) = toldiff_tokens.reconstruct_line(dat,itoken,nguides)
            fout.write("> "+line+"\n")
        xi1  = xi3
        yj1  = yj3
        lxi1 = lxi3
        lyj1 = lyj3
  except IOError, e:
    (errno,errmsg) = e
    try:
      ferr.write("toldiff: error writing diff output\n")
      ferr.write("toldiff: error message: ")
      ferr.write(errmsg)
      ferr.write("\n")
    except IOError, e:
      pass
    sys.exit(5)

def diff_to_lcs(Nb,Ne,Mb,Me,fin,ferr):
  """This routine reads the difference file resulting from diff-ing the
     tokenized reference and data files. It extracts the token numbers that 
     differ and reconstructs the Longest Common Subsequence (the actual tokens
     do not need to be extracted as they are already available in memory from
     the tokenisation of the reference and data files).

     It returns lcs where:
     - lcs: the list of snakes defining the longest common subsequence
     """
  try:
    xi1 =  1
    yj1 =  1
    lcs = [ ]
    line = fin.readline()
    while (line):
      #
      # Analyse the type of the difference as well as the token numbers
      # This analysis is slightly complicated by the clever way diff handles
      # the line numbers. Different valid forms are:
      #   323d525
      #   324d525,632
      #   324,456c525,623
      #   456a345
      # etc. So to discover what is there we need to:
      # 1. Find the character on the line
      # 2. Split the line at the character
      # 3. Split the two line number sections to learn the token numbers for 
      #    the reference file and data file.
      #
      if string.find(line,"a") > -1:
        type = "a"
      elif string.find(line,"d") > -1:
        type = "d"
      elif string.find(line,"c") > -1:
        type = "c"
      else:
        try:
          ferr.write("toldiff: diff_to_lcs: cannot interpret line: %s\n" % line)
        except IOError, e:
          pass
        exit(2)
      (refstr,datstr) = string.split(line,type)
      ref_list = string.split(refstr,",")
      dat_list = string.split(datstr,",")
      #
      # Sort out which token ranges to read for the reference file
      #
      ref_begin = int(ref_list[0])
      if len(ref_list) == 2:
        ref_end = int(ref_list[1])
      else: 
        if type == "a":
          #
          # The reference file difference part is empty
          #
          ref_end = ref_begin - 1
        else:
          ref_end = ref_begin
      #
      # Sort out which token ranges to read for the data file
      #
      dat_begin = int(dat_list[0])
      if len(dat_list) == 2:
        dat_end = int(dat_list[1])
      else: 
        if type == "d":
          #
          # The data file difference part is empty
          #
          dat_end = dat_begin - 1
        else:
          dat_end = dat_begin
      #
      # Read the reference file lines
      #
      if (type == "d") or (type == "c"):
        n = ref_begin
        while (n <= ref_end):
          line = fin.readline()
          n = n + 1
      #
      # Skip the separator line
      #
      if type == "c":
        fin.readline()
      #
      # Read the data file lines
      #
      if (type == "a") or (type == "c"):
        n = dat_begin
        while (n <= dat_end):
          line = fin.readline()
          n = n + 1
      #
      # Add to the LCS
      #
      # The starting point of the snake was sorted in the previous pass, so
      # here we need to sort out the end point of the snake.
      #
      if type == "d":
        #
        # The data file snake end point is the line mentioned in the diff.
        #
        xi2 = dat_begin
        #
      else:
        #
        # The data file snake end point is the line before the difference
        # begins.
        #
        xi2 = dat_begin - 1
        #
      if type == "a":
        #
        # The reference file snake end point is the line mentioned in the
        # diff.
        #
        yj2 = ref_begin
        #
      else:
        #
        # The reference file snake end point is the line before the difference
        # begins.
        #
        yj2 = ref_begin - 1
      if (xi2 >= xi1) and (yj2 >= yj1):
        # we have a snake of length >= 1 (this allows for empty snakes)
        lcs.append((xi1,yj1,xi2,yj2,1))
      #
      # Now address the beginning of the next snake.
      #
      if type == "d":
        #
        # The next snake will start on the line after the one mentioned in 
        # data file difference part.
        #
        xi1 = dat_begin + 1
        #
      else:
        #
        # The next snake will start on the line after the end of the data file
        # difference part.
        #
        xi1 = dat_end + 1
        #
      if type == "a":
        #
        # The next snake will start on the line after the one mentioned in 
        # reference file difference part.
        #
        yj1 = ref_begin + 1
        #
      else:
        #
        # The next snake will start on the line after the end of the reference
        # file difference part.
        #
        yj1 = ref_end + 1
        #
      line = fin.readline()
    #
    # Add the final snake if any
    #
    xi2 = Me
    yj2 = Ne
    if (xi2 >= xi1) and (yj2 >= yj1):
      # we have a snake of length >= 1 (this allows for empty snakes)
      lcs.append((xi1,yj1,xi2,yj2,1))
  except IOError, e:
    (errno,errmsg) = e
    try:
      ferr.write("toldiff: error reading difference file from external diff\n")
      ferr.write("toldiff: error message: ")
      ferr.write(errmsg)
      ferr.write("\n")
    except IOError, e:
      pass
    sys.exit(5)
  return lcs

