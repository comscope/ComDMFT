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

import toldiff_lcs
import sys
import copy
import math
import string

def show_tolerance(out_fp,ref_txt,Nb,Ne,chg_txt,add_txt,del_txt,ferr,nguides):
  """Prints out the reference file with the tolerances marked on it.
     This function is essentially an informational aid to help the user
     understand what toldiff is doing to his data.

     As the reference file is available only in a tokenised form it is 
     necessary to provide information per token as to whether this token
     may be deleted or followed by a number of extra tokens. For this 
     purpose if there is a delete or add tolerance on a token it is prepended
     with:
     1. An 'X' if the token may be deleted, '' otherwise
     2. The number of additional tokens if any, '' otherwise
     3. An underscore '_'
     In case some additional tokens may appear in front of the first token this
     will be represented as above but with an empty token string (see
     example 5).
     Examples:
     1. X_token  : indicates the token may not be there
     2. 5_token  : indicates that the token may be followed by 5 additional
                   tokens.
     3. X2_token : indicates the token may be deleted or followed by 2 
                   additional tokens.
     4. token    : indicates that there are no add or delete tolerances on this
                   token.
     5. 5_       : may occur only at the front of the file to indicate that
                   5 additional tokens may be inserted in front of the first
                   token in the reference file.

     For the change tolerances a different approach is needed. First of all it
     has to be recognized that the different types of tokens need different 
     representations of the tolerances:
     1. Strings : Tolerable differences in strings are indicated with the '#'
                  character. To indicate which characters may change the
                  original characters are replaced with '#'.
                  Example: "John" is printed as "J####" if this token always 
                           starts with 'J' and is at most 5 characters long.
     2. Numbers : Tolerable differences in numbers are different from those in
                  strings. So the best thing to do is to use the scientific
                  notation for tolerances. I.e. give the reference value first,
                  followed by '+/-' and the tolerance value. In case of complex
                  numbers the tolerance only applies to the norm of the complex
                  value. So the complex value will be printed as a norm.
                  Examples:
                  1. 123+/-5          : tolerance on integers
                  2. 1.23e2+/-5.0e0   : tolerance on floating point values
                  3. |1.0+2.1j|+/-5.0 : tolerance on complex values

     Finally with respect to guides these should be ignored and only the
     tolerances on the true tokens displayed.

note: +/- is ascii character 177.
        """
  try:
    #
    # Define the +/- character. In the ASCII character set this is number 177 or
    # B1 in hexadecimal
    #
    plusminus = "\xB1"
    underscore = "_"
    #
    # First find out whether there is an add tolerance for token zero. If so
    # a line zero has to be generated.
    #
    if add_txt.has_key(0):
      (ntoken,lnt,tnt) = add_txt[0]
      ntoken = ntoken/(nguides+1)
      line = str(ntoken)+underscore+"\n"
      out_fp.write(line)
    #
    # Now loop over all tokens in the reference file, find out what tolerances
    # are associated with it, and generate the appropriate string.
    #
    iline = 1
    ntoken = len(ref_txt.token)
    itoken = -nguides
    while (itoken < ntoken-nguides):
      itoken = itoken + nguides + 1
      (tr,dr,lnr,tnr) = ref_txt.token[itoken]
      while (iline < lnr):
        out_fp.write("\n")
        iline = iline + 1
      token = ""
      if del_txt.has_key(itoken):
        token = "X"
      if add_txt.has_key(itoken):
        (ntadd,lna,tna) = add_txt[itoken]
        ntadd = ntadd/(nguides+1)
        token = token+str(ntadd)
      if token != "":
        token = token+underscore
      if chg_txt.has_key(itoken):
        (tt,dt,lnt,tnt) = chg_txt[itoken]
        if tt[0] != tr:
          print "huh? tr != tt"
          print "tt = "+str(tt)+" tr = "+str(tr)
        if tr == "s":
          #
          # Express the tolerances on the string
          #
          text = copy.deepcopy(dr)
          text = toldiff_lcs.tol_decode(dt,text)
          token = token + text
          #
        elif tr == "i":
          #
          # Express the tolerances using the number notation
          #
          if tt == "ia":
            token = token + "|"+str(dr)+"|"
          else:
            token = token + str(dr)
          token = token + plusminus + str(dt)
          #
        elif tr == "f":
          #
          # Express the tolerances using the number notation
          #
          if tt == "fa":
            token = token + "|"+str(dr)+"|"
          else:
            token = token + str(dr)
          token = token + plusminus + str(dt)
          #
        elif tr == "c":
          #
          # Express the tolerances using the number notation
          #
          if tt == "ca":
            token = token + "||"+str(dr)+"||"
          else:
            token = token + "|"+str(dr)+"|"
          token = token + plusminus + str(dt)
          #
      else:
        token = token+str(dr)
      token = token+" "
      out_fp.write(token)
    out_fp.write("\n")
  except IOError, e:
    (errno,errmsg) = e
    try:
      ferr.write("toldiff: error writing the marked reference file\n")
      ferr.write("toldiff: error message: ")
      ferr.write(errmsg)
      ferr.write("\n")
    except IOError, e:
      pass
    sys.exit(5)


def show_tolerance_old(out_fp,ref_txt,Nb,Ne,chg_txt,add_txt,del_txt,ferr):
  """Prints out the reference with the tolerance marked on it. This function
     is essentially an informational aid to help the user understand what
     toldiff is doing to his data.

     For every line we print
     1. the number of lines that may be added after it if that exceeds zero
     2. a 'X' if the line may be deleted and a ' ' otherwise
     3. the line of the reference with the characters that may change replaced
        by a '#'"""
  try:
    #
    # First find out how many characters are needed to print the numbers of
    # potentially added lines. Also find out whether the first addition might
    # appear after line zero.
    #
    line_zero = add_txt.has_key(Nb-1)
    maxno = 1
    keys = add_txt.keys()
    while (len(keys) > 0):
      key = keys.pop()
      if add_txt[key] > maxno:
        maxno = add_txt[key]
    numdigit = int(1.0+math.log10(maxno))
    #
    # Use numdigit to create the right format string
    #
    fmtn = "%"+str(numdigit)+"d%c%s"
    fmte = ""
    i = 0
    while (i < numdigit):
      fmte = fmte + " "
      i = i + 1
    fmte = fmte + "%c%s"
    if line_zero:
      nadd = add_txt[Nb-1]
      cdel = " "
      line = "\n"
      out_fp.write(fmtn % (nadd,cdel,line))
    #
    # Loop over all the lines in the reference and sort it out
    #
    i = Nb
    while (i <= Ne):
      if chg_txt.has_key(i):
        line = copy.deepcopy(ref_txt[i])
        line = toldiff_lcs.tol_decode(chg_txt[i],line)
        if string.find(line,"\n") == -1:
          line = line + "\n"
      else:
        line = ref_txt[i]
      if del_txt.has_key(i):
        cdel = "X"
      else:
        cdel = " "
      if add_txt.has_key(i):
        nadd = add_txt[i]
        out_fp.write(fmtn % (nadd,cdel,line))
      else:
        out_fp.write(fmte % (cdel,line))
      i = i + 1
  except IOError, e:
    (errno,errmsg) = e
    try:
      ferr.write("toldiff: error writing the marked reference file\n")
      ferr.write("toldiff: error message: ")
      ferr.write(errmsg)
      ferr.write("\n")
    except IOError, e:
      pass
    sys.exit(5)
