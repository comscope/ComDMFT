# Copyright (C) 2008 Huub van Dam, Science and Technology Facilities Council,
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
import math
import toldiff_lcs

class tokenized_file:
    """
    This class is meant as a container for the data of tokenized reference
    or data files. It has to hold a dictionary of the actual tokens as well
    as a dictionary to translate line numbers to token numbers. I.e. 
    for a given line number it contains the number of the first token on that
    line.
    """

    def __init__(self):
        """Initialisation routine that creates two empty dictionaries."""
        self.token = {}
        self.line2token = {}

def change_tol(list):
    """
    This function takes a list of strings and returns a single change 
    tolerance token. I.e. a token of the form (type,tolerance,line no,
    token no in line).
    """
    type   = list[0]
    tol    = list[1]
    lineno = list[2]
    tokno  = list[3]
    if   type[0] == "t":
        pass
    elif type[0] == "f":
        tol = float(tol)
    elif type[0] == "i":
        tol = int(tol)
    elif type[0] == "c":
        tol = complex(tol)
    else:
        # Something has gone seriously wrong here... Error message to be added.
        pass
    lineno = int(lineno)
    tokno  = int(tokno)
    return (type,tol,lineno,tokno)

def string2token(t,nl,nt):
    """
    This function takes a string and returns a token. A token is a tuple 
    where the first element specifies the type of the data stored in the 
    second element. 

    In this case the data types are limited to numbers, either integer, real
    or complex, and strings. The types a denoted as follows:
    i - integer
    f - float/real
    c - complex
    s - string

    For navigational purposes two more elements added to identify the line
    number (nl) the token was on, and the token number (nt) within the line.
    """
    try:
       i_a = int(t)
       #
       # Toldiff should recognise that -0 and 0 are the same, however, in 
       # a text based comparison that is not automatic so we have to force this.
       #
       if i_a == 0:
         i_a = 0
       token = ("i",i_a,nl,nt)
    except ValueError:
       #
       # In Fortran double precision constants are often printed with a
       # "D" for the exponent rather than an "E", i.e. 1.0E+01 might be
       # printed as 1.0D+01 in Fortran. Python is not aware of this convention
       # so we need to replace any potential "D"-s to obtain valid floating
       # values.
       #
       z = t.replace("d","e")
       z = z.replace("D","e")
       try:
          i_f = float(z)
          #
          # Toldiff should recognise that -0.0 and 0.0 are the same, however,
          # in a text based comparison that is not automatic so we have to
          # force this.
          #
          if i_f == 0.0:
            i_f = 0.0
          token = ("f",i_f,nl,nt)
       except ValueError:
          #
          # The handling of complex numbers is unlikely to work in practice
          # as in most cases complex numbers are printed as (1.0,2.0)
          # rather than 1.0+2.0j. Therefore it is impossible to reliably 
          # distinguish between a complex number and a list of 2 real numbers.
          #
          try:
             i_c = complex(z)
             #
             # Toldiff should recognise that x-0.0j and x+0.0j and that
             # -0.0+y*j and 0.0+y*j are the same, however, in a text based
             # comparison that is not automatic so we have to force this.
             #
             if i_c.real == 0.0:
               i_c = complex(0.0,i_c.imag)
             if i_c.imag == 0.0:
               i_c = complex(i_c.real,0.0)
             token = ("c",i_c,nl,nt)
          except ValueError:
             token = ("s",t,nl,nt)
    return token

def line2strlist(l,separators):
    """
    This routine breaks a line stored in l up and produces a list of 
    strings.
    """
    #separators = ["=","(",")","{","}","[","]",",","*","%",":",";"]
    a = l.split()
    if len(a) == 0:
       #
       # We have found a blank line. Introduce a special token to cope with
       # this. If this token is not introduced it is impossible to say whether
       # blank lines match or not. As a result the diffs would change 
       # significantly.
       #
       # In practice we cannot add a dummy token as it leads to unexpected 
       # results. For example if we delete a line just before a white space
       # line the diff procedure will match the newline token from the deleted
       # line up with the newline token of the remain whitespace line. As a
       # result deleting a single line will appear in the output as a change
       # on 2 lines. Clearly this is very confusing.
       #
       # The alternative of adding dummy tokens on white space lines only
       # turns out to lead to strange results as well. In particular because
       # it is then not clear whether a line with a single token contains 
       # the dummy token or it is just a normal line with a single token on.
       #
       # So ultimately these dummy tokens only cause problems. Therefore they
       # have to be avoided and the whitespace lines have to be dealt with in
       # the token-to-line snake list conversion somehow.
       #
       #a = ["#newline#"]
       pass
    else:
       nseps = len(separators)
       isep  = 0
       while (isep < nseps):
          sep = separators[isep]
          b = []
          while (len(a) > 0):
             tmp = a.pop(0)
             n = tmp.count(sep)
             elm3 = tmp
             while n > 0:
                (elm1,elm2,elm3) = tmp.partition(sep)
                if elm1 != "":
                   b.append(elm1)
                if elm2 != "":
                   b.append(elm2)
                tmp = elm3
                n = n - 1
             if elm3 != "":
                b.append(elm3)
          a = b
          isep = isep + 1
       # Do not do dummy tokens, see above.
       #a.append("#newline#")
    return a

def line2tokens(l,nline,separators):
    """
    This routine takes a line and returns a list of tokens.
    The separators are characters other than whitespace at which
    strings will be split.
    """
    a = line2strlist(l,separators)
    b = []
    ntok = 0
    while (len(a) > 0):
      ntok = ntok + 1
      b.append(string2token(a.pop(0),nline,ntok))
    return b

def compare_tokens(ref,dat,tol,feps,ieps):
    """
    Compare two tokens taking a potential tolerance into account.
    The number returned is the number of characters that are the same.
    """
    tmp = str(tol)
    (tr,dr,lnr,tnr) = ref
    (td,dd,lnd,tnd) = dat
    #result = -2*math.log10(feps)
    result = 0
    if tmp == "":
       if tr == td:
          if tr == "s":
             length = min(len(dr),len(dd))
             i = 0
             while (i < length):
                if dr[i] == dd[i]:
                   result = result + 1
                i = i + 1
          elif (tr == "f") or (tr == "c"):
             denom  = abs(dr)
             #
             # The error enume is divided by 2 to introduce a bonus for matching
             # signs. So that if everything else is equal matching signs will
             # be preferred.
             #
             enume = abs(dr-dd)/2.0
             enuma = abs(abs(dr)-abs(dd))
             enum  = min(enume,enuma)
             if enum <= 0.0:
               inverr = 1.0/feps
             else:
               inverr = denom/enum
               if inverr <= 0.0:
                 inverr = 1.1
             result = result + min(-math.log10(feps),max(math.log10(inverr),0))
          elif tr == "i":
             #
             # The factor 10.0 is there to ensure a non-zero if the reference
             # number is exactly zero.
             #
             # The factor 5.0 is there to ensure that the result is
             # 0 if the difference is 1 order of magnitude larger than the
             # reference value.
             #
             denom  = max(float(abs(dr)),10.0*ieps)*5.0
             #
             # The error enume is divided by 2 to introduce a bonus for matching
             # signs. So that if everything else is equal matching signs will
             # be preferred.
             #
             enume  = max(float(abs(dr-dd)),ieps)/2.0
             enuma  = max(float(abs(abs(dr)-abs(dd))),ieps)
             inverr = max(denom/enume,denom/enuma)
             result = result + max(math.log10(inverr),0)
             #
          else:
             #
             # This must be a guide so if they match they match exactly
             #
             result = result + -math.log10(feps)
       else:
          result = -1
    else:
       (tt,dt,lnt,tnt) = tol
       if tr != tt[0]:
          # the type of the tolerance and the reference token do not match!?
          sys.stdout.write("error mechanism needed here!\n")
       if tr == td:
          if tr == "s":
             tmpref = toldiff_lcs.tol_decode(dt,dr)
             tmpdat = toldiff_lcs.tol_decode(dt,dd)
             length = min(len(tmpref),len(tmpdat))
             i = 0
             while (i < length):
                if tmpref[i] == tmpdat[i]:
                   result = result + 1
                i = i + 1
          elif (tr == "f") or (tr == "c"):
             denom = abs(dr)
             if tt[1] == "a":
                #
                # Compare absolute values
                #
                enum = abs(abs(dr)-abs(dd))
             else:
                #
                # Compare normal values
                #
                # The divide by 2 is to introduce a bonus for matching signs.
                #
                enume = abs(dr-dd)/2.0
                enuma = abs(abs(dr)-abs(dd))
                enum  = min(enume,enuma)
             if enum <= 0.0:
               inverr = 1.0/feps
             else:
               inverr = denom/enum
               if inverr <= 0.0:
                 inverr = 1.1
             result = result + min(-math.log10(feps),max(math.log10(inverr),0))
          elif tr == "i":
             #
             # The factor 10.0 is there to ensure a non-zero if the reference
             # number is exactly zero.
             #
             # The factor 5.0 is there to ensure that the result is
             # 0 if the difference is 1 order of magnitude larger than the
             # reference value.
             #
             denom  = max(float(abs(dr)),10.0*ieps)*5.0
             if tt[1] == "a":
                # compare absolute values
                enume = max(float(abs(abs(dr)-abs(dd))),ieps)
             else:
                # compare normal values
                # the additional term ieps introduces a small penalty for
                # ignoring the sign change so that if everything else is equal
                # the signs will tend to match up
                enume = max(float(abs(dr-dd)),ieps)/2.0
                enuma = max(float(abs(abs(dr)-abs(dd))),ieps)
                enume = min(enume,enuma)
             inverr = denom/enume
             result = result + max(math.log10(inverr),0)
             #
          else:
             #
             # This must be a guide so if they match they match exactly
             #
             result = result + -math.log10(feps)
       else:
          result = -1
    return result

def tokens_match(ref,dat,tol,feps,ieps):
    """
    Compare two tokens taking a potential tolerance into account.
    The value returned is "true" if the tokens match and false otherwise.
    """
    true = (0 == 0)
    false = not(true)
    tmp = str(tol)
    (tr,dr,lnr,tnr) = ref
    (td,dd,lnd,tnd) = dat
    if tmp == "":
       if tr == td:
          if tr == "s":
             length = min(len(dr),len(dd))
             result = len(dr) == len(dd)
             i = 0
             while (i < length):
                result = result and dr[i] == dd[i]
                i = i + 1
          elif (tr == "f") or (tr == "c"):
             #denom  = max(abs(dr),feps)
             #enume  = abs(dr-dd)
             err    = abs(dr-dd)
             result = err <= 0.0
          elif tr == "i":
             err    = abs(dr-dd)
             result = err == 0
          else:
             #
             # This must guide and when the types match the guides must match
             #
             result = true
       else:
          result = false
    else:
       (tt,dt,lnt,tnt) = tol
       if tr != tt[0]:
          # the type of the tolerance and the reference token do not match!?
          sys.stdout.write("error mechanism needed here!\n")
       if tr == td:
          if tr == "s":
             tmpref = toldiff_lcs.tol_decode(dt,dr)
             tmpdat = toldiff_lcs.tol_decode(dt,dd)
             length = min(len(tmpref),len(tmpdat))
             result = len(tmpref) == len(tmpdat)
             i = 0
             while (i < length):
                result = result and tmpref[i] == tmpdat[i]
                i = i + 1
          elif (tr == "f") or (tr == "c"):
             #denom = max(abs(dr),feps)
             if tt[1] == "a":
                # compare absolute values
                #enume = abs(abs(dr)-abs(dd))
                err = abs(abs(dr)-abs(dd))
             else:
                # compare normal values
                #enume = abs(dr-dd)
                err = abs(dr-dd)
             #err   = max(enume/denom,feps)
             result = err <= dt
          elif tr == "i":
             #denom  = max(abs(dr),ieps)
             if tt[1] == "a":
                # compare absolute values
                #enume  = abs(abs(dr)-abs(dd))
                err = abs(abs(dr)-abs(dd))
             else:
                # compare normal values
                #enume  = abs(dr-dd)
                err = abs(dr-dd)
             result =  err <= dt
          else:
             #
             # This must guide and when the types match the guides must match
             #
             result = true
       else:
          result = false
    return result

def tolerance(ref,dat,tol,feps,ieps,itol_scale,ftol_scale,ctol_scale):
    """
    This function generates the tolerance needed to tolerate the difference
    between the reference and the data value, taking any pre-existing 
    tolerances into account.

    The tolerance may be scaled by a scale factor Xtol_scale where X refers
    to the type of tolerance.
    """
    tmp = str(tol)
    (tr,dr,lnr,tnr) = ref
    (td,dd,lnd,tnd) = dat
    result = ""
    #
    # Increase the value for the precision to ensure that tolerances are
    # rounded up to guarantee that accepted values are within the tolerances
    #
    if tr == td:
       if tmp == "":
          if tr == "s":
             nmin = min(len(dr),len(dd))
             nmax = max(len(dr),len(dd))
             i = 0
             tol_ln = ""
             while (i < nmin):
                if dr[i] != dd[i]:
                   tol_ln = tol_ln + "#"
                else:
                   tol_ln = tol_ln + " "
                i = i + 1
             while (i < nmax):
                tol_ln = tol_ln + "#"
                i = i + 1
             dt = toldiff_lcs.tol_encode(tol_ln)
             if dt != "":
                result = ("s",dt,lnr,tnr)
          elif tr == "f":
             #denom = max(abs(dr),feps)
             enumn = abs(dr-dd)*(1.0+10.0*feps)*ftol_scale
             enuma = abs(abs(dr)-abs(dd))*(1.0+10.0*feps)*ftol_scale
             if max(enuma,enumn) > 0.0:
                if enuma < 0.9*enumn:
                   #err = max(enuma/denom*(1.0+feps),feps)
                   result = ("fa",enuma,lnr,tnr)
                else:
                   #err = max(enumn/denom*(1.0+feps),feps)
                   result = ("fd",enumn,lnr,tnr)
          elif tr == "i":
             diffa = int(abs(abs(dr)-abs(dd))*itol_scale+0.5)
             diffn = int(abs(dr-dd)*itol_scale+0.5)
             if max(diffa,diffn) > 0:
                if diffa < 0.9*diffn:
                   result = ("ia",diffa,lnr,tnr)
                else:
                   result = ("id",diffn,lnr,tnr)
          elif tr == "c":
             #denom = max(abs(dr),feps)
             enumn = abs(dr-dd)*(1.0+10.0*feps)*ctol_scale
             enuma = abs(abs(dr)-abs(dd))*(1.0+10.0*feps)*ctol_scale
             if max(enuma,enumn) > 0.0:
                if enuma < 0.9*enumn:
                   #err = max(enuma/denom*(1.0+feps),feps)
                   result = ("ca",enuma,lnr,tnr)
                else:
                   #err = max(enumn/denom*(1.0+feps),feps)
                   result = ("cd",enumn,lnr,tnr)
       else:
          (tt,dt,lnt,tnt) = tol
          if tr == "s":
             nmin = min(len(dr),len(dd))
             nmax = max(len(dr),len(dd))
             i = 0
             tol_ln = toldiff_lcs.tol_decode(dt,"")
             while (i < nmin):
                if dr[i] != dd[i]:
                   tol_ln = tol_ln[:i] + "#" + tol_ln[i+1:]
                i = i + 1
             while (i < nmax):
                tol_ln = tol_ln[:i] + "#" + tol_ln[i+1:]
                tol_ln = tol_ln + "#"
                i = i + 1
             result = ("s",toldiff_lcs.tol_encode(tol_ln),lnt,tnt)
          elif tr == "f":
             #denom = max(abs(dr),feps)
             enumn = abs(dr-dd)*(1.0+10.0*feps)
             enuma = abs(abs(dr)-abs(dd))*(1.0+10.0*feps)
             if enuma < 0.9*enumn or tt == "fa":
                err = enuma
                if err > dt:
                   err = err*ftol_scale
                   result = ("fa",err,lnt,tnt)
                else:
                   result = ("fa",dt,lnt,tnt)
             else:
                err = enumn
                if err > dt:
                   err = err*ftol_scale
                   result = ("fd",err,lnt,tnt)
                else:
                   result = ("fd",dt,lnt,tnt)
          elif tr == "i":
             diffa = abs(abs(dr)-abs(dd))
             diffn = abs(dr-dd)
             if diffa < 0.9*diffn or tt == "ia":
                if diffa > dt:
                   diffa = int(diffa*itol_scale+0.5)
                   result = ("ia",diffa,lnt,tnt)
                else:
                   result = ("ia",dt,lnt,tnt)
             else:
                if diffn > dt:
                   diffn = int(diffn*itol_scale+0.5)
                   result = ("id",diffn,lnt,tnt)
                else:
                   result = ("id",dt,lnt,tnt)
          elif tr == "c":
             #denom = max(abs(dr),feps)
             enumn = abs(dr-dd)*(1.0+10.0*feps)
             enuma = abs(abs(dr)-abs(dd))*(1.0+10.0*feps)
             if enuma < 0.9*enumn or tt == "ca":
                err = enuma
                if err > dt:
                   err = err*ctol_scale
                   result = ("ca",err,lnt,tnt)
                else:
                   result = ("ca",dt,lnt,tnt)
             else:
                err = enumn
                if err > dt:
                   err = err*ctol_scale
                   result = ("cd",err,lnt,tnt)
                else:
                   result = ("cd",dt,lnt,tnt)
    return result

def reconstruct_line(dat,tokno,nguides):
    """
    As all lines have been broken down into tokens it is non-trivial to
    produce the line that holds a particular token. 

    This routine reconstructs as best as possible the line that holds a 
    particular token. Of course this is limited by the fact that all the 
    information about the white space has been lost in the tokenisation.

    In the reconstruction all guides have to suppressed of course, so that
    only the token in the original file are reproduced.

    Returns the reconstructed line and the token number of the first token
    on the next line.
    """
    if len(dat.token) == 0:
        return ("",-1)
    (type,token,lineno,tmp) = dat.token[tokno]
    ntb = lineno2tokenno(dat,lineno)
    nte = lineno2tokenno(dat,lineno+1)-1
    line = ""
    while (ntb <= nte):
        (type,token,linenum,tmp) = dat.token[ntb]
        line = line + " " + str(token)
        ntb = ntb + nguides + 1
    ntb = lineno2tokenno(dat,lineno+1)
    return (line,ntb)

def reconstruct_line_old(dat,tokno):
    """
    As all lines have been broken down into tokens it is non-trivial to
    produce the line that holds a particular token. 

    This routine reconstructs as best as possible the line that holds a 
    particular token. Of course this is limited by the fact that all the 
    information about the white space has been lost in the tokenisation.

    Returns the reconstructed line and the token number of the first token
    on the next line.
    """
    (type,token,lineno,tmp) = dat.token[tokno]
    ntb = dat.line2token[lineno]
    nte = len(dat.token)
    (type,token,linenum,tmp) = dat.token[ntb]
    line = str(token)
    while (lineno == linenum):
        ntb = ntb + 1
        if ntb <= nte:
            (type,token,linenum,tmp) = dat.token[ntb]
        else:
            linenum = 0
        if lineno == linenum:
            line = line + " " + str(token)
    return (line,ntb)


def tokenno2lineno(dat,tokenno):
    """
    This routine takes a token number and returns the corresponding line
    number.

    This routine is needed to produce the diff output as this is represented
    in terms of line numbers whereas the Longest Common Subsequence is given
    in terms of tokens.
    """
    ntoken = len(dat.token)
    if ntoken == 0:
        lineno = 0
    elif tokenno <= 0:
        list = dat.line2token.keys()
        list.sort()
        lineno = list[0]-1
        # lineno = 0
    elif tokenno > ntoken:
        list = dat.line2token.keys()
        list.sort()
        nlist = len(list)
        lineno = list[nlist-1]+1
        # (type,token,lineno,tokno) = dat.token[ntoken]
        # lineno = lineno+1
    else:
        (type,token,lineno,tokno) = dat.token[tokenno]
    return lineno


def lineno2tokenno(dat,lineno):
    """
    This routine takes a line number and returns the token number of the 
    token corresponding to the first token on that line.
    """
    nlines = len(dat.line2token)
    if lineno == 0:
        tokenno = 0
    elif lineno > nlines:
        tokenno = len(dat.token)
        tokenno = tokenno + 1
    else:
        tokenno = dat.line2token[lineno]
    return tokenno
    

def max(a,b):
    """
    Return the maximum of A and B.
    """
    if a > b:
       return a
    else:
       return b

def min(a,b):
    """
    Return the minimum of A and B.
    """
    if a > b:
       return b
    else:
       return a
