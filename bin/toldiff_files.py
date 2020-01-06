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
import string
import toldiff_tokens

def max(a,b):
  """Return the maximum of a and b."""
  if a >= b:
    r = a
  else:
    r = b
  return r

def mod(a,b):
  """Return the modulus of a with respect to b."""
  c = int(a/b)
  result = a - b*c
  return result

# now follows a random string that can be changed to force subversion to
# update the version number of this file (sigh...): aaa

def version_toldiff(fp,errfp):
  """Print out the version information to the specified file object."""
  revision = "$Revision: 75 $"
  date     = "$Date: 2008-07-28 18:21:53 -0400 (Mon, 28 Jul 2008) $"
  rlist    = string.split(revision)
  dlist    = string.split(date)
  version  = "toldiff version "+rlist[1]+" of "+dlist[1]+"\n"
  try:
    fp.write(version)
  except IOError, e:
    (errno,errmsg) = e
    try:
      errfp.write("toldiff: error writing version information\n")
      errfp.write("toldiff: error message: ")
      errfp.write(errmsg)
      errfp.write("\n")
    except IOError, e:
      pass
    sys.exit(5)

def load_plain_text(fp,data,nline,ntoken,separators,nguides):
  """Reads a plain text file refered to by the file descriptor 'fp'. The 
     contents of the file are stored in the dictionary 'data'. The number of
     lines read are returned in 'nline'. The number of tokens read are
     returned in 'ntoken'.
     The token number is used as the key in the dictionary, the lines themselves
     are tokenized and the tokens are added. Finally the first real line is
     stored as token number 1.
     With respect to the guides when the tokens are stored an additional 
     nguides guide tokens are introduced. These tokens are of the type string
     and contain a string depending on the type of the token. The last token
     on the line has an a modified guide token to go along with it. The guides
     follow the true token.
     """
  #
  # IOError exceptions are handled by the calling routine.
  #
  efftoken = 0
  line = fp.readline()
  while line:
    nline = nline + 1
    data.line2token[nline] = (nguides+1)*(ntoken+1)-nguides
    toklist = toldiff_tokens.line2tokens(line,nline,separators)
    if len(toklist) > 0:
      while len(toklist) > 0:
        ntoken = ntoken + 1
        efftoken = (nguides+1)*ntoken-nguides
        token = toklist.pop(0)
        data.token[efftoken] = token
        #
        # If needed insert the guides
        #
        if nguides > 0:
          (type,value,line,tokenno) = token
          if type == "s":
            guide = ("sg","_s_",line,tokenno)
          elif type == "f":
            guide = ("fg",123123.123123,line,tokenno)
          elif type == "i":
            guide = ("ig",123123,line,tokenno)
          elif type == "c":
            guide = ("cg",123123.123123+123123.123123j,line,tokenno)
          efftoken = efftoken + 1
          while (efftoken <= (nguides+1)*ntoken):
            data.token[efftoken] = guide
            efftoken = efftoken + 1
      #
      # For the last token on the line modify the guide to incorporate
      # the newline character.
      #
      if nguides > 0:
        if type == "s":
          guide = ("sn","_n_",line,tokenno)
        elif type == "f":
          guide = ("fn",369369.369369,line,tokenno)
        elif type == "i":
          guide = ("in",369369,line,tokenno)
        elif type == "c":
          guide = ("cn",369369.369369+369369.369369j,line,tokenno)
        efftoken = efftoken - nguides
        while (efftoken <= (nguides+1)*ntoken):
          data.token[efftoken] = guide
          efftoken = efftoken + 1
        efftoken = efftoken - 1
    line = fp.readline()
  return (data,nline,efftoken)

def load_tolerances(fp,separators,nguides):
  """Reads a tolerance file from the file refered to by the file descriptor
     'fp'. The contents of the file are stored in three sections:

     Section 1: The change section
     -----------------------------

     Contains a specification of tokens that are allowed to differ and how
     they are allowed to differ. Per line a single token is stored as:

     "<token number> <token type> <tolerance> <line no> <token no in line>"


     Section 2: The add section
     --------------------------

     Contains a specification of the number of tokens that might be added at
     a particular point in the reference file. The information will be stored
     in a dictionary with the reference file token numbers as key. Each entry
     records the maximum number of tokens that may be added at that point. 
     The format of this section is:

     "<token number> <max. number of added tokens> <line no> <token no in line>"

     Section 3: The delete section
     -----------------------------

     Contains a specification of the tokens in reference file that may be
     missing from the data file. The information will once again be stored in 
     a dictionary, but the entries in the dictionary will be empty. Instead
     the fact that a reference file token is mentioned in the dictionary alone
     means that it might be deleted in the data file. The format of this 
     section simply is:

     "<token number> <line no> <token no in line>"

     Section 4: The settings section
     -------------------------------

     This section contains information that is essential to unambiguously 
     establish the relationship between the reference file and the tolerance
     file. At present the only information contained in this section is
     the separator character list. However in future there could be more
     items stored here. All items are stored as name and value pairs.

     The separators are stored following the keyword "separators" and the 
     separator characters themselves are specified in a double quote delimited
     list separated by white space. E.g.

       separators "* % " ,"

     
     Finally, to separate the four sections the following constructs are used:
    
     <section type>
     ...
     end

     Where section type can be any of the keywords "change", "add", "delete"
     or "settings".

     Apart from the tolerance sections every tolerance file will also contain
     a version line. This line specifies the version of the toldiff script
     that last wrote the tolerance file. The version line is produced by
     the function version_toldiff and follows the associated format.

     Note on the handling of guides: Guides are additional dummy tokens to
     direct the diff procedure in favour of matching tokens. The number of
     these additional tokens is specified by nguides. The inclusion of these
     guides has consequences for the way the tolerances will be stored. In 
     particular if the true token number (trutoken) is the one stored in the 
     tolerance file and the effective token (efftoken) number is the number of
     the corresponding token in the internal data structures and the number
     of guides (nguides) is given then:
     - efftoken = (nguides+1)*trutoken-nguides -- is the original token
     - efftoken+1 .. efftoken+nguides          -- are the dummy guide tokens
     In this context the tolerances will be handled as:
     - delete tolerances: 
         if a token may be deleted then the same is true for its guides, so
         mark the token itself and all its guides accordingly.
     - add tolerances:
         if tokens may be added after a token then the same may happen after
         the guide immediately following it. The exception is token 0 which is
         a true virtual token which cannot have a following guide.
     - change tolerances:
         these only ever apply to the true tokens and never to the guides.
     Similar considerations apply when writing the tolerance file.
     """
  #
  # IOError exceptions are handled by the calling routine.
  #
  change = { }
  add    = { }
  delete = { }
  line = fp.readline()
  while line:
    list = string.split(line)
    if len(list) != 0:
      if   list[0] == "change":
        line = fp.readline()
        while (line != "end\n"):
          list             = string.split(line)
          trutoken         = list.pop(0)
          trutoken         = int(trutoken)
          tolerances       = toldiff_tokens.change_tol(list)
          efftoken         = (nguides+1)*trutoken-nguides
          change[efftoken] = tolerances
          line = fp.readline()
        line = fp.readline()
      elif list[0] == "add":
        line = fp.readline()
        while (line != "end\n"):
          list          = string.split(line)
          trutoken      = list[0]
          numtoks       = list[1]
          lineno        = list[2]
          tokno         = list[3]
          trutoken      = int(trutoken)
          numtoks       = int(numtoks)*(nguides+1)
          lineno        = int(lineno)
          tokno         = int(tokno)
          efftoken      = (nguides+1)*trutoken-nguides
          if trutoken > 0:
            add[efftoken] = (numtoks,lineno,tokno)
          else:
            add[0]        = (numtoks,lineno,tokno)
          if trutoken > 0 and nguides > 0:
            iguide = 1
            while (iguide <= nguides):
              add[efftoken+iguide] = (numtoks,lineno,tokno)
              iguide = iguide + 1
          line = fp.readline()
        line = fp.readline()
      elif list[0] == "delete":
        line = fp.readline()
        while (line != "end\n"):
          list           = string.split(line)
          trutoken       = list[0]
          lineno         = list[1]
          tokno          = list[2]
          trutoken       = int(trutoken)
          lineno         = int(lineno)
          tokno          = int(tokno)
          efftoken       = (nguides+1)*trutoken-nguides
          number         = efftoken
          while (number <= efftoken+nguides):
            delete[number] = (lineno,tokno)
            number = number + 1
          line = fp.readline()
        line = fp.readline()
      elif list[0] == "settings":
        line = fp.readline()
        while (line != "end\n"):
          list           = string.split(line)
          if list[0] == "separators":
            ibgn = string.find(line,"\"")
            iend = string.rfind(line,"\"")
            seps = line[ibgn+1:iend]
            separators = string.split(seps)
          line = fp.readline()
        line = fp.readline()
      elif list[0] == "toldiff":
        # this is the version line 
        version = list[2]
        version = int(version)
        if version < 67:
          print line
          print "Unsupported tolerance file format encountered\n"
          print "Tolerance file too old\n"
          exit(11)
        elif version > 75:
          print line
          print "Unsupported tolerance file format encountered\n"
          print "Tolerance file too new\n"
          exit(12)
        line = fp.readline()
      elif list[0] == "#":
        line = fp.readline()
      else:
        print line
        print "Unknown tolerance file section encountered\n"
        exit(10)
    else:
      line = fp.readline()
  return (change,add,delete,separators)

def save_tolerances(fp,change,add,delete,errfp,separators,nguides):
  """Saves the tolerances specified in the change, add and delete dictionaries.
     The file pointer fp is used to write the data to. See load_tolerances
     for a specification of the file format.

     To handle the guides a few modifications are needed. What they are 
     depends on the nature of the tolerance:
     - delete tolerances:
         with any delete tolerances there must be a true token that is deleted
         as well. So we store only the delete tolerances associated with the
         true tokens.
     - add tolerances:
         add tolerances may be associated with the true token or the guides that
         immediately follow it. So the add tolerance is computed as the 
         maximum number of tokens that may be added after the true token or
         the guide immediately following it.
     - change tolerances:
         these are stricktly associated with the true tokens only. So only
         store change tolerances for the true tokens.
     """
  version_toldiff(fp,errfp)
  #
  # Beyond this point IOError exceptions are handled by the calling routine
  # above this one.
  #
  fp.write("settings\n")
  nseps = len(separators)
  isep  = 0
  fp.write("  separators \"")
  while (isep < nseps-1):
    fp.write(separators[isep]+" ")
    isep = isep + 1
  if isep < nseps:
    fp.write(separators[isep])
  fp.write("\"\n")
  fp.write("end\n")
  fp.write("change\n")
  list = change.keys()
  list.sort()
  for number in list:
    if mod(number+nguides,nguides+1) == 0:
      (type,tol,lineno,tokno) = change[number]
      trutoken = (number+nguides)/(nguides+1)
      fp.write("  %d %s %s %d %d\n" % (trutoken,type,str(tol),lineno,tokno))
  fp.write("end\n")
  fp.write("add\n")
  #
  # Transfer add tolerances onto the true token
  #
  list = add.keys()
  list.sort()
  for number in list:
    (numtoks,lineno,tokno) = add[number]
    trutoken = (number+nguides)/(nguides+1)
    efftoken = (nguides+1)*trutoken-nguides
    if add.has_key(efftoken):
      (numtoks1,lineno,tokno) = add[efftoken]
      numtoks = max(numtoks,numtoks1)
    add[efftoken] = (numtoks,lineno,tokno)
  #
  # Now store the add tolerances. We only need to look at the true tokens now.
  # Rounding up of the number of add tokens is needed if the guides mismatch
  # but the token does not.
  #
  list = add.keys()
  list.sort()
  for number in list:
    if number == 0:
      (numtoks,lineno,tokno) = add[number]
      numtoks = (numtoks+nguides)/(nguides+1)
      trutoken = 0
      fp.write("  %d %d %d %d\n" % (trutoken,numtoks,lineno,tokno))
    elif mod(number+nguides,nguides+1) == 0:
      (numtoks,lineno,tokno) = add[number]
      trutoken = (number+nguides)/(nguides+1)
      numtoks = (numtoks+nguides)/(nguides+1)
      fp.write("  %d %d %d %d\n" % (trutoken,numtoks,lineno,tokno))
  fp.write("end\n")
  fp.write("delete\n")
  #
  # Transfer delete tolerances onto the true token 
  # This may seem unexpected but is needed when e.g. an additional newline
  # causes the guides to mismatch but the token actually does not. 
  #
  list = delete.keys()
  list.sort()
  for number in list:
    (lineno,tokno) = delete[number]
    trutoken = (number+nguides)/(nguides+1)
    efftoken = (nguides+1)*trutoken-nguides
    delete[efftoken] = (lineno,tokno)
  #
  # Now store the delete tolerances. We only need to look at the true tokens
  # now.
  #
  list = delete.keys()
  list.sort()
  for number in list:
    if mod(number+nguides,nguides+1) == 0:
      (lineno,tokno) = delete[number]
      trutoken = (number+nguides)/(nguides+1)
      fp.write("  %d %d %d\n" % (trutoken,lineno,tokno))
  fp.write("end\n")

def save_tokenized(fp,refdat,errfp):
  """Save tokenized reference or data file for the external diff program
     to compare. It assumed that the file opening and closing is handled by
     the calling routine."""
  Ib = 1
  Ie = len(refdat.token)
  i  = Ib-1
  while (i < Ie):
    i = i + 1
    (type,info,line,tokno) = refdat.token[i]
    info = str(info)
    fp.write("%s\n" % info)
