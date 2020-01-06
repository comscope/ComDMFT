#!/usr/bin/env python
#
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

import os
import sys
import string
import toldiff_files
import toldiff_lcs
import toldiff_diff
import toldiff_update
import toldiff_transfer
import toldiff_show
import toldiff_tokens

def max(a,b):
  """Return the maximum value of the two arguments"""
  if a >= b:
    result = a
  else:
    result = b
  return result

def license_toldiff(fp,errfp):
  """Print out the license information to the specified file object."""
  try:
    fp.write("""
    Copyright (C) 2006 Huub van Dam, Science and Technology Facilities Council,
    Daresbury Laboratory.
    All rights reserved.

    Developed by:        Huub van Dam
                         Science and Technology Facilities Council
                         Daresbury Laboratory
                         Computational Science and Engineering Department
                         Computational Chemistry Group
                         http://www.cse.clrc.ac.uk/ccg

    Permission is hereby granted, free of charge, to any person obtaining 
    a copy of this software and associated documentation files (the "Software"),
    to deal with the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense, 
    and/or sell copies of the Software, and to permit persons to whom the 
    Software is furnished to do so, subject to the following conditions: 

    Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimers. 
    Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimers in the documentation 
    and/or other materials provided with the distribution. 
    Neither the names of the Science and Technology Facilities Council,
    Daresbury Laboratory, the Computational Science and Engineering Department,
    the Computational Chemistry Group, nor the names of its contributors may be
    used to endorse or promote products derived from this Software without
    specific prior written permission. 

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS WITH THE SOFTWARE. 
    \n""")
    sys.exit(1)
  except IOError, e:
    (errno,errmsg) = e
    try:
      errfp.write("toldiff: error writing license information\n")
      errfp.write("toldiff: error message: ")
      errfp.write(errmsg)
      errfp.write("\n")
    except IOError, e:
      pass
    sys.exit(5)

def usage_toldiff(fp,errfp):
  """Print out the usage information to the specified file object."""
  try:
    fp.write("""
  Usage:

     toldiff [[--diff] <reference file> <data file>]
             [--update <reference file> <data file>]
             [--transfer <reference file> <new reference file>]
             [--show <reference file>]
             [--tolerance <tolerance file name>]
             [--new-tolerance <new tolerance file name>]
             [--diff-exe <diff executable>]
             [--output full|summary|none]
             [--summary <identical>:<equivalent>:<different>]
             [--exit <identical>:<equivalent>:<different>]
             [--[no]exact] [--[no]tolerant] [--[no]best]
             [--itol-scale <integer tolerance scale factor>]
             [--ftol-scale <floating point tolerance scale factor>]
             [--ctol-scale <complex tolerance scale factor>]
             [--separators <separator character list>]
             [--guides <number of guides>]
             [--[no]backtrack]
             [--help] [--license] [--version]

  Toldiff is a script that compares two files allowing for tolerable
  differences. Tolerable differences often arise in meta data like who ran
  the test and on which date, timing data, and which machines and how many
  processors were used. In scientific/technical codes additional variations
  may result from numerical accuracy limitations.

  Toldiff is designed to assist in software testing by suppressing tolerable
  or trivial differences and highlighting only the significant ones. This
  facilitates checking whether an output a program has just produced matches
  the reference result obtained in the past.

  The toldiff script knows of the following files:

  A. The reference file:
     - THE correct file
  B. The data file: 
     - A file the correctness of which is to be tested against the reference
       file. Once its correctness has been established it may be used to update
       the tolerances.
  C. The tolerance file:
     - This file records where all allowed differences can occur, if any.
  D. The new reference file:
     - The file that is to replace the reference file after a change has
       taken place that outdates the reference file
  E. The new tolerance file:
     - This file records where allowed differences can occur relative to the
       new reference file instead of the current reference file.
    
  The script offers three processes:

  1. The diff process:
     - This process reports all differences between the reference file and
       the data file that are not explicitly tolerated.
  2. The update process:
     - This process updates the tolerances file adding all differences between 
       the reference file and the data file that were not tolerated before.
  3. The transfer process:
     - If the current reference file needs to be replaced by a new one this 
       process will carry as many as possible known tolerances relative to the
       current reference file over to the new reference file.

  There are various command line options to control toldiff. In cases where
  environment variables can be used as an alternative to command line options
  the precedence is handled as:
  - environment variables take precedence over default settings
  - command line options take precedence over environment variables.

  There are three categories of options this script will recognise:

  1. Process options:

    1.1 --diff <reference file name> <data file name>

        This triggers the script to perform the default diff process of
        comparing the data file against the reference file.

    1.2 --update <reference file name> <data file name>

        This requests the update process to be performed updating the
        tolerance file to allow for any differences between the reference and
        data files.

        During this process the new tolerances computed can be scaled by a
        factor that is equal to or larger than one. This may be useful when
        the expected fluctuations are larger than the current differences.
        Separate scale factors may be set for each of the three different
        numerical data types supported, i.e. integer, floating point, and
        complex. The scale factors are always floating point numbers but
        after scaling the tolerance the result is rounded where appropriate.

       1.2.1 --itol-scale <integer tolerance scale factor>

         Sets the scale factor for integer tolerances.

       1.2.2 --ftol-scale <floating point tolerance scale factor>

         Sets the scale factor for floating point tolerances.

       1.2.3 --ctol-scale <complex tolerance scale factor>

         Sets the scale factor for complex tolerances.

    1.3 --tolerance <tolerance file name>

        This option allows explicit specification of the tolerance file name.
        If omitted the script will construct a name for the tolerance file 
        from the name of the reference file.

    1.4 --transfer <reference file name> <new reference file name>

        This option invokes the transfer process to migrate as many tolerances
        as possible from the current reference file over to the new one.

    1.5 --new-tolerance <new tolerance file name>

        This option allows for the explicit specification of the name of the
        new tolerance file. If this is omitted the script will construct a name
        for the new tolerance file from the new reference file name.

    1.6 --diff-exe <diff executable>

        This option enbles replacing some of the Python diff implementation
        by invoking a binary diff program. This greatly improves the
        performance without changing the functionality. As an alternative 
        mechanism the environment variable TOLDIFF_EXE may be set to specify
        the diff program. In case both the command line option and the
        environment variable are provided the command line option has 
        precedence.

    1.7 --output full|summary|none

        This option controls the amount of output toldiff produces. The default
        setting "full" results in printing a full diff output. The setting
        "summary" suppresses the diff output and replaces it with a short 
        string for files being identical, equivalent or different. The values
        of these strings can be specified with the --summary option. Finally,
        setting "none" suppresses all output. Other than the --output option
        setting the TOLDIFF_OUTPUT environment variable does the same.

    1.8 --summary <identical>:<equivalent>:<different>

        This option allows the specification of short results for toldiff. The
        first string is reported if the reference file and data file are 
        identical. The second string is reported if the reference and data files
        are not identical but all differences are tolerated. The last string
        is reported if there are differences that are not tolerated. The
        default strings are "identical", "equivalent", and "different". Finally,
        these settings can be specified by setting the TOLDIFF_SUMMARY
        environment variable. In both ways the values are colomn separated.

    1.9 --exit <identical>:<equivalent>:<different>

        This option specifies the exit codes for toldiff. The first value is
        reported if the reference file and data file are identical. The second
        value is reported if the reference and data files are not identical but
        all differences are tolerated. The last value is reported if there are
        differences that are not tolerated. The default values are 0, 0, and 1.
        Finally, these settings can be specified by setting the TOLDIFF_EXIT
        environment variable. In both ways the values are colomn separated.

    1.10 --separators <separator character list>

        Toldiff splits the data in the reference file and the data file into
        tokens. It always uses white space to separate tokens. However it may
        be necessary to break the tokens up further. It uses any characters
        in the separator character list for that purpose. As the tolerances
        depend on the separator character list this list can only be specified
        when the tolerance file is created. In all other instances specifying
        this list will be ignored.

        Of course there is the potential to discover that the current set of
        separator characters stored in the tolerance file is not optimal.
        In that case the transfer process can be used to create a new tolerance
        file based on a new set of separators. The specified separator list
        will be used to create the new tolerance file.

        The separator character list is specified as a white space separated
        list of characters, e.g.

           --separators "% = ,"

        Alternatively the separator character list may be specified using
        the environment variable TOLDIFF_SEPARATORS.

    1.11 --guides <number of guides>

        Tokens are typically short character sequences. As a result if a token
        has changed there is a significant chance it will accidentally match
        another token. This results in rather unexpected tolerances. Guides
        are dummy tokens that direct the diff process to match tokens correctly
        even if the tokens do not match exactly. The number of guides used 
        determines strict this enforcement is, 0 means no enforcement, 2 means
        maximum enforcement. Alternatively the environment variable
        TOLDIFF_GUIDES may be used.

    1.12 --[no]backtrack

        Another way to deal with the issue discussed under --guides is to let
        the tolerant diff procedure re-analyse some of the differences found
        initially. Initially a traditional diff procedure is used that finds
        exact matches. As this cannot take tolerances into account suboptimal
        matches may result. Rather than rigidly adhering to the matches the
        initial diff has found the --backtrack option extends the differences
        to the nearest number of whole lines. These whole line sections are
        then re-analysed using the tolerant diff procedure, thus allowing
        matches to be found that the initial diff by design cannot find.
        The environment variable TOLDIFF_BACKTRACK may be used instead of the
        command line flag.

    Both the --guides and --backtrack options are designed to deal with the
    situation where adjacent tokens have overlapping ranges of valid values.
    However, even in these situations unintended matches are unlikely unless
    the values have very few relevant digits. I.e. is the tolerance is such 
    that only 1 digit may change then the chance of accidently matching a
    neighbouring number is 1 in 10, if 3 digits may change then the chance is
    1 in 1000. As a result one may want to check whether the extra expense of
    using the --guides and --backtrack options is justified given the 
    associated risk.

  2. Information options:

    2.1 --help

        Print this information on how to use this scripts.

    2.2 --show <reference file name>

        Prints the reference file marking all the known tolerances on it.
        This allows checking how the program has resolved differences through
        the tolerances chosen.
        The tolerances are marked on each line in the following order:
        1. The number of lines that may be inserted after this line.
        2. Whether this line may be deleted in which case it will be marked by
           a 'X', otherwise white space indicates that the line has to be 
           present.
        3. The contents of the line are shown with those characters that may
           change replaced by '#'.

    2.3 --version

        Print the version number of the toldiff script you are using.

    2.4 --license

        Print the license conditions under which this script is distributed.

  3. Debug options:

    These options are normally set automatically based on the requirements of
    the selected process. The default settings aim to complete the selected
    process with the highest efficiency. However, for debugging purposes it
    is possible to override these settings. You are free to try them to your
    own peril.

    3.1 --[no]exact

        Enable or disable the file differencing procedure that is based on
        exact line matches.

    3.2 --[no]tolerant

        Enable or disable the file differencing procedure that uses a line
        comparison which allows for tolerable differences between lines.

    3.3 --[no]best

        Enable or disable the file differencing procedure that matches lines
        based on maximum similarity.

  Copyright 2006, Huub van Dam, Science and Technology Facilities Council,
  Daresbury Laboratory\n""")
    sys.exit(1)
  except IOError, e:
    (errno,errmsg) = e
    try:
      errfp.write("toldiff: error writing usage information\n")
      errfp.write("toldiff: error message: ")
      errfp.write(errmsg)
      errfp.write("\n")
    except IOError, e:
      pass
    sys.exit(5)

def load_file(filename,err_fp,separators,nguides):
  """Open and load a file. Returns the file text and the number of lines.
     The routine also handles I/O errors. I.e. it reports the error to the
     user and terminates the program.

     When the file is read the appropriate number of guides are inserted
     as specified by nguides.
     """
  text = toldiff_tokens.tokenized_file()
  lines = 0
  tokens = 0
  try:
    file_fp = open(filename,"r")
    (text,lines,tokens) = toldiff_files.load_plain_text(file_fp,text,lines,tokens,separators,nguides)
    file_fp.close()
  except IOError, e:
    (errno,errmsg) = e
    try:
      err_fp.write("toldiff: I/O error on file: ")
      err_fp.write(filename)
      err_fp.write("\n")
      err_fp.write("toldiff: I/O error message: ")
      err_fp.write(errmsg)
      err_fp.write("\n")
    except IOError, e:
      pass
    sys.exit(10)
  return (text,lines,tokens)

def store_tolerance(tol_fnm,chg_txt,add_txt,del_txt,err_fp,separators,nguides):
  """Open and write the tolerance file. The routine handles any I/O errors.
     I.e. it reports the error to the user and terminates the program."""
  try:
    tol_fp = open(tol_fnm,"w")
    toldiff_files.save_tolerances(tol_fp,chg_txt,add_txt,del_txt,err_fp,separators,nguides)
    tol_fp.close()
  except IOError, e:
    (errno,errmsg) = e
    try:
      err_fp.write("toldiff: I/O error encountered attempting to write: ")
      err_fp.write(tol_fnm)
      err_fp.write("\n")
      err_fp.write("toldiff: I/O error message: ")
      err_fp.write(errmsg)
      err_fp.write("\n")
    except IOError, e:
      pass
    sys.exit(30)

def run_diff(diff_exe,ref_fnm,dat_fnm,ref,dat,fp):
  """This routine starts off an external diff program. 

     As the tokenized versions of the reference and data files do not exist
     these have to be written first. Next the diff program is started.
     Both the stdout and stderr file descriptors are returned as due file
     buffer space the diff program cannot complete if stdout is not read. 
     So only after reading stdout to drive diff to completion can stderr be
     checked to see if diff ran successfully.

     If an error is reported on stderr this should be passed on to the user
     and the program should terminate.

     After diff has run the tokenized files should be deleted.

     - diff_exe - the path of the diff executable
     - ref_fnm  - the filename for the temporary tokenized reference file
     - dat_fnm  - the filename for the temporary tokenized data file
     - ref      - the tokenized reference
     - dat      - the tokenized data
     - fp       - a file descriptor for error reporting
     """
  cmd = diff_exe+" "+ref_fnm+" "+dat_fnm

  try:
    ref_fp = open(ref_fnm,"w")
    toldiff_files.save_tokenized(ref_fp,ref,fp)
    ref_fp.close()
  except IOError, e:
    (errno,errmsg) = e
    try:
      fp.write("toldiff: I/O error on tokenized reference file\n")
      fp.write("toldiff: I/O error message: ")
      fp.write(errmsg)
      fp.write("\n")
    except IOError, e:
      pass
    sys.exit(25)

  try:
    dat_fp = open(dat_fnm,"w")
    toldiff_files.save_tokenized(dat_fp,dat,fp)
    dat_fp.close()
  except IOError, e:
    (errno,errmsg) = e
    try:
      fp.write("toldiff: I/O error on tokenized data file\n")
      fp.write("toldiff: I/O error message: ")
      fp.write(errmsg)
      fp.write("\n")
    except IOError, e:
      pass
    sys.exit(25)

  try:
    (in_fp,out_fp,err_fp) = os.popen3(cmd)
  except IOError, e:
    (errno,errmsg) = e
    try:
      fp.write("toldiff: I/O error on external diff standard error file\n")
      fp.write("toldiff: I/O error message: ")
      fp.write(errmsg)
      fp.write("\n")
    except IOError, e:
      pass
    sys.exit(25)
  in_fp.close()
  return (out_fp,err_fp)
      
def find_overall_lcs(lexact,ltol,lbest,tol,ref_fnm,dat_fnm,diff_exe,feps,ieps,err_fp,separators,nguides,snake_trim,update):
  """Find the overall LCS including the tolerances. The general procedure is
     simply to establish the exact LCS, then try to resolve as much of the
     mismatches by considering the tolerances, then try to match the remaining
     differences to minimize the mismatches.

     This routine will read in the reference file and the data file as well.
     The reason for this is that this is more efficient in case an external
     diff program is used for the first phase.

     The routine returns the overall LCS, the reference file text, the data
     file text and beginning and ending token numbers of both files.

     This routine allows each phase to be disabled explicitly through a
     flag passed in as an argument:
     - lexact: if false skip the exact matching
     - ltol  : if false skip the tolerant matching
     - lbest : if false skip the minimal difference matching.

     The number of guides is specified in nguides. This is used in reading
     in the reference and data files.
     """
  lcs = [ ]
  if lexact:
    if diff_exe == "":
      Nb  = 1
      Ntb = 1
      Mb  = 1
      Mtb = 1
      (ref,Ne,Nte) = load_file(ref_fnm,err_fp,separators,nguides)
      (dat,Me,Mte) = load_file(dat_fnm,err_fp,separators,nguides)
      lcs = toldiff_lcs.find_lcs1(ref,Ntb,Nte,dat,Mtb,Mte)
    else:
      error = false
      Nb  = 1
      Ntb = 1
      Mb  = 1
      Mtb = 1
      (ref,Ne,Nte) = load_file(ref_fnm,err_fp,separators,nguides)
      (dat,Me,Mte) = load_file(dat_fnm,err_fp,separators,nguides)
      #
      # Construct temporary file names
      #
      pid = os.getpid()
      # The extra "a" and "b" ensure unique file names even if the reference
      # and data file names are the same.
      tmp_ref_fnm = ref_fnm+"a"+str(pid)
      tmp_dat_fnm = dat_fnm+"b"+str(pid)
      #
      # Construct temporary files, invoke diff and parse diff output
      #
      (diff_out_fp,diff_err_fp) = run_diff(diff_exe,tmp_ref_fnm,tmp_dat_fnm,ref,dat,err_fp)
      lcs = toldiff_diff.diff_to_lcs(Ntb,Nte,Mtb,Mte,diff_out_fp,err_fp)
      diff_out_fp.close()
      #
      # Delete temporary files
      #
      os.remove(tmp_ref_fnm)
      os.remove(tmp_dat_fnm)
      #
      # Check whether the diff program detected any errors
      #
      try:
        line = diff_err_fp.readline()
        while line:
          error = true
          err_fp.write("toldiff:"+line)
          line = diff_err_fp.readline()
        diff_err_fp.close()
      except IOError, e:
        (errno,errmsg) = e
        try:
          err_fp.write("toldiff: I/O error on external diff standard error file\n")
          err_fp.write("toldiff: I/O error message: ")
          err_fp.write(errmsg)
          err_fp.write("\n")
        except IOError, e:
          pass
        sys.exit(25)
      if error:
        sys.exit(20)
  else:
    Nb  = 1
    Ntb = 1
    Mb  = 1
    Mtb = 1
    (ref,Ne,Nte) = load_file(ref_fnm,err_fp,separators,nguides)
    (dat,Me,Mte) = load_file(dat_fnm,err_fp,separators,nguides)
  #Snake trimming may only be used here!
  if (snake_trim and ((ltol and len(tol) > 0 ) or (update and lbest))):
    lcs = toldiff_lcs.trim_snakes(lcs,ref,Ntb,Nte,dat,Mtb,Mte)
  #
  if (len(tol) <= 0) or (not ltol):
    #
    # No tolerances were specified or this phase is explicitly suppressed
    #
    pass
    #
  else:
    #
    # Consider all the differences and try to resolve as many as possible.
    #
    if (len(lcs) <= 0):
      #
      # Then the new LCS is simply the result of the tolerant diff
      #
      lcs = toldiff_lcs.find_lcs2(tol,ref,Ntb,Nte,dat,Mtb,Mte,feps,ieps)
      #
    else:
      #
      # First consider whether there is anything to compare before the first
      # snake
      #
      lcs1 = lcs
      (xbot1,ybot1,xtop1,ytop1,type1) = lcs1.pop(0)
      if (xbot1 > Mtb) and (ybot1 > Ntb):
        lcs = toldiff_lcs.find_lcs2(tol,ref,Ntb,ybot1-1,dat,Mtb,xbot1-1,feps,ieps)
      else:
        lcs = [ ]
      xtop0 = xtop1
      ytop0 = ytop1
      lcs.append((xbot1,ybot1,xtop1,ytop1,type1))
      while (len(lcs1) > 0 ):
        (xbot1,ybot1,xtop1,ytop1,type1) = lcs1.pop(0)
        if (xbot1 > xtop0+1) and (ybot1 > ytop0+1):
          lcs2 = toldiff_lcs.find_lcs2(tol,ref,ytop0+1,ybot1-1,dat,xtop0+1,xbot1-1,feps,ieps)
          lcs = lcs + lcs2
        xtop0 = xtop1
        ytop0 = ytop1
        lcs.append((xbot1,ybot1,xtop1,ytop1,type1))
      if (Nte >= ytop0+1) and (Mte >= xtop0+1):
        #
        # The some more stuff at the end left to do
        #
        lcs2 = toldiff_lcs.find_lcs2(tol,ref,ytop0+1,Nte,dat,xtop0+1,Mte,feps,ieps)
        lcs = lcs + lcs2

  if (not lbest):
    #
    # This phase is explicitly suppressed
    #
    pass
    #
  else:
    #
    # Consider all the differences and try to match different lines as best as
    # possible minimizing the number of differences.
    #
    #Snake trimming does not work here as the lcs3 may pair tokens up in a way
    #that is different from what lcs2 would do. The result of this inconsistency
    #is that some differences will never be tolerated! Clearly this breaks
    #toldiff.
    #lcs = toldiff_lcs.trim_snakes(lcs,ref,Ntb,Nte,dat,Mtb,Mte)
    if (len(lcs) <= 0):
      #
      # Then the new LCS is simply the result of the best match diff,
      # which will probably hurt as this will get very expensive.
      #
      lcs = toldiff_lcs.find_lcs3(tol,ref,Ntb,Nte,dat,Mtb,Mte,feps,ieps)
      #
    else:
      #
      # First consider whether there is anything to compare before the first
      # snake
      #
      lcs1 = lcs
      (xbot1,ybot1,xtop1,ytop1,type1) = lcs1.pop(0)
      if (xbot1 > Mtb) and (ybot1 > Ntb):
        lcs = toldiff_lcs.find_lcs3(tol,ref,Ntb,ybot1-1,dat,Mtb,xbot1-1,feps,ieps)
      else:
        lcs = [ ]
      xtop0 = xtop1
      ytop0 = ytop1
      lcs.append((xbot1,ybot1,xtop1,ytop1,type1))
      while (len(lcs1) > 0 ):
        (xbot1,ybot1,xtop1,ytop1,type1) = lcs1.pop(0)
        if (xbot1 > xtop0+1) and (ybot1 > ytop0+1):
          lcs2 = toldiff_lcs.find_lcs3(tol,ref,ytop0+1,ybot1-1,dat,xtop0+1,xbot1-1,feps,ieps)
          lcs = lcs + lcs2
        xtop0 = xtop1
        ytop0 = ytop1
        lcs.append((xbot1,ybot1,xtop1,ytop1,type1))
      if (Nte >= ytop0+1) and (Mte >= xtop0+1):
        #
        # There is some more stuff at the end left to do
        #
        lcs2 = toldiff_lcs.find_lcs3(tol,ref,ytop0+1,Nte,dat,xtop0+1,Mte,feps,ieps)
        lcs = lcs + lcs2
  return (lcs,ref,Ntb,Nte,dat,Mtb,Mte)

def construct_tolerance_filename(ref_fnm,dat_fnm,tol_fnm):
  if tol_fnm == "":
    #
    # No tolerance file name given so try and construct one
    #
    i = string.rfind(ref_fnm,".")
    if i == -1:
      #
      # The reference file name has no extension.
      #
      ref_ext = ""
      i = len(ref_fnm)+1
      #
    else:
      #
      # The reference file name has an extension extract it.
      #
      ref_ext = ref_fnm[i:]
      #
    j = string.rfind(dat_fnm,".")
    if j == -1:
      #
      # The data file name has no extension.
      #
      dat_ext = ""
      #
    else:
      #
      # The data file name has an extension extract it.
      #
      dat_ext = dat_fnm[j:]
      #
    tol_ext = ".tol"
    if (tol_ext == ref_ext) or (tol_ext == dat_ext):
      tol_ext = ".tlr"
      if (tol_ext == ref_ext) or (tol_ext == dat_ext):
        tol_ext = ".tlc"
    tol_fnm = ref_fnm[:i]+tol_ext
  return tol_fnm

true = (0 == 0)
false = not true
#
# Set up the default comparison options
#
diff     = 1
update   = 2
transfer = 3
show     = 4
process  = diff
lexact   = true
ltol     = true
lbest    = false
#
# Set up default comparison results
#
identical  = 1
equivalent = 2
different  = 3
#
# Set up default comparison exit codes
#
exit_identical  = 0
exit_equivalent = 0
exit_different  = 1
#
# Set up default comparison summary texts
#
text_identical  = "identical"
text_equivalent = "equivalent"
text_different  = "different"
#
# Set up output options and default output option
#
output_full    = 3
output_summary = 2
output_none    = 1
output         = output_full
#
# Set up the default list of separator characters for the reference and
# data file tokenisation. In addition to these characters whitespace will
# be used as token separator as well. Note that a long list of separators
# deteriorates the performance significantly.
#
# Separators is the list of additional separator characters used for the
# reference file and the data file.
# Separators_new is the list of additional separator characters used for the
# new reference file in case of a transfer operation.
#
separators = []
separators_new = []
#
# Set the default snake trimming behaviour
#
snake_trim = false
#
# Set the default number of guides
#
nguides = 0
#
# Set up default precisions for floating point and integer numbers
#
feps = 1.0e-12
ieps = 0.1
#
# Set up default scale factors for new tolerances
#
tol_scale  = 1.0
itol_scale = tol_scale
ftol_scale = tol_scale
ctol_scale = tol_scale

lcs = [ ]
diff_exe    = ""
tol_fnm     = ""
tol_new_fnm = ""
ref_fnm     = ""
dat_fnm     = ""
narg = len(sys.argv)
iarg = 1

if os.environ.has_key("TOLDIFF_EXE"):
  diff_exe = os.environ["TOLDIFF_EXE"]

if os.environ.has_key("TOLDIFF_OUTPUT"):
  output = os.environ["TOLDIFF_OUTPUT"]
  if output == "FULL" or output == "full":
    output = output_full
  elif output == "SUMMARY" or output == "summary":
    output = output_summary
  elif output == "NONE" or output == "none":
    output = output_none

if os.environ.has_key("TOLDIFF_EXIT"):
  exit_codes = os.environ["TOLDIFF_EXIT"]
  exit_codes = string.split(exit_codes,":")
  if len(exit_codes) == 3:
    exit_identical  = int(exit_codes[0])
    exit_equivalent = int(exit_codes[1])
    exit_different  = int(exit_codes[2])

if os.environ.has_key("TOLDIFF_SUMMARY"):
  text_summaries = os.environ["TOLDIFF_SUMMARY"]
  text_summaries = string.split(text_summaries,":")
  if len(text_summaries) == 3:
    text_identical  = text_summaries[0]
    text_equivalent = text_summaries[1]
    text_different  = text_summaries[2]

if os.environ.has_key("TOLDIFF_ITOLSCALE"):
  itol_scale = max(tol_scale,float(os.environ["TOLDIFF_ITOLSCALE"]))

if os.environ.has_key("TOLDIFF_FTOLSCALE"):
  ftol_scale = max(tol_scale,float(os.environ["TOLDIFF_FTOLSCALE"]))

if os.environ.has_key("TOLDIFF_CTOLSCALE"):
  ctol_scale = max(tol_scale,float(os.environ["TOLDIFF_CTOLSCALE"]))

if os.environ.has_key("TOLDIFF_SEPARATORS"):
  separators     = string.split(os.environ["TOLDIFF_SEPARATORS"])
  separators_new = string.split(os.environ["TOLDIFF_SEPARATORS"])

if os.environ.has_key("TOLDIFF_GUIDES"):
  nguides = max(0,int(os.environ["TOLDIFF_GUIDES"]))

if os.environ.has_key("TOLDIFF_BACKTRACK"):
  tmptxt = os.environ["TOLDIFF_BACKTRACK"]
  tmptxt = tmptxt.lower()
  if   tmptxt == "yes" or tmptxt == "y":
    snake_trim = true
  elif tmptxt == "no"  or tmptxt == "n":
    snake_trim = false
  else:
    try:
      sys.stderr.write("toldiff: invalid value for TOLDIFF_BACKTRACK should be \"yes\" or \"no\"\n")
    except IOError, e:
      pass
    sys.exit(5)
    

if narg == 1:
  usage_toldiff(sys.stdout,sys.stderr)
while iarg < narg:
  if sys.argv[iarg]   == "--exact":
    lexact = true
  elif sys.argv[iarg] == "--noexact":
    lexact = false
  elif sys.argv[iarg] == "--tolerant":
    ltol   = true
  elif sys.argv[iarg] == "--notolerant":
    ltol   = false
  elif sys.argv[iarg] == "--best":
    lbest  = true
  elif sys.argv[iarg] == "--nobest":
    lbest  = false
  elif sys.argv[iarg] == "--tolerance":
    iarg = iarg + 1
    if iarg < narg:
      tol_fnm = sys.argv[iarg]
    else:
      try:
        sys.stderr.write("toldiff: missing tolerance file name\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--new-tolerance":
    iarg = iarg + 1
    if iarg < narg:
      tol_new_fnm = sys.argv[iarg]
    else:
      try:
        sys.stderr.write("toldiff: missing new tolerance file name\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--diff-exe":
    iarg = iarg + 1
    if iarg < narg:
      diff_exe = sys.argv[iarg]
    else:
      try:
        sys.stderr.write("toldiff: missing diff executable specification\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--diff":
    process = diff
  elif sys.argv[iarg] == "--update":
    process = update
  elif sys.argv[iarg] == "--transfer":
    process = transfer
  elif sys.argv[iarg] == "--show":
    process = show
  elif sys.argv[iarg] == "--version":
    toldiff_files.version_toldiff(sys.stdout,sys.stderr)
    sys.exit(0)
  elif sys.argv[iarg] == "--help":
    usage_toldiff(sys.stdout,sys.stderr)
  elif sys.argv[iarg] == "--license":
    license_toldiff(sys.stdout,sys.stderr)
  elif sys.argv[iarg] == "--exit":
    iarg = iarg + 1
    if iarg < narg:
      exit_codes = sys.argv[iarg]
      exit_codes = string.split(exit_codes,":")
      if len(exit_codes) == 3:
        exit_identical  = int(exit_codes[0])
        exit_equivalent = int(exit_codes[1])
        exit_different  = int(exit_codes[2])
    else:
      try:
        sys.stderr.write("toldiff: missing exit codes specification\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--summary":
    iarg = iarg + 1
    if iarg < narg:
      text_summaries = sys.argv[iarg]
      text_summaries = string.split(text_summaries,":")
      if len(text_summaries) == 3:
        text_identical  = text_summaries[0]
        text_equivalent = text_summaries[1]
        text_different  = text_summaries[2]
    else:
      try:
        sys.stderr.write("toldiff: missing summaries specification\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--output":
    iarg = iarg + 1
    if iarg < narg:
      output = sys.argv[iarg]
      if output == "FULL" or output == "full":
        output = output_full
      elif output == "SUMMARY" or output == "summary":
        output = output_summary
      elif output == "NONE" or output == "none":
        output = output_none
      else:
        sys.stderr.write("toldiff: unknown output specification: %s\n" % output)
        sys.exit(5)
    else:
      try:
        sys.stderr.write("toldiff: missing output specification\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--itol-scale":
    iarg = iarg + 1
    if iarg < narg:
      itol_scale = max(tol_scale,float(sys.argv[iarg]))
    else:
      try:
        sys.stderr.write("toldiff: missing integer tolerance scale factor\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--ftol-scale":
    iarg = iarg + 1
    if iarg < narg:
      ftol_scale = max(tol_scale,float(sys.argv[iarg]))
    else:
      try:
        sys.stderr.write("toldiff: missing floating point tolerance scale factor\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--ctol-scale":
    iarg = iarg + 1
    if iarg < narg:
      ctol_scale = max(tol_scale,float(sys.argv[iarg]))
    else:
      try:
        sys.stderr.write("toldiff: missing complex tolerance scale factor\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--separators":
    iarg = iarg + 1
    if iarg < narg:
      separators     = string.split(sys.argv[iarg])
      separators_new = string.split(sys.argv[iarg])
      i = 0
      n = len(separators)
      while (i < n):
        if len(separators[i]) != 1:
          sys.stderr.write("toldiff: separator character list is not a list of single characters\n")
          sys.stderr.write("toldiff: --separators \""+sys.argv[iarg]+"\"\n")
          sys.exit(5)
        i = i + 1
  elif sys.argv[iarg] == "--guides":
    iarg = iarg + 1
    if iarg < narg:
      nguides = max(0,int(sys.argv[iarg]))
    else:
      try:
        sys.stderr.write("toldiff: missing number of guides\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--backtrack":
    snake_trim = true
  elif sys.argv[iarg] == "--nobacktrack":
    snake_trim = false
  else:
    argstr = sys.argv[iarg]
    if (process < show) and (iarg == narg-2):
      ref_fnm = sys.argv[iarg]
      iarg = iarg + 1
      dat_fnm = sys.argv[iarg]
    elif (process == show) and (iarg == narg-1):
      ref_fnm = sys.argv[iarg]
    elif argstr[0:1] == "-":
      try:
        sys.stderr.write("toldiff: unknow option encountered: ")
        sys.stderr.write(argstr)
        sys.stderr.write("\n")
      except IOError, e:
        pass
      sys.exit(8)
    else:
      sys.stderr.write("toldiff: missing reference or data files?\n")
      sys.exit(9)
  iarg = iarg + 1


if ref_fnm == "":
  sys.stderr.write("toldiff: error: no reference filename given\n")
  sys.exit(5)

if (process < show) and (dat_fnm == ""):
  sys.stderr.write("toldiff: error: no data filename given\n")
  sys.exit(6)

tol_fnm = construct_tolerance_filename(ref_fnm,dat_fnm,tol_fnm)
if process == transfer:
  tol_new_fnm = construct_tolerance_filename(dat_fnm,ref_fnm,tol_new_fnm)

ref_txt   = { }
dat_txt   = { }
chg_txt   = { }
add_txt   = { }
del_txt   = { }
ref_lines = 0
dat_lines = 0

try:
  tol_fp  = open(tol_fnm,"r")
  (chg_txt,add_txt,del_txt,separators) = toldiff_files.load_tolerances(tol_fp,separators,nguides)
  tol_fp.close()
except IOError, e:
  #
  # If an exception was thrown it is assumed that there is no valid 
  # tolerance file present. Hence proceed as if there is no tolerance 
  # information.
  #
  pass

if process == diff:
  (lcs,ref_txt,Ntb,Nte,dat_txt,Mtb,Mte) = find_overall_lcs(lexact,ltol,lbest,chg_txt,ref_fnm,dat_fnm,diff_exe,feps,ieps,sys.stderr,separators,nguides,snake_trim,false)
  lcs = toldiff_lcs.filter_lcs(lcs,Ntb,Nte,Mtb,Mte,add_txt,del_txt)
  analysis = toldiff_diff.lcs_analysis(Ntb,Nte,Mtb,Mte,lcs,identical,equivalent,different)
  if output == output_full:
    (line_lcs,Nlb,Nle,Mlb,Mle) = toldiff_lcs.lcs_tokens2lines(lcs,ref_txt,Ntb,Nte,dat_txt,Mtb,Mte,nguides)
    toldiff_diff.lcs_to_diff(ref_txt,Nlb,Nle,dat_txt,Mlb,Mle,line_lcs,sys.stdout,sys.stderr,nguides)
  elif output == output_summary:
    if   analysis == identical:
      sys.stdout.write("%s" % text_identical)
    elif analysis == equivalent:
      sys.stdout.write("%s" % text_equivalent)
    elif analysis == different:
      sys.stdout.write("%s" % text_different)
    else:
      sys.stderr.write("illegal value of analysis")
  elif output == output_none:
    pass
  else:
    sys.stderr.write("illegal value of output")
  if   analysis == identical:
    sys.exit(exit_identical)
  elif analysis == equivalent:
    sys.exit(exit_equivalent)
  elif analysis == different:
    sys.exit(exit_different)
  else:
    sys.stderr.write("illegal value of analysis")
elif process == update:
  (lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines) = find_overall_lcs(true,true,true,chg_txt,ref_fnm,dat_fnm,diff_exe,feps,ieps,sys.stderr,separators,nguides,snake_trim,true)
  chg_txt = toldiff_update.lcs_to_change(lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines,chg_txt,feps,ieps,itol_scale,ftol_scale,ctol_scale)
  add_txt = toldiff_update.lcs_to_add(lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines,add_txt)
  del_txt = toldiff_update.lcs_to_delete(lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines,del_txt)
  store_tolerance(tol_fnm,chg_txt,add_txt,del_txt,sys.stderr,separators,nguides)
elif process == transfer:
  (lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines) = find_overall_lcs(true,true,false,chg_txt,ref_fnm,dat_fnm,diff_exe,feps,ieps,sys.stderr,separators,nguides,snake_trim,false)
  (chg_new,add_new,del_new) = toldiff_transfer.transfer_tol(lcs,Nb,ref_lines,Mb,dat_lines,chg_txt,add_txt,del_txt)
  store_tolerance(tol_new_fnm,chg_new,add_new,del_new,sys.stderr,separators_new,nguides)
elif process == show:
  Nb  = 1
  Ntb = 1
  (ref_txt,Ne,Nte) = load_file(ref_fnm,sys.stderr,separators,nguides)
  toldiff_show.show_tolerance(sys.stdout,ref_txt,Nb,Ne,chg_txt,add_txt,del_txt,sys.stderr,nguides)
else:
  try:
    sys.stderr.write("toldiff: internal error: invalid process")
  except IOError, e:
    pass
  sys.exit(999)
