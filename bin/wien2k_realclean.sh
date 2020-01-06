#!/bin/bash
#
# Wien2K generates tonnes of files that if left lying around can
# get reused when running the code. Hence it may be difficult to
# establish how a result was obtained. This script deletes all
# files that Wien2K has generated to leave a completely clean 
# directory, ready to rerun a Wien2K calculation without 
# interference from artifacts from prior calculations.
#
path=`pwd`
case=`basename $path`
rm $case.b* $case.d* $case.i* $case.k* $case.n* $case.o* $case.r* $case.s* $case.v* *.def *.error
