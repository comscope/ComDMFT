# DMFT_MatDeLab bin

This directory contains executable binaries as well as executable scripts 
that are associated with the DMFT_MatDeLab project. While most of these
codes are meant to be available to regular code users some of them are
utilities for the build, test, and installation process.

- User scripts
  - cif2matdelab.py -- a script to generate a valid input file given a .cif file
  - matdelab_plot.py -- a script to plot band structures and density of states
- Developer utilities
  - build_CDMFTS_Mat_De_Lab_src_distro -- a bash script to generate a source code tar-ball as well as the web-pages
  - ini2qaini.py -- a script to transform a regular ini file into a QA ini file
  - test_perform.py -- a script to run and validate the test cases in the test directory
    - test_batch.py -- interacts with batch job systems
    - test_joblist.py -- collects test cases from multiple files and returns a list
  - toldiff.py -- a script to validate the test job outputs <https://sourceforge.net/projects/toldiff/>
    - toldiff_diff.py -- write a diff style output
    - toldiff_files.py -- handles the data files and tolerance files
    - toldiff_lcs.py -- compute the longest common sequence
    - toldiff_show.py -- write the reference file augmented with the tolerance data
    - toldiff_tokens.py -- tokenize data file content
    - toldiff_transfer.py -- transfer the tolerance data from an old reference file to a new one
    - toldiff_update.py -- update the tolerance data
    

