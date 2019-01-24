#!/usr/bin/env python

from pyglib.run.comrisb import read_comrisb_ini, \
        initial_file_directory_setup, dft_risb, close_h_log


control,wan_hmat,imp = read_comrisb_ini()
initial_file_directory_setup(control)
dft_risb(control, wan_hmat, imp)
close_h_log(control)
