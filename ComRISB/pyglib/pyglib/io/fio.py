'''
Help functions for io.
'''
import os


def file_exists(fname):
    '''check whether a non-empty file exists.
    '''
    return os.path.isfile(fname) and os.path.getsize(fname) > 0
