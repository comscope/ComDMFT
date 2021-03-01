#!/usr/bin/env python

import numpy as np

def dual_vector(matrix_in):
    matrix_out=np.zeros((3,3))
    vol=np.dot(np.cross(matrix_in[:,0], matrix_in[:,1]), matrix_in[:,2])
    matrix_out[:,0]=np.cross(matrix_in[:,1], matrix_in[:,2])/vol
    matrix_out[:,1]=np.cross(matrix_in[:,2], matrix_in[:,0])/vol
    matrix_out[:,2]=np.cross(matrix_in[:,0], matrix_in[:,1])/vol    
    return matrix_out

def vec_coeff(vector1, basis1, basis2):
    return np.dot(np.transpose(dual_vector(basis2)), np.dot(basis1, vector1))

def vec_norm(basis, vector):
    tempmat=np.zeros((3,3))
    for ii in range(3):
        for jj in range(3):
            tempmat[ii,jj]=np.dot(basis[:,ii], basis[:,jj])
    return np.dot(np.dot(vector, tempmat), vector)

