import pyglib.math.matrix_basis as mb

m_basis = mb.hermitian_csc_matrix_basis(4)
for a in m_basis:
    print a.todense()
