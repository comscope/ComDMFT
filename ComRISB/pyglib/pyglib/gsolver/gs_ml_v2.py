from __future__ import print_function
#Gutzwiller embedding Hamiltonian solver using machine learning.

import h5py, time, pickle, sys, numpy
import pyglib.gsolver as gsolver


def get_input_soc_only(imp=1, l=3):
    '''get the embedding Hamiltonian parameters of impurity with index imp,
    given the orbital angular momentum l
    and assuming local spin-orbit interaction dominates.
    '''
    j2 = 2*l
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        d1 = f['/D'][0, 0]
        d2 = f['/D'][j2, j2]
        e1 = f['/H1E'][0, 0]
        e2 = f['/H1E'][j2, j2]
        l1 = f['/LAMBDA'][0,0]
        l2 = f['/LAMBDA'][j2, j2]

    # to real
    phase1 = d1/numpy.sqrt(d1.real**2 + d1.imag**2)
    phase2 = d2/numpy.sqrt(d2.real**2 + d2.imag**2)
    # to positive
    if d1*numpy.conj(phase1).real > 0:
        phase1 *= -1.0
    if d2*numpy.conj(phase2).real > 0:
        phase2 *= -1.0

    d1 *= numpy.conj(phase1)
    d2 *= numpy.conj(phase2)

    delta = (e1+e2+l1+l2)/4.0
    delta1 = (e1-e2)/2.0
    delta2 = (l1-l2)/2.0
    eshift = l1*2*l + l2*(2*l + 2)
    x_list = numpy.array([[d1, d2, delta, delta1, delta2]])
    if numpy.max(numpy.abs(x_list.imag)) > 1.e-8:
        print(' warning: imaginary part detected in ml input!')
    x_list = x_list.real
    return x_list, eshift, [phase1, phase2]


def driver_gs_ml(l=3):
    '''driver for the machines learning based Gutzwiller solver.
    '''
    start_time = time.clock()

    imp = 1
    if '-i' in sys.argv:
        imp = int(sys.argv[sys.argv.index('-i')+1])
    if '-l' in sys.argv:
        l = int(sys.argv[sys.argv.index('-l')+1])
    if l == 3:
        subd = 'f_so_v2'
    else:
        raise ValueError(' l = {} not available!'.format(l))

    print(' solving emb. hamil. for impurity {}\n'.format(imp) +\
            ' with l = {} based on machine learning.'.format(l))

    dpath = gsolver.__path__[0]
    inp, eshft, phase = get_input_soc_only(imp=imp, l=l)
    with open('{}/{}/kr_fso.pkl'.format(dpath, subd), \
                'rb') as f:
        kr = pickle.load(f)
    res = kr.predict(inp)[0, :]

    # density matrix
    na2 = (2*l+1)*2
    na4 = na2*2
    dm = numpy.zeros([na4, na4], dtype=numpy.complex)
    for i in range(2*l):
        dm[i, i] = res[0]
        dm[i + na2, i + na2] = res[2]
        dm[i, i + na2] = res[4]*numpy.conj(phase[0])
        dm[i + na2, i] = numpy.conj(res[4])*phase[0]
    for i in range(2*l, na2):
        dm[i, i] = res[1]
        dm[i + na2, i + na2] = res[3]
        dm[i, i + na2] = res[5]*numpy.conj(phase[1])
        dm[i + na2, i] = numpy.conj(res[5])*phase[1]
    etot = kr.predict(inp)[0, 6]
    # convert to molecule energy
    etot -= eshft

    with h5py.File('EMBED_HAMIL_RES_{}.h5'.format(imp), 'w') as f:
        f['/DM'] = dm.T
        f['/emol'] = [etot.real]

    end_time = time.clock()

    print(' total time: {}'.format(end_time - start_time))
    sys.exit(0)



if __name__ == '__main__':
    driver_gs_ml()

