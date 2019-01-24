import numpy


def primitive_vec2_reciprocal_primvec(avec):
    '''
    given primitive vector, return reciprocal primitive vector.
    '''
    bvec = []
    ac = numpy.cross(avec[1], avec[2])
    bvec.append(numpy.pi*2*ac/(avec[0].dot(ac)))
    ac = numpy.cross(avec[2], avec[0])
    bvec.append(numpy.pi*2*ac/(avec[1].dot(ac)))
    ac = numpy.cross(avec[0], avec[1])
    bvec.append(numpy.pi*2*ac/(avec[2].dot(ac)))
    return numpy.array(bvec)


def get_k_label(kpt, br_latt):
    '''
    given the bravais lattice br_latt and the k-point in fractions
    of primitive reciprocal lattice vectors, return the k-label.
    Ref. W. Setyawan, S. Curtarolo, comp. mater. sci. 49, 299.
    '''
    if br_latt == 'op':
        if numpy.allclose(kpt, [0., 0., 0.]):
            return 'GAMMA'
        elif numpy.allclose(kpt, [0.5, 0., 0.5]):
            return 'U'
        elif numpy.allclose(kpt, [0.5, 0.5, 0.5]):
            return 'R'
        elif numpy.allclose(kpt, [0.5, 0., 0.]):
            return 'X'
        elif numpy.allclose(kpt, [0.5, 0.5, 0.0]):
            return 'S'
        elif numpy.allclose(kpt, [0., 0.5, 0.]):
            return 'Y'
        elif numpy.allclose(kpt, [0., 0.5, 0.5]):
            return 'T'
        elif numpy.allclose(kpt, [0., 0., 0.5]):
            return 'Z'
        else:
            return ' '
    else:
        # to be implemented.
        return ' '
