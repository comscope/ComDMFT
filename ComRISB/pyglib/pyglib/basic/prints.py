def print_2d_array(dm, msg):
    print msg
    print 'real part'
    for row in dm:
        print ''.join('%9.4f' % (x.real) for x in row)
    print 'imag part'
    for row in dm:
        print ''.join('%9.4f' % (x.imag) for x in row)
