import numpy

'''Special models.
'''


class semicircular(object):
    '''semi-circular DOS.
    '''
    def __init__(self):
        '''define dos and cumulative dos function.
        '''
        self.dos = lambda e: 2./numpy.pi * numpy.sqrt(1-e**2)
        self.cdos = lambda e: (e*numpy.sqrt(1-e**2) \
                + numpy.arcsin(e)) / numpy.pi + 0.5

    def get_e_list_of_uniform_wt(self, nmesh=5000):
        '''Get the energy mesh with uniform weight.
        '''
        cdos_list = numpy.linspace(0,1,nmesh+1)
        from scipy.optimize import bisect
        e_list = [bisect(lambda x: self.cdos(x)-a, -1 ,1) \
                for a in cdos_list]
        e_list = numpy.asarray(e_list)
        e_list = (e_list[1:] + e_list[0:-1])/2
        return e_list
