# -*- coding: utf-8 -*-
"""
Created on Tue Sep 04 15:29:48 2018

@author: sebastian
"""

import numpy as np
import scipy.constants as spc
import scipy.interpolate as spi

@np.vectorize
def s13_2_Qc(s13, n = 1, length = 'quarterwave'):
    '''
    Input:
        s13 in magnitude form
    Output:
        Qc = pi/2/abs(s13)**2
    '''
    fac = {'quarterwave' : 1/2.,
           'halfwave' : 1.,}
    return n*fac[length]*spc.pi/np.abs(s13)**2

def gen_get_Lc(sonnetsweep, sweepparam = 'Lc', kid_length = 'quarterwave'):

    f_arr = sonnetsweep[0].f
    Qc_mat= np.vstack((s13_2_Qc(loc_ts.s31, length = kid_length) for loc_ts in sonnetsweep))
    Lc_arr = [v*1e-6 for v in sonnetsweep.sweep_values(sweepparam)]
    # define actual function
    @np.vectorize
    def get_Lc(Qc_design, F):
        """
        Input:
            Qc_design
            F [Hz]
        Output:
            Lc[m]
        """
        Lc_Q = [spi.splrep(np.log(Qc_mat[:,i]), Lc_arr, s = 0) for i in range(len(f_arr))]
        Lc_Q_eval = lambda x : [spi.splev(np.log(x), loc_Lc_Q, der = 0) for loc_Lc_Q in Lc_Q]
        Lc_Qdesign = Lc_Q_eval(Qc_design)
        Lc_F = spi.splrep(f_arr, Lc_Qdesign)
        return spi.splev(F, Lc_F, der = 0)
    return get_Lc