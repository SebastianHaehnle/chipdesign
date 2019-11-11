# -*- coding: utf-8 -*-
"""
Created on Thu May 03 11:13:15 2018

@author: sebastian
"""

import numpy as np
import scipy.constants as spc
from .touchstone import Touchstone, CST_Plotdata
from .abcd import ABCDmatrix
from ..util import find_nearest_id
import scipy.special
import os
from glob import glob
import matplotlib.pyplot as plt
import time

class FabryPerot(object):
    def __init__(self, F, ts_coupler, Z0_outside, Z0_line, length, eeff, alpha):
        self.F = F
        self.ts = ts_coupler
        self.Z0 = Z0_outside
        self.Z0_line = Z0_line
        self.length = length
        self.eeff = eeff
        self.alpha = alpha
        self.coupler = ABCDmatrix()
        
        self.load_coupler(self.ts)
        self.calc_line()
        
    def load_coupler(self, ts):
        self.ts = ts
        self.coupler.setMatrix_ts(self.ts, F = self.F)
    
    @property
    def f0(self):
        return spc.c/(2*self.length*np.sqrt(self.eeff))
        
    @classmethod
    def load_csv(cls, F, filename, part = 0):
        if not part:
            fef, value = np.loadtxt(filename, skiprows = 2, delimiter = ',', unpack = True)
        else:
            cur_part = 1
            got_numbers = False
            with open(filename, 'r') as ff:
                lines = ff.readlines()
                fef = []
                value = []
                for line in lines:
                    values = line.split(',')
                    if len(values) == 1:
                        if got_numbers:
                            cur_part += 1
                            got_numbers = False
                        continue
                    else:
                        got_numbers = True
                        if cur_part == part:
                            fef.append(float(values[0]))
                            value.append(float(values[1]))
                        elif cur_part < part:
                            continue
                        else:
                            break                            
        fef = np.array(fef)                   
        tck = scipy.interpolate.splrep(fef*1e9, value, k = 3)
        value_interp = scipy.interpolate.splev(F, tck, der=0)
        return value_interp
    
    @classmethod
    def spline_interp(cls, fi, fo, data):
        tck = scipy.interpolate.splrep(fi, data, k = 3)
        value_interp = scipy.interpolate.splev(fo, tck, der=0)
        return tck,value_interp
    
    def calc_line(self, z0_line = -1, length=-1, eeff=-1, alpha=-1):
        if length >=0:
            self.length = length
        if np.any(z0_line != -1):
            self.Z0_line = z0_line
        self.arr_n_hr = 2*self.length*np.sqrt(self.eeff)*self.F/spc.c
        arr_nint_hr = self.arr_n_hr.astype(int)
        arr_nint_u_hr = np.unique(arr_nint_hr)
        self.mask_n_hr = [np.where(arr_nint_hr==v)[0][0] for v in arr_nint_u_hr]
    
    def calc_Q(self, which = 'ci', squeezed = True, smooth = 0, exp = True):
        s21c = self.coupler.getS21(self.Z0, otype = 'linear')
        if exp:
            s21c = s21c*(1-np.abs(self.Z0-self.Z0_line)/(self.Z0+self.Z0_line))
            self.correction_s21c =(1-np.abs(self.Z0-self.Z0_line)/(self.Z0+self.Z0_line))
        else:
            self.correction_s21c = np.ones_like(self.Z0)
        self.s21c = s21c
        arr_Q_highres = self.arr_n_hr*np.pi/(s21c**2)
        if squeezed:
            delF_init = self.F/arr_Q_highres
            F2_init = self.F + delF_init/2 # Get 3dB freqs
            F1_init = self.F - delF_init/2
            phi2 = 2*np.pi*F2_init*np.sqrt(self.eeff)/spc.c * self.length # get phase at 3dB freqs for each center frequency -> this calculation is performed at the id corresponding to continuous peak centers
            phi1 = 2*np.pi*F1_init*np.sqrt(self.eeff)/spc.c * self.length
            # Get array with correct F,eeff->phi indexing  
            phi_correct = 2*np.pi*self.F*np.sqrt(self.eeff)/spc.c * self.length
            F2 = self.F[[find_nearest_id(phi_correct[self.mask_n_hr[ii]:self.mask_n_hr[ii+2]], phi2[v]) for ii,v in enumerate(self.mask_n_hr[1:-1])]]
            F1 = self.F[[find_nearest_id(phi_correct[self.mask_n_hr[ii]:self.mask_n_hr[ii+2]], phi1[v]) for ii,v in enumerate(self.mask_n_hr[1:-1])]]
            mask_Q_hr = self.mask_n_hr[1:-1]
            arr_Qcorr = self.F[mask_Q_hr]/(F2-F1)
            if smooth:
                arr_Qcorr = scipy.signal.savgol_filter(arr_Qcorr, smooth, 3)
            self.Q = arr_Qcorr
            self.mask_Q = mask_Q_hr
            self.smooth = smooth
            self.squeezed = squeezed
            return arr_Qcorr, mask_Q_hr
        else:
            self.Q = arr_Q_highres
            self.mask_Q = self.mask_n_hr
            self.smooth = smooth
            self.squeezed = squeezed
            return arr_Q_highres, self.mask_n_hr

    def get_Q_at_F(self, F_val, Q_arr = []):
        if Q_arr == []:
            Q_arr = self.Q
        assert F_val >= self.F[0] and F_val <= self.F[-1] # make sure that we are looking in the proper array
        return Q_arr[find_nearest_id(self.F[self.mask_Q], F_val)]

    def do_abcd(self):
        A = self.coupler
        B = ABCDmatrix()
        B.setMatrix_line(self.length, self.Z0_line, self.alpha, beta = -1, F = self.F, eps_eff = self.eeff)
        C = self.coupler
        self.A = A
        self.B = B
        self.C = C
        self.M = A.multiply(B.multiply(C))
        self.M.length = self.B.length
        return self.M
        
    
class SimPeaks(object):
    def __init__(self):
        self.f = []
        self.s21 = []
        self.peak_mask = []
        self.peak_ids = []
        self.peak_max = []
        self.peak_freq = []
        self.mask_3db = []
        self.mask_3db_left = []
        self.mask_3db_right = []
        self.peak_Q = []
        self.peak_phaseDelta = []
        self.peak_phaseLeft = []
        self.peak_phaseRight = []
        self.fit_masks = []
        
    def applyPrimary(self, func):
        self.peak_mask = func(self.peak_mask)
        self.peak_max = func(self.peak_max)
        self.peak_ids = func(self.peak_ids[0])
        self.peak_freq =func(self.peak_freq)
        self.mask_3db =func(self.mask_3db)
        self.mask_3db_left =func(self.mask_3db_left)
        self.mask_3db_right =func(self.mask_3db_right)
        
def get_peaks_sim(s21tot):
    peak_mask =  np.r_[False, s21tot[1:] > s21tot[:-1]] & np.r_[s21tot[:-1] > s21tot[1:], False]
    peak_ids = np.where(peak_mask)[0]
    return peak_mask, peak_ids
    
# def get_Q_sim(s21tot, f, peak_mask, peak_ids):

def get_Q_ofMatrix(matrix, Z0):
    s21tot = matrix.getS21(Z0)
    s21_linear = matrix.getS21(Z0, otype='linear')
    
    peak_mask, peak_ids = get_peaks_sim(s21tot)
    if type(peak_ids) == tuple:
        peak_ids = peak_ids[0]
    
    f = matrix.F
    peak_max = s21tot[peak_mask]
    peak_freq = f[peak_mask]
    num_peaks = len(peak_max)

    sp = SimPeaks()
    sp.peak_mask = peak_mask
    sp.peak_ids = peak_ids
    sp.f = f
    sp.s21= s21tot
    sp.peak_max = peak_max
    sp.peak_freq = peak_freq
    print 'masking'

    leftedge = lambda x: 0 if x == 0 else peak_ids[x-1]
    rightedge = lambda x: len(s21tot)-1 if x == num_peaks-1 else peak_ids[x+1]
    id_left = []
    id_right  = []
    Q_fitted = []
    for i,v in enumerate(peak_ids):
        jl = v-1
        jr = v+1
        left_id = None
        right_id  = None
        ref = peak_max[i]-3

        while not left_id or not right_id:
            if (s21tot[jl] <= ref or jl == leftedge(i)) and not left_id:
                left_id = jl
                id_left.append(jl)
            else:
                jl -= 1
            if (s21tot[jr] <= ref or jr == rightedge(i)) and not right_id:
                right_id = jr
                id_right.append(jr)
            else:
                jr += 1
        bw0 =  sp.f[right_id]-sp.f[left_id]
        p0 = [ np.pi*bw0/2 , sp.peak_freq[i], bw0]
        fit_mask = slice(left_id - 3*(v-left_id), right_id + 3*(right_id-v))
#         popt, pcov = scipy.optimize.curve_fit(cdt.fabry_perot.lorentz, sp.f[fit_mask], s21_linear[fit_mask], p0)
        sp.fit_masks.append(fit_mask)
#         Q_fitted.append(popt[1]/popt[2])
        if i > num_peaks-1:
            fig, ax = plt.subplots()
            ax.scatter(sp.f[fit_mask]*1e-9, s21_linear[fit_mask]**2, s = 5)
#             ax.plot(sp.f[fit_mask]*1e-9, cdt.fabry_perot.lorentz(sp.f[fit_mask], *popt), color = 'green')
            ax.plot(sp.f[fit_mask]*1e-9, lorentz(sp.f[fit_mask], *p0), color = 'orange')
            ax.axvline(f[id_right[-1]]*1e-9)
            ax.axvline(f[id_left[-1]]*1e-9)
            ax.set_title(str(i)+' '+str(v))
    matrix.eps_eff= np.full(matrix.F.shape, matrix.eps_eff)
#    sp.peak_n = 2*matrix.length*np.sqrt(matrix.eps_eff[peak_mask])*sp.peak_freq/spc.c
#    sp.peak_phaseDelta = (2*np.pi*f[id_right]*np.sqrt(matrix.eps_eff[id_right])/spc.c - 2*np.pi*f[id_left]*np.sqrt(matrix.eps_eff[id_left])/spc.c) * matrix.length
#    sp.peak_phaseLeft = 2*np.pi*f[id_left]*np.sqrt(matrix.eps_eff[id_left])/spc.c * matrix.length
#    sp.peak_phaseRight = 2*np.pi*f[id_right]*np.sqrt(matrix.eps_eff[id_right])/spc.c * matrix.length
    sp.id_left = id_left
    sp.id_right = id_right
    id_left = np.array(id_left)
    id_right = np.array(id_right)
    x1 = f[id_right]
    x2 = f[id_right-1]
    y1 = s21tot[id_right]
    y2 = s21tot[id_right+1]
    f_right = (((sp.peak_max-3-y1)/(y2-y1))*(x2-x1))+x1

    x1 = f[id_left]
    x2 = f[id_left+1]
    y1 = s21tot[id_left]
    y2 = s21tot[id_left+1]
    f_left = (((sp.peak_max-3-y1)/(y2-y1))*(x2-x1))+x1

    Q_abcd = sp.peak_freq/(f_right - f_left)
    sp.Q = Q_abcd
    return sp, Q_abcd    
        
def interpolate_s(fo, fi, s):
    """
    s input form: [freq, portA, portB], i.e. s21 = s[:,1, 0], s32 = s[:,2,1]

    returns interpolated complex s array
    """
    n_ports = s.shape[-1]

    output = np.zeros((len(fo), n_ports, n_ports), dtype = complex)

    for i in xrange(n_ports):
        for j in xrange(n_ports):
            mag = np.abs(s[:,i,j])
            angle = np.angle(s[:,i,j])
            new_mag = np.interp(fo, fi, mag)
            new_angle = np.interp(fo, fi, angle)
            output[:,i,j] = new_mag*np.exp(1j*new_angle)
            if any(output[:,i,j].real) == 0:
                print "what"

    return output


def abcd_fabryperot_inline(F, S, Z0cpw, Z0char, Z01, eps_eff, alpha, length_use, debug = True):
    # Set propagation constant and loss
    b = 2*np.pi*(F/spc.c)*np.sqrt(eps_eff) # propagation constant
    a = alpha*(F/(F[-1]+F[0])*2)**3 # alpha, set positive value in dB/mm for lossy case

    # Interpolate to desired frequency range
    # this part is only to make the following section easier writable/readable
    s11 = S[:,0,0]
    s21 = S[:,1,0]
    s12 = S[:,0,1]
    s22 = S[:,1,1]
    # transpose matrix,
    s21t = s12
    s12t = s21
    s11t = s11
    s22t = s22


    s21_tot = np.zeros_like(s11)

    # stupid debug section WARNING

    #%%
    M1t = np.array([[((1+s11)*(1-s22)+s12*s12)/(2*s12),
                         Z01*(((1+s11)*(1+s22)-(s12*s21))/(2*s21))],
                         [(1/Z01)*(((1-s11)*(1-s22)-(s12*s21))/(2*s21)),
                          (((1-s11)*(1+s22)+(s12*s21))/(2*s21))]])

    M2t = np.array([[np.cosh((a*length_use)+(1j*b*length_use)),
                         Z0cpw*np.sinh((a*length_use)+(1j*b*length_use))],
        [(1/Z0cpw)*np.sinh((a*length_use)+(1j*b*length_use)),
         np.cosh((a*length_use)+(1j*b*length_use))]])

    M3t = np.array([[((1+s11t)*(1-s22t)+s12t*s21t)/(2*s21t),
                         Z01*(((1.+s11t)*(1.+s22t)-(s12t*s21t))/(2.*s21t))],
        [(1/Z01)*(((1-s11t)*(1-s22t)-(s12t*s21t))/(2*s21t)),
         (((1-s11t)*(1+s22t)+(s12t*s21t))/(2*s21t))]])

    if debug:
        time_1 = time.time()
    Mtot = np.einsum('jki, kli -> ijl', np.einsum('jki, kli -> jli', M1t, M2t), M3t)
#    Mtot = np.array([M1t[:,:,i].dot(M2t[:,:,i]).dot(M3t[:,:,i]) for i in xrange(len(F))])
#    print Mtot == Mtest1
#    print Mtot.shape
#    Mtot = np.tensordot(M1t, M2t)

    s21_tot = np.array([2/(Mtot[i, 0,0] + (Mtot[i, 0,1]/Z0char[i]) + Mtot[i,1,0]*Z0char[i] + Mtot[i,1,1]) for i in xrange(len(F))])
    if debug:
       print "%.5f s for ABCD matrix multiplication" % (time.time() - time_1)
    return s21_tot

def abcd(touchstone_file, alpha, FP_length, Coup_ref_length, Freq_resolution = 0.001e9, **kwargs):
    """
    alpha in Nb/m
    Return [F, s21]
    """
    length_use = FP_length - 2*Coup_ref_length

    # Extract line parameters from sonnet simulation (touchstone file)
    ts = Touchstone(touchstone_file)
#    F = np.arange(ts.f[0], ts.f[-1] + Freq_resolution, Freq_resolution) # frequency range
    F = np.arange(300e9, 400e9 + Freq_resolution, Freq_resolution)
    Z0cpw = np.full(F.shape, np.mean(ts.fz0[:,1].real))
    Z0char = np.full(F.shape, np.mean(ts.fz0[:,0].real))
    Z01 = 50. # normalizing impedance for touchstone file
    assert Z01 == np.mean(ts.z0) # Check that touchstone file is normalized properly
    eps_eff = kwargs.pop('epseff', np.mean(ts.feps[:,1].real))

    # Set propagation constant and loss
    b = 2*np.pi*(F/spc.c)*np.sqrt(eps_eff) # propagation constant
    a = alpha*(F/(F[-1]+F[0])*2)**3 # alpha, set positive value in dB/mm for lossy case

    # Interpolate to desired frequency range
    s = interpolate_s(F, ts.f, ts.s)
    # this part is only to make the following section easier writable/readable
    s11 = s[:,0,0]
    s21 = s[:,1,0]
    s12 = s[:,0,1]
    s22 = s[:,1,1]
    # transpose matrix,
    s21t = s12
    s12t = s21
    s11t = s11
    s22t = s22


    s21_tot = np.zeros_like(s11)

    # stupid debug section WARNING
    #%%
    M1t = np.array([[((1+s11)*(1-s22)+s12*s12)/(2*s12),
                         Z01*(((1+s11)*(1+s22)-(s12*s21))/(2*s21))],
                         [(1/Z01)*(((1-s11)*(1-s22)-(s12*s21))/(2*s21)),
                          (((1-s11)*(1+s22)+(s12*s21))/(2*s21))]])

    M2t = np.array([[np.cosh((a*length_use)+(1j*b*length_use)),
                         Z0cpw*np.sinh((a*length_use)+(1j*b*length_use))],
        [(1/Z0cpw)*np.sinh((a*length_use)+(1j*b*length_use)),
         np.cosh((a*length_use)+(1j*b*length_use))]])

    M3t = np.array([[((1+s11t)*(1-s22t)+s12t*s21t)/(2*s21t),
                         Z01*(((1.+s11t)*(1.+s22t)-(s12t*s21t))/(2.*s21t))],
        [(1/Z01)*(((1-s11t)*(1-s22t)-(s12t*s21t))/(2*s21t)),
         (((1-s11t)*(1+s22t)+(s12t*s21t))/(2*s21t))]])

    Mtot = np.array([M1t[:,:,i].dot(M2t[:,:,i]).dot(M3t[:,:,i]) for i in xrange(len(F))])
    s21_tot = np.array([2/(Mtot[i, 0,0] + (Mtot[i, 0,1]/Z0char[i]) + Mtot[i,1,0]*Z0char[i] + Mtot[i,1,1]) for i in xrange(len(F))])
    lineparam = {'eps' : eps_eff,
                 'Zfp' : Z0cpw,
                 'Zchar' : Z0char,
                 'ts' : ts,}

    return F, s21_tot, lineparam

def get_peaks(s21_tot, F, bandwidth):
    """
    stupid peak finder, REQUIRES NO NOISE IN DATA
    """
    s21_max = np.max(np.abs(s21_tot))
    peak_range = int(bandwidth/(F[1]-F[0]))
    # Find peaks by checking local maxima and their relative height (the second check is to exluce peaks located at minima)
    peak_locs = np.where((np.abs(s21_tot[1:-1])> 0.01*s21_max)*(np.abs(s21_tot[1:-1]) > np.abs(s21_tot[0:-2]))*(np.abs(s21_tot[1:-1]) > np.abs(s21_tot[2:])))[0]
    bad_peaks = []
    for i, p in enumerate(peak_locs):
        if i > 0:
            comp_low = np.abs(s21_tot[p]) > 0.5*np.abs(s21_tot[peak_locs[i-1]])
#            print np.abs(s21_tot[p]),  np.abs(s21_tot[peak_locs[i-1]])
        else:
            comp_low = True
        if i < len(peak_locs)-1:
            comp_high = np.abs(s21_tot[p]) > 0.5*np.abs(s21_tot[peak_locs[i+1]])
#            print np.abs(s21_tot[p]), np.abs(s21_tot[peak_locs[i+1]])
        else:
            comp_high = True
#        print i, comp_high, comp_low
        if not(comp_high*comp_low):
            bad_peaks.append(i)
    print bad_peaks
    peak_locs = np.delete(peak_locs, bad_peaks)
    peak_freqs = F[peak_locs]
    peak_slices = [slice(max(loc-peak_range, 0), min(loc+peak_range, len(F)-1)) for loc in peak_locs]
    # Remove peaks with less than 3db bandwidth
    bad_peaks = []
    for i, p in enumerate(peak_slices):
        sl =  20*np.log10(np.abs(s21_tot[p]))
        if sl[0] >= max(sl)-3:
            bad_peaks.append(i)
            continue
        if sl[-1] >= max(sl)-3:
            bad_peaks.append(i)
    peak_locs = np.delete(peak_locs, bad_peaks)
    peak_freqs = np.delete(peak_freqs, bad_peaks)
    peak_slices = np.delete(peak_slices, bad_peaks)

    return peak_locs, peak_freqs, peak_slices, peak_range

def get_peaks_sim(s21tot):
    peak_mask =  np.r_[True, s21tot[1:] > s21tot[:-1]] & np.r_[s21tot[:-1] > s21tot[1:], True]
    peak_ids = np.where(peak_mask)
    return peak_mask, peak_ids
    
def get_Q_sim(s21tot, f, peak_mask, peak_ids):
    
    peak_max = s21tot[peak_mask]
    peak_freq = f[peak_mask]
    
    mask_3db = [np.where(s21tot < v-3)[0] for v in peak_max]
    right_3db = np.array([mask_3db[i][v:][0]+v for i, v in peak_ids])
    left_3db = np.array([mask_3db[i][:v][-1] for i, v in peak_ids])
    assert len(right_3db) == len(left_3db)
    
    x1 = f[right_3db]
    x2 = f[right_3db+1]
    y1 = s21tot[right_3db]
    y2 = s21tot[right_3db+1]
    f_right = (((peak_max-3-y1)/(y2-y1))*(x2-x1))+x1

    x1 = f[left_3db]
    x2 = f[left_3db-1]
    y1 = s21tot[left_3db]
    y2 = s21tot[left_3db+1]
    f_left = (((peak_max-3-y1)/(y2-y1))*(x2-x1))+x1

    Q = peak_freq/(f_right - f_left)
    
    return Q

def get_Q(s21_tot, F, peak_slices):
    """
    Get Q-value for all peaks in FP-transmission by evaluating the FWHM. The exact 3dB point is determined by linear interpolation for values close to the 3dB point. This results in precise Q-values for lower frequency resolution as long as the maximum of the peak is properly resolved.

    Lorentzian fit was tried, but did not work properly. TODO: search error at some point

    Return:
        Q : array of length len(peak_slices)
    """
    def Q_from_BW(fc, bw):
        return int(fc/bw)

    Q_out = np.zeros_like(peak_slices)
    peakmax_out = np.zeros_like(peak_slices)
    bw_points = np.zeros_like(peak_slices)
    for i, peak in enumerate(peak_slices):
        peak_range = (peak.stop-peak.start)/2

        s21_slice_db = 20*np.log10(np.abs(s21_tot[peak]))
        peakval = np.max(s21_slice_db)
        right_3db = np.where(s21_slice_db[peak_range:] <= peakval-3)[0][0:2]+peak_range
        left_3db = np.where(s21_slice_db[:peak_range] <= peakval-3)[0][-2:]
        bw_points[i] = right_3db[0] - left_3db[1]

        x1 = F[peak][right_3db[0]]
        x2 = F[peak][right_3db[1]]
        y1 = s21_slice_db[right_3db[0]]
        y2 = s21_slice_db[right_3db[1]]
        f_right = (((peakval-3-y1)/(y2-y1))*(x2-x1))+x1

        x1 = F[peak][left_3db[0]]
        x2 = F[peak][left_3db[1]]
        y1 = s21_slice_db[left_3db[0]]
        y2 = s21_slice_db[left_3db[1]]
        f_left = (((peakval-3-y1)/(y2-y1))*(x2-x1))+x1

        Q_out[i] = Q_from_BW(F[peak][np.argmax(s21_slice_db)], f_right-f_left)
        peakmax_out[i] = peakval
    return Q_out, 10**(peakmax_out/20), bw_points

def lorentz(x, A, fc, fwhm):
#    fwhm = fwhm*2
    return 2*A/np.pi*fwhm / ((fwhm)**2 + 4*(x - fc)**2)

#@np.vectorize
def frankel(f, w, s, eps_r, eps_eff = 0):
    """
    Returns (alpha [Np/m], eps_eff, Z0 [Ohm])
    """
    d = 350e-6
    k = s / (s + 2*w)
    eps_q = (eps_r + 1)/2.

    if eps_eff == 0:
        f_TE = spc.c/(4*d*np.sqrt(eps_r-1))
        q = np.log10(s/float(d))
        u = 0.54-0.649*q+0.015*(q**2)
        nu = 0.43-0.864*q+0.54*(q**2)

        a=10**((u*np.log(s/w))+nu)# % a parameter in effective permittivity
        # Calculate effective permittivity
        eps_eff = (np.sqrt(eps_q)+(np.sqrt(eps_r)+np.sqrt(eps_q))/(1+a*((f/f_TE)**(-1.8))))**2

    alpha = ((np.pi/2)**5)*2*(((1-(eps_eff/eps_r))**2)/np.sqrt(eps_eff/eps_r))*((((s+2*w)**2)*(eps_r**(3./2))*(f**3))/((spc.c**3)*scipy.special.ellipk((np.sqrt(1-k**1))**2)*scipy.special.ellipk(k**2)))

    beta = 2*np.pi*(f/spc.c)*np.sqrt(eps_eff)

    Z0CPW = ((120*np.pi)/np.sqrt(eps_eff))*(scipy.special.ellipk(np.sqrt(1-(k**2)))/(4.*scipy.special.ellipk(k)))
    return alpha, beta, eps_eff, Z0CPW

def s21min(mode, Qc):
    T = mode*np.pi/Qc

    return T**2/(2-T)**2


def s21max(Qc, Ql):
    return Ql/Qc

def s21_from_Qc_mode(mode, Qc):
    return np.sqrt(mode*np.pi/Qc)

def closest_index(array, value):
    if type(array) == list:
        array = np.array(array)
    return np.argmin(np.abs(array - value))
