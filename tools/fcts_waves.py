# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 17:15:56 2016

beta2vph(beta,f)
beta2lambda(beta)

f2w(f)

vph2eeff(vph)
@author: sebastian
"""

import numpy as np
import scipy.constants as spc
from scipy.optimize import *
import matplotlib.pyplot as plt

#==============================================================================
# General conversions
#==============================================================================

def beta2lambda(beta):
    return 2*spc.pi/beta

def beta2vph(beta, f):
    return f2w(f)/beta

def f2w(f):
    return 2 * spc.pi * f
    
def vph2eeff(vph):
    return (spc.c/vph)**2
    
def f_2_lambda(f, eeff):
    return spc.c/f/np.sqrt(eeff)
    
def lambda2_2_f(lambda2, eeff):
    return spc.c/2/lambda2/np.sqrt(eeff)

@np.vectorize
def magnitude_2_db(value):
    return 20*np.log10(value)    


#==============================================================================
#  3port total filter 
#==============================================================================
    
@np.vectorize
def s13_2_Qc(s13):
    '''
    Input:
        s13 in magnitude form
    Output:
        Qc = pi/2/abs(s13)**2
    '''
    return spc.pi/2./np.abs(s13)**2

def filter_s13(f, f0, Ql, Qi):
    temp = Ql / (np.sqrt(2 * Qi**2 * Ql**2 / (Qi - Ql)**2) * (1 + 2j*Ql *(f-f0)/f0))    
    s13sq = np.abs(temp)#**2
    return np.sqrt(s13sq)

def filter_s13_fit_QiFix(s13, f0, Ql, Qi, df = 3e9):
    mask = np.where((s13[0] < f0+df)*(s13[0] > f0-df))
    p0 = [f0, Ql]      
    popt, pocv = curve_fit(lambda x, a, b: filter_s13(x, a, b, Qi),  s13[0][mask], s13[1][mask], p0 = p0)
    return popt    
    
def filter_s13_fit(s13, f0, Ql, Qi, df = 3e9):
    mask = np.where((s13[0] < f0+df)*(s13[0] > f0-df))
    p0 = [f0, Ql, Qi]    
    popt, pocv = curve_fit(filter_s13, s13[0][mask], s13[1][mask], p0 = p0)
    return popt

def filter_s13_Qc(Ql, Qi):
    return Qi*Ql/(2*(Qi-Ql))
    
def filter_s13_print(s13fit):
    print 'F = ', s13fit[0]*1e-9, ' GHz; Qi =', s13fit[2], '; Qc = ', filter_s13_Qc(s13fit[1], s13fit[2]), '; Ql = ', s13fit[1]





#==============================================================================
#  2 port resonator
#==============================================================================
@np.vectorize
def halfwave_s12(f, f0, Ql, Qi):
    ff = (f-f0)/f0
    temp = (Ql/Qi + 2j*Ql*ff)/(1+2j*Ql*ff)
    return np.abs(temp)

def halfwave_s12_Qc(s12, f0, Qi, df = 5e9):
    imin = np.where(s12[1] == s12[1].min())
    Q = Qi * s12[1][imin]
    p0 = [f0, Q, Qi]
    mask = np.where((s12[0] < f0+df)*(s12[0] > f0-df))
    popt, pocv = curve_fit(halfwave_s12, s12[0][mask], s12[1][mask], p0 = p0)
    return popt[1]*popt[2] / (popt[2] - popt[1])
    

if __name__ == '__main__':
    
    print lambda2_2_f(68.64*1e-6, 35.5)*1e-9
    print lambda2_2_f(66.29*1e-6, 35.5)*1e-9
    print lambda2_2_f(63.92*1e-6, 35.5)*1e-9
    print f_2_lambda(369e9, 35.5)/2*1e6