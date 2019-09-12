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

def beta2eps(beta, f):
    return(vph2eeff(beta2vph(beta, f)))

def eps2vph(eps):
    return spc.c/np.sqrt(eps)


def eps2beta(eps, f):
    return 2*np.pi*f*np.sqrt(eps)/spc.c

def f2w(f):
    return 2 * spc.pi * f

def vph2eeff(vph):
    return (spc.c/vph)**2

def vph2beta(vph, f):
    return f2w(f)/vph

def f_2_lambda(f, eeff):
    '''
    LEGACY FUNCTION, USE f2lambda INSTEAD
    '''
    return spc.c/f/np.sqrt(eeff)

def f2lambda(f, eeff):
    return spc.c/f/np.sqrt(eeff)

def lambda_2_f(lamb, eeff):
    '''
    LEGACY FUNCTION USE lambda2f INSTEAD
    '''
    return spc.c/lamb/np.sqrt(eeff)

def lambda2f(lamb, eeff):
    return spc.c/lamb/np.sqrt(eeff)


@np.vectorize
def mag2db(value):
    """
    Deprecated, use mag2db20
    """
    value = float(value)
    return 20*np.log10(value)

@np.vectorize
def mag2db1(value):
    """
    Deprecated, use mag2db10
    """
    value = float(value)
    return 10*np.log10(value)

@np.vectorize
def mag2db10(value):
    value = float(value)
    return 10*np.log10(value)

@np.vectorize
def mag2db20(value):
    value = float(value)
    return 20*np.log10(value)



@np.vectorize
def db2mag(value):
    value = float(value)
    return 10**(value/20.0)


@np.vectorize
def db2mag1(value):
    value = float(value)
    return 10**(value/10.0)

def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return float('.'.join([i, (d+'0'*n)[:n]]))

#==============================================================================
# general microwaves
#==============================================================================

def reflection(Z0, Zl):
    Z0 = float(Z0)
    Zl = float(Zl)
    return abs(Z0-Zl)/(Z0+Zl)

#==============================================================================
#  3port total filter
#==============================================================================

def filter_s13(f, f0, Ql, Qi):
    temp = Ql / (np.sqrt(2 * Qi**2 * Ql**2 / (Qi - Ql)**2) * (1 + 2j*Ql *(f-f0)/f0))
    s13sq = np.abs(temp)#**2
    return s13sq

def filter_s13_Qc(s13, f0, Ql, Qi, df = 5e9):
    imin = np.where(s13[1] == s13[1].min())
    mask = np.where((s12[0] < f0+df)*(s12[0] > f0-df))
    p0 = [f0, Ql, Qi]
    popt, pocv = curve_fit(filter_s13, s13[0][mask], s13[1][mask], p0 = p0)
    return popt

def s13max_2_Qi(s13max, Qc):
    return Qc / (np.sqrt(2)/s13max - 2)



#==============================================================================
#  2 port resonator
#==============================================================================

def s21min_2_Qi(s21db, Qc):
    s21mag = db2mag1(s21db)
    return Qc*(1/s21mag - 1)

def Qi_2_s21min(Qi, Qc):
    s21mag = float(Qc)/(Qi+Qc)
    return mag2db1(s21mag)

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


def Qloaded(Qarray):
    Qarray = np.array(Qarray)
    return 1./np.sum(1./Qarray)

def QlQi_2_Qc(Ql, Qi):
    return Ql*Qi / (Qi - Ql)


#==============================================================================
# 3 port coupler
#==============================================================================

@np.vectorize
def s13mag_2_Qc(s13):
    '''
    Input:
        s13 in magnitude form
    Output:
        Qc = pi/2/abs(s13)**2
    '''
    return spc.pi/2./np.abs(s13)**2

@np.vectorize
def s13db_2_Qc(s13):
    s13 = np.abs(db2mag1(s13))
    return spc.pi/2./s13

if __name__ == '__main__':
    s13old = -47.56 #db
#    s13old = -46.28
    s13new = -44.66 #db
    Qold = s13db_2_Qc(s13old)
    Qnew = s13db_2_Qc(s13new)
    print Qold, Qnew, Qnew/Qold

    #6.5 GHz comparison
    s13 = -46 #db
    Qs13 = s13db_2_Qc(s13)
    Qfullkid = 6.578875 / (abs(6.578827 - 6.578922))
#    Qfullkid = 6.577942 / (6.577869 - 6.5778014)
    print Qs13, Qfullkid
#    s11db = -0.0112 # db
##    s11 = db2mag(s11)
#    s21 = (1-db2mag1(s11db)) *0.09#mag
#    s21db = mag2db1(s21)
#    Qkid = s13db_2_Qc(s21db)
#    print s21db, Qkid