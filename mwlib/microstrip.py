# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 10:36:55 2016

Microstrip class python 2.7

@author: sebastian
"""

import numpy as np
import scipy.constants as spc
from superconductor import *
from fcts_waves import *
from fcts_general import coth, csch

#==============================================================================
# Define general functions first
#==============================================================================

# conversion factor inch to cm

class Options(object):
    '''
    Helper class to be loaded in other classes to store options
    Options for specific applications are then accessed by self.Options.application
    '''
    def __init__(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)

    def change(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)

    def __getitem__(self, i):
        return self.__dict__[i]

class Microstrip(Superconductor):
    def __init__(self, w, h, epsr, f = 350e9, sc = None, opt = None, **kws):
        # Some control on options
        # Set default options if none are provided
        self.opt = Options(impedance = 'YW', epseff = 0)
        if opt != None:
            self.opt.change(**opt)
        self.options_ztot = {'YW' : self.Z_Yassin_Withington,
                             'P' : self.Z_Pozar}
        # TODO: NO SUPERCONDUCTOR => PEC functionality
        self.sc = sc
        if self.sc == None:
            self.PEC = True

        self.w = w
        self.h = h
        self.epsr = epsr
        # DEFAULT = 350 GHz
        self.f = f

        # RETURN PARAMETERS
        self.Ztot = None
        self.Cl = None
        self.Ll = None
        self.kinfrac = None
        self.gamma = 0 + 1j*0
#        self.kwargs_in(kws)

        # CALC ON INIT
        self.impedance(option = 'YW')

    @staticmethod
    def capacitance(w, h, t, epsr, mode = ''):
        return 0

#==============================================================================
#   Define dependent properties
#==============================================================================

    @property
    def eeff(self):
        if self.opt.epseff == 0:
            #
            _eeff = (self.epsr+1)/2 + (self.epsr-1)/2*(1+10*self.h/self.w)**(-0.555)
        elif self.opt.epseff == 1:
            _eeff = (self.epsr+1)/2 + (self.epsr-1)/2 / np.sqrt(1+12*self.h/self.w)
        else:
            print 'ERROR: INVALID OPTION FOR eeff GIVEN'
        return _eeff

    @property
    def eeff_sc(self):
        return self.beta2eeff()

    def beta2eeff(self):
        return beta2eps(self.beta, self.f)

    @property
    def alpha(self):
        return self.gamma.real

    @property
    def beta(self):
        return self.gamma.imag


#==============================================================================
#     Define Helper Functions
#==============================================================================
    def kwargs_in(self, **kwargs):
        for k,v in kwargs.items():
            if k in self.opt.__dict__.keys():
                setattr(self.opt, k, v)
            else:
                setattr(self, k, v)

    def loss_dielectric(self, tandelta, **kwargs):
        '''
        Dielectric loss according to Pozar in dB/mm
        '''
        self.kwargs_in(**kwargs)
        Npm = self.beta/2*self.epsr*(self.eeff-1)*tandelta/(self.epsr-1)
        return 20 / np.log(10) * Npm * 1e-3


    def impedance(self, **kwargs):
        '''
        check Microstrip.options_tot for possible implementations
        '''
        self.kwargs_in(**kwargs)
#        if option != None:
#            self.opt.impedance = option
        return self.options_ztot[self.opt.impedance]()

    def calc_impedance(self, **kwargs):
        self.Z_Yassin_Withington(**kwargs)
        return self

    def Z_Yassin_Withington(self, **kwargs):
        """
        Implementation of Yassin/Withington analytic model for superconducting microstrip. Compare Pieter Matlab code.
        Uses: f, w, h, t, Zs, epsr (only for efm calc)
        """
        # accept new input
        self.kwargs_in(**kwargs)
        # TODO:
        # check whether surface impedance Zs was already calculated for given set of parameters, otherwise recalculate

        # TODO: check whether sc exists => PEC
        self.PEC = False
        # get angular frequency
        omega = self.omega
        # Calculate effective dielectric constant
        efm = self.eeff
        # calculate b, p, eta, deltapn
        b  = 1 + self.sc.t/self.h
        p = 2*b**2 - 1 + 2*b*np.sqrt(b**2-1)
        eta = np.sqrt(p)*(spc.pi/2*self.w/self.h + (p+1)/(2*np.sqrt(p))*(1 + np.log(4/(p-1))) - 2*np.arctanh(1/np.sqrt(p)))
        if eta >= p:
            deltapn = eta
        else:
            deltapn = p
        # calculate ra, rb, Kf
        rb0 = eta + (p+1)/2 * np.log(deltapn)
        if self.w/self.h >= 5:
            rb = rb0
        else:
            rb = rb0 - np.sqrt(((rb0-1)*(rb0-p))) + (p+1)*np.arctanh(np.sqrt(((rb0-p)/(rb0-1)))) - 2*np.sqrt(p)*np.arctanh(np.sqrt(((rb0-p)/(p*(rb0-1))))) + spc.pi/2*self.w/self.h*np.sqrt(p)
        ra = np.exp(-1-spc.pi/2*self.w/self.h-(p+1)/np.sqrt(p)*np.arctanh(1/np.sqrt(p))-np.log((p-1)/(4*p)))
        Kf = self.h/self.w*2/spc.pi*np.log(2*rb/ra)
        # calculate ra1, ra2, is1, is2, rb1, rb2, ig1, ig2, chi
        ra1 = (1-ra) * (p-ra)
        is1 = np.log((2*p-(p+1)*ra + 2*np.sqrt(p*ra1))/(ra*(p-1)))
        rb1 = (rb-1) * (rb-p)
        is2 = - np.log(((p+1)*rb - 2*p - 2*np.sqrt(p*rb1 ))/(rb*(p-1)))
        rb2 = (rb+1) * (rb+p)
        ig1 = - np.log(((p+1)*rb+2*p+2*np.sqrt(p*rb2))/(rb*(p-1)))
        ra2 = (ra+1) * (ra+p)
        ig2 = np.log(((p+1)*ra+2*p+2*np.sqrt(p*ra2))/(ra*(p-1)))
        if self.w/self.h <= 2:
            chi = (is1+is2+ig1+ig2+spc.pi) / (2*np.log(rb/ra))
        else:
            chi = (is1+is2+ig1+ig2+spc.pi) / (2*np.log(2*rb/ra))
        # calculate characteristic impedance and some related quantities
        self.Ztot = 1./np.sqrt(efm)*self.h / (self.w*Kf*spc.epsilon_0*spc.c) * np.sqrt(1-2j*self.sc.Zs*chi*spc.epsilon_0*spc.c**2/(omega*self.h)).real # with superconductivity\
        print 1./np.sqrt(efm)*self.h / (self.w*Kf*spc.epsilon_0*spc.c) * np.sqrt(1-2j*self.sc.Zs*chi*spc.epsilon_0*spc.c**2/(omega*self.h))
        self.Z0 = 1./np.sqrt(efm)*self.h / (self.w*Kf*spc.epsilon_0*spc.c) # without superconductivity
        self.Z = 2*self.sc.Zs*chi/(self.w*Kf) + 1j*omega*spc.mu_0*self.h/(self.w*Kf)
        self.kinfrac = 1 - self.Z0**2/self.Ztot**2
        self.L = self.Z.imag/omega
        self.Lkpm = self.kinfrac*self.L
        self.Lg = (1-self.kinfrac)*self.L
        self.Y = 1j*(omega*spc.epsilon_0)*efm/self.h*(self.w*Kf)
        self.C = self.Y.imag/omega
        # calculate propagation constant
        self.gamma = 1j*omega * np.sqrt(efm)/spc.c * np.sqrt(1-2j*self.sc.Zs*chi*spc.epsilon_0*spc.c**2/(omega*self.h))

        # save variables in class; update self.properties to confirm new calculated values; return Ztot as fct-result
        self.Qimb = self.beta / (2*self.alpha)
        return self.Ztot

    def Z_Chang(self, sc2 = None, **kwargs):
        self.kwargs_in(**kwargs)
        self.PEC = False
        if sc2 == None:
            sc2 = self.sc

        def tanh1(x):
#            return 1/np.tanh(x)
            return np.arctan(x)

        # Recipe from paper: Chang (1979) - The inductance of a superconducting strip transmission line
        beta = 1 + self.sc.t/self.h
        p = 2*beta**2 - 1 + np.sqrt((2*beta**2 - 1)**2 - 1)
        eta = np.sqrt(p) * (np.pi*self.w/(2.*self.h) + (p+1)/(2.*np.sqrt(p))*(1+np.log(4./(p-1))) - 2*tanh1(np.sqrt(p)))
        delta = max(p, eta)
        rb0 = eta + (p+1)/2.*np.log(delta)
        lnra = -1 - np.pi*self.w/(2*self.h) - (p+1)/np.sqrt(p)*tanh1(1/np.sqrt(p)) - np.log((p-1)/(4*p))
        ra = np.exp(lnra)
        if self.w/self.h >= 5:
            rb = rb0
        else:
            rb = rb0 - np.sqrt((rb0-1)*(rb0-p)) + (p+1)*tanh1(np.sqrt((rb0-p)/(rb0-1))) - 2*np.sqrt(tanh1(np.sqrt((rb0-p)/(p*(rb0-1))))) + np.pi*self.w/(2*self.h)*np.sqrt(p)
        if self.w/self.h < 1:
            print 'Warning: impedance calculation outside of valid geometry'

        K = self.h/self.w * 2/np.pi * np.log(2*rb/ra)

        self.L = spc.mu_0/(self.w*K)*(self.h + self.sc.lambda0*(coth(self.sc.t/self.sc.lambda0) + 2*np.sqrt(p)/rb*csch(self.sc.t/self.sc.lambda0)) + sc2.lambda0*coth(sc2.t/sc2.lambda0))
        self.kinfrac = 1 - (spc.mu_0/(self.w*K)*(self.h)/self.L)
        self.C = self.eeff*spc.epsilon_0*self.w/self.h * K

        self.Z0 = np.sqrt(self.L / self.C)
        vph = 1/np.sqrt(self.L*self.C)
        self.gamma = 0 + 1j*vph2beta(vph, self.f)
        self.Ztot = self.Z0
        return self.Ztot

    def Z_Pozar(self, **kwargs):
        self.PEC = True
        self.kwargs_in(**kwargs)
        wh = self.w/self.h
        if wh <= 1:
            self.Z0 = 60/np.sqrt(self.eeff)*np.log(8/wh + wh/4)
        else:
            self.Z0 =120*spc.pi / (np.sqrt(self.eeff)*(wh + 1.393 + 0.667*np.log(wh+1.444)))
        self.Ztot = self.Z0
        return self.Ztot

    def inductance():
        return 0


if __name__ == '__main__':
    f = 350e9
    NbTiNstatic = Superconductor(14.8, 100e-8, f = f, t = 40e-9)
    NbTiN40 = Superconductor(14.8, 100e-8, f = f, t = 40e-9)
    NbTiN120 = Superconductor(14.6, 110e-8, f = f, t = 100e-9)
    NbTiN300 = Superconductor(15.4, 190e-8, f = f, t = 300e-9)
    Al50 = Superconductor(1.28, 1.450e-8, f = f, t = 50e-9)
    ms1 = Microstrip(1.5e-6, 0.25e-6, 10, f = f, sc = NbTiN120)
    ms2 = Microstrip(3e-6, 1e-6, 10, f = f, sc = NbTiN120)
    ms3 = Microstrip(1.5e-6, 1e-6, 10, f = f, sc = NbTiN120)

#    ms2 = Microstrip(1.4e-6, 1e-6, 10, f = f, sc = NbTiN40)
#    print ms1.Ztot, ms1.kinfrac
#    print ms1.Z_Chang(NbTiN300), ms1.kinfrac
##    print ms1.L, ms1.C
###    print ms1.Z_Pozar()
##    print NbTiN120.Ls
##    print NbTiN300.Ls
#
#    hlist = np.array([0.25])*1e-6
#    print "hsweep"
#    ms = ms2
#    for h in hlist:
#        ms.Z_Chang(NbTiN300, h = h)
#        print 'Ch', ms.h, ms.Ztot, np.sqrt(ms.beta)/(2*np.pi), ms.beta2eeff()
#        ms.Z_Yassin_Withington( h = h)
#        print 'YW', ms.h, ms.Ztot, np.sqrt(ms.beta)/(2*np.pi), ms.beta2eeff()
#
#    wlist = np.array([1.0, 1.4, 2.0, 3.0])*1e-6
#    print "wsweep"
#    ms = ms2
#    for w in wlist:
#        NbTiNstatic.update(t = 30e-9)
#        ms.Z_Chang(NbTiN300, w = w, sc = NbTiNstatic)
#        print 'Ch', ms.w, ms.Ztot, np.sqrt(ms.beta)/(2*np.pi), ms.beta2eeff()
#        ms.Z_Yassin_Withington( w = w, sc = NbTiNstatic)
#        print 'YW', ms.w, ms.Ztot, np.sqrt(ms.beta)/(2*np.pi), ms.beta2eeff()
#
#
#
#    tlist = np.array([40, 100])*1e-9
#    print "Ls sweep"
#    for t in tlist:
#        NbTiNstatic.update(t = t)
#        ms.sc = NbTiNstatic
#        ms.Z_Chang(NbTiN300, h = 0.25e-6)
#        print 'Ch', ms.h, ms.Ztot, np.sqrt(ms.beta)/(2*np.pi), ms.beta2eeff(), NbTiNstatic.Ls
#
