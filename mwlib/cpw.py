# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 14:58:57 2017

Holloway script, python edition

@author: sebastian
"""

import superconductor as sc
import fcts_general as fct

import numpy as np
import scipy.constants as spc
import scipy.special
from copy import deepcopy as copy


class CPW(object):
    def __init__(self, w, s, epsr, scgnd, PEC = False, scline = None, name = ''):
        """
        w = linewidth
        s = slotwidth
        """
        self.name = name
        self.w = float(w)
        self.s = float(s)
        self.epsr = float(epsr)
        self.scgnd = scgnd # THIS IS GND PLANE FOR HYBRID
        self.PEC = PEC
        if scline != None:
            self.scline = scline
        else:
            self.scline = scgnd # use same material for central line

    def kwargs_in(self, **kwargs):
        for k,v in kwargs.items():
            if k == 'f':
                self.scgnd.update(f = v)
                self.scline.update(f = v)
            else:
                setattr(self, k, v)

    @property
    def f(self):
        return self.scgnd.f

    @f.setter
    def f(self, f):
        self.scgnd.update(f = f)
        self.scline.update(f = f)

    @property
    def beta(self):
        try:
            return 2*np.pi*self.scgnd.f*np.sqrt(self.eeff_sc)/spc.c
        except:
            print "eeff_sc in CPW not defined, returning beta for (epsr+1)/2"

    @property
    def eeff(self):
        return (self.epsr + 1)/2.

    @property
    def lambda_sc(self):
        return spc.c/self.scgnd.f/np.sqrt(self.eeff_sc)

    def holloway(self, **kwargs):
        self.kwargs_in(**kwargs)
        Lsline = kwargs.pop('Lsline', self.scline.Ls)
        Lsgnd = kwargs.pop('Lsgnd', self.scgnd.Ls)
        if self.PEC:
            Lsline = 0
            Lsgnd = 0
        k1 = float(self.w / (self.w + 2*self.s))
        self.k1 = k1
        k2 = float(np.sqrt(1-k1**2))
        ellip1 = scipy.special.ellipk(k1**2)
        ellip2 = scipy.special.ellipk(k2**2)
        self.Lg = spc.mu_0/4. * ellip2/ellip1
        self.Cl = 4*spc.epsilon_0 * self.eeff * ellip1/ellip2

        gc = 1. / (4.*self.w*(1-k1**2)*ellip1**2)*(spc.pi + np.log(4.*spc.pi*self.w/self.scline.t) - k1*np.log((1+k1)/(1-k1)))
        self.gc = gc
        gg = k1 / (4.*self.w*(1-k1**2)*ellip1**2)*(spc.pi + np.log(4.*spc.pi*(self.w+2*self.s)/self.scgnd.t) - (1/k1)*np.log((1+k1)/(1-k1)))
        self.gg = gg
        self.Lkline = Lsline*gc
        self.Lkgnd = Lsgnd*gg

        self.linepartRatio = self.Lkline / (self.Lkline + self.Lkgnd)
        self.Lk = self.Lkline + self.Lkgnd
        self.alphak = self.Lk / (self.Lg + self.Lk)
        self.alphakline = self.Lkline / (self.Lg + self.Lk)
        self.alphakgnd = self.Lkgnd / (self.Lg + self.Lk)
        self.Z0=np.sqrt((self.Lg+self.Lk)/self.Cl)
        self.vphase = np.sqrt(1./((self.Lk+self.Lg)*self.Cl))
        self.eeff_sc = (spc.c/self.vphase)**2

        return self.alphak, self.eeff_sc, self.Z0

    def frankel(self, f, w, s, eps_r, eps_eff):
        d = 350e-6
        k = s / (s + 2*w)
        eps_q = (eps_r + 1)/2.

        alpha = ((np.pi/2)**5)*2*(((1-(eps_eff/eps_r))**2)/np.sqrt(eps_eff/eps_r))*((((s+2*w)**2)*(eps_r**(3./2))*(f**3))/((spc.c**3)*scipy.special.ellipk(np.sqrt(1-(k**2)))*scipy.special.ellipk(k)))

        beta = 2*np.pi*(f/spc.c)*np.sqrt(eps_eff)

        Z0CPW = ((120*np.pi)/np.sqrt(eps_eff))*(scipy.special.ellipk(np.sqrt(1-(k**2)))/(4.*scipy.special.ellipk(k)))
        return alpha, beta, Z0CPW


if __name__ == '__main__':
    f = 350e9
    rhoN = 100e-8
    shuttledllsAM = sc.shuttledLLS(t = 500e-9, f = f, rhoN = rhoN)
    shuttledllsMO = sc.shuttledLLS(t = 120e-9, f = f, rhoN = rhoN, Tc = 15.0)
    staticllsMO = sc.staticLLS(t = 300e-9, f = f)
    nordico = sc.Superconductor(14.6, 109e-8, f = f, t = 100e-9)
    al50 = sc.aluminum(t = 40e-9, f = f)
#    bondpad = CPW(600e-6, 225e-6, 11.44, shuttledllsMO)
#    bondpadAM = CPW(400e-6, 220e-6, 11.44, shuttledllsAM)
#    tlineAM = CPW(18e-6, 9e-6, 11.44, shuttledllsMO)
#    tlineMO = CPW(19e-6, 8e-6, 11.44, shuttledllsMO)
    tline = CPW(2e-6, 2e-6, 11.44, nordico)
#    kidwide = CPW(6e-6, 16e-6, 11.44, shuttledllsMO)
#    hybrid = CPW(3.2e-6, 1.4e-6, 11.44, shuttledllsMO, scline = al50)

    # Rainier paper
    al = sc.Superconductor(1.28, 1.1e-8, t = 80e-9,f = 6.6e9)
    nbtin = sc.Superconductor(14.5, 1.3e-6, t = 300e-9, f = 350e9)
    nbtin_thin = sc.Superconductor(15, 1.2e-6, t = 100e-9, f = 350e9)
    hybrid = CPW(3e-6, 3e-6, 11.44, scgnd = nbtin, scline = al)
    wide = CPW(5.4e-6, 23.7e-6, 11.44, scgnd = nbtin, scline = nbtin)
    hybrid.holloway()
    wide.holloway()

#    lwide = 4e-3
#    lhyb = 1e-3
#    foo = lhyb*hybrid.Lk /(lhyb*(hybrid.Lg + hybrid.Lk) + lwide*(wide.Lg + wide.Lk))

    # pieter paper example
    al = sc.Superconductor(1.28, 2.2e-8, t = 50e-9, f = 6.6e9)
    algnd = sc.Superconductor(1.28, 0.28e-8, t = 100e-9, f = 6.6e9)
#    nbtin = sc.Superconductor()
    cpw = CPW(1e-6, 1e-6, 10.5, scline = nbtin_thin, scgnd = nbtin_thin)
    cpw.holloway()
#    print tline.holloway()

#    # Measured LT149
#    badNbTiN = sc.Superconductor(11.55, 157e-8, t = 100e-9, f = 350e9)
#    w = 5e-6
#    s = w
#    cpw = CPW(w, s, 10.3, nbtin_thin)
#    cpw_overetch = CPW(w-0.1e-6, s+0.1e-6, 10.3, nbtin_thin)
#    print nbtin_thin.Ls
#    print cpw.holloway()
#    print cpw_overetch.holloway()
