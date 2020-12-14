# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 16:00:33 2016

@author: sebastian
"""

import numpy as np
import scipy as sp
import scipy.constants as spc
import scipy.integrate
import matplotlib.pyplot as plt

class Superconductor(object):
    def __init__(self, Tc, rhoN, name = 'default', T = 0.001, f = 350e9, N0 = None, **kwargs):
        # INIT PARAMETERS
        self.Tc = Tc
        self.rhoN = rhoN
        self.name = name
        # default temperature is 300 mK
        self.T = T
        # default frequency is 350GHz
        self.f = f
        # Thickness of superconductor. TODO: if None assume t >> lambda
        self._t = None
        # Single spin density of states, default is Al value
        self.use_deltaT = kwargs.pop('use_deltaT', False)
        self.Tdb = kwargs.pop('Tdb', 0)
        self.N0 = N0 or (1.7e10 * 1e6**3 / spc.e)
        self.etapb = kwargs.pop('etapb', 0.57)

        # RETURN PARAMETERS
        self._lambda0 = None
        self.sigma1 = None
        self.sigma2 = None
        self.Zs = None
        # Set eventual additional stuff
        for k,v in kwargs.items():
            setattr(self, k, v)
        # RUN IMMEDIATELY ON INIT
        self.update()
#==============================================================================
#   Define dependent properties
#==============================================================================
    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, t):
        self._t = t

    @property
    def sigmaN(self):
        return 1/self.rhoN

    @property
    def delta(self):
        """
        delta in Joule
        """
        return 1.76 * spc.k * self.Tc

    @delta.setter
    def delta(self, delta):
        """
        delta needs to be input in Joule
        """
        self.Tc = delta / (1.76 * spc.k)


    @property
    def deltaGHz(self):
        """
        delta in GHz
        """
        return 1.76*spc.k*self.Tc / spc.h *1e-9

    def solve_deltaT_lowlvl(self, N0V, Tdb, Delta, T, wtdaccur):
        """
        After matlab code from jochem (energygap 5). Is used inside calc_deltaT for actual solution. calc_deltaT has cleaner input
        N0 is the density of states at the fermi surface/10^23/eV/cm^3
        V is the scattering parameter /10^23eVcm^3
        %scrip needs N0V=N0*V
        %Tdb is the Debeye temperature in K
        %Delta is guessed delta and should be accurate within a factor of 2 and is a guess for
        %T=0 in meV should be nonzero!
        %Tcalc is the temperature at which you want to calculate the gap
        %wtdaccur is the accuracy used for calculating the integral and for the
        %convergence of iteration.
        %
        %OUTPUT
        %calculatedDelta is a single value in J
        """
        NVinv0 = float(1/(N0V))
        Delta1 = Delta
        algconst = 0.99
        k = 0
        
        diff3 = 1
        diffr = 1
        
        
        if Delta1 == 0:
            return 0
        else:
            wtdaccur = wtdaccur*NVinv0
            while diffr > wtdaccur:
                k = k+1
                integral_fun = lambda eps: np.tanh(np.sqrt(eps**2+Delta1**2)/(2*spc.k*T))/(2*np.sqrt(eps**2+Delta1**2))
                NVinv = sp.integrate.quad(integral_fun, -Tdb*spc.k, Tdb*spc.k, epsabs = NVinv0*wtdaccur )[0]
                diff1 = NVinv - NVinv0
                
                Delta2 = Delta1*(1+1e-6)
                integral_fun = lambda eps: np.tanh(np.sqrt(eps**2+Delta2**2)/(2*spc.k*T))/(2*np.sqrt(eps**2+Delta2**2))
                NVinv = sp.integrate.quad(integral_fun, -Tdb*spc.k, Tdb*spc.k, epsabs = NVinv0*wtdaccur )[0]
                diff2 = NVinv - NVinv0           
                deriv = (diff2-diff1)/(Delta2-Delta1)
                # calc new guess
                Deltanew = Delta1 - algconst*diff1/deriv
                diffr = abs(diff1)
                Delta1 = Deltanew
        if np.isnan(Delta1):
            return 0
        else:
            return Delta1
        

    def calc_deltaT(self, T, Tdb, **kwargs):
        """
        Tdb is debye temperature of the material
        T is temperature, can be array
        """
        T = float(T)
        Tdb = float(Tdb)
        # Calculate correct value of N0V
        Delta0 = self.delta*spc.electron_volt*1e3 # delta0 in meV
        N0p1 = np.log((spc.k*Tdb + np.sqrt((spc.k*Tdb)**2-self.delta**2))/self.delta)
        
        N0V = 1/N0p1
        self.N0V = N0V
        T = np.atleast_1d(T)
        Delta = np.zeros_like(T)
        for ii, iT in enumerate(T):
            Delta[ii] = self.solve_deltaT_lowlvl(N0V, Tdb, self.delta, T[ii], 1e-12)
        return Delta
        

    @property
    def omega(self):
        return 2 * spc.pi * self.f

    @property
    def lambda0(self):
        if self._lambda0 == None:
            self.pendepth()
        return self._lambda0

    @property
    def Ls(self):
        return self.Zs.imag / self.omega

    @property
    def Rs(self):
        return self.Zs.real

    def Lsf(self, f):
        self.update(f = f)
        return self.Ls

    def Qi(self, alpha_k, **kwargs):
        self.kwargs_in(**kwargs)
        self.update()
        assert self.t != None
        beta = 2 * self.t / self.lambda0 / np.sinh(2*self.t / self.lambda0)
        Qi = 2* self.sigma2 / self.sigma1 / alpha_k / beta
        return Qi


    def nqp(self, unit = 'm', **kwargs):
        '''
        Theoretical Quasiparticle density, multiply by volume to get absolute
        Can change initial parameters if needed:
            T :: SC temperature
            N0 :: Single spin density
        '''
        self.kwargs_in(**kwargs)
        nqp = 2*self.N0*np.sqrt(2*spc.pi*spc.k*self.T*self.delta)*np.exp(-self.delta/(spc.k*self.T))
        if unit == 'um':
            nqp *= 1e-6**3
        return nqp

    def infodump(self):
        print '\nSuperconductor impedance properties'
        print 'name: ', self.name, '\nf = ', self.f*1e-9 , 'GHz; Tc = ', \
        self.Tc, 'K; rhoN = ', self.rhoN , 'Ohmm; t = ', self.t*1e9, 'nm'
        print 'Zs = ', self.Zs \
        , '\nLs = ', self.Ls, 'pH/Sq\n'
#==============================================================================
#   Define Helper functions
#==============================================================================
    def update(self, **kwargs):
        self.kwargs_in(**kwargs)
        if self.use_deltaT:
            if not self.Tdb:
                print "NO DEBYE TEMPERATURE SET, ABORTING DELTA(T) CALCULATION"
                
        self.sigma12()
        self.pendepth()
        self.surfaceimpedance()
        return self

    def kwargs_in(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)

    def dist_fd(self, E, mu = 0):
        """
        Fermi-dirac distribution
        """
        return 1 / ( np.exp((E - mu) / (spc.k * self.T)) + 1)

    def dist_be(self, E, mu = 0):
        """
        Bose-Einstein distribution
        """
        return 1 / ( np.exp((E - mu) / (spc.k * self.T)) - 1)

#==============================================================================
#   Define Major functions
#==============================================================================
    def pendepth(self, thermalsigma = False, **kwargs):
        if thermalsigma:
            
            self._lambda0 =  1/np.sqrt(spc.mu_0 * self.omega * self.sigma2_thermal())
        else:
            self._lambda0 =  1/np.sqrt(spc.mu_0 * self.omega * self.sigma2)
        return self._lambda0

    @property
    def lambdaApprox(self):
        return 105*1e-9* np.sqrt(self.rhoN*1e8/self.Tc)

    def Zs_thickfilm(self):
       return np.sqrt(1j*spc.mu_0*self.f*2*np.pi/(self.sigma1-1j*self.sigma2))

    def Ls_thickfilm(self):
        return np.sqrt(spc.mu_0/(self.f*np.pi*2*self.sigma2))

    def sigma1_thermal(self, **kwargs):
        self.kwargs_in(**kwargs)
        hf = spc.h*self.f
        kbt = spc.k*self.T
        a = 4*self.delta/(hf)
        b = np.exp(-self.delta/kbt)
        c = np.sinh(hf/(2*kbt))*sp.special.k0(hf/(2*kbt))
        return self.sigmaN*a*b*c

    def sigma2_thermal(self, **kwargs):
        self.kwargs_in(**kwargs)
        hf = spc.h*self.f
        kbt = spc.k*self.T
        a = np.exp(-self.delta/kbt)
        b = np.exp(-hf/(2*kbt))
        c = sp.special.i0(hf/2*kbt)
        return self.sigmaN*np.pi*self.delta/hf*(1-2*a*b*c)

    def sigma12(self, **kwargs):
        self.kwargs_in(**kwargs)
        if self.use_deltaT and self.Tdb:
            delta = self.calc_deltaT(self.T, self.Tdb)
        else:
            delta = self.delta
        hf = spc.h * self.f
        #define functions for integration
        funcg1 = lambda E : (E**2 + delta**2 + hf*E) / ((np.sqrt(E**2-delta**2)) * np.sqrt((E+hf)**2 - delta**2))
        funcg2 = lambda E : (E**2 + delta**2 + hf*E) / ((np.sqrt(delta**2-E**2)) * np.sqrt((E+hf)**2 - delta**2))
        integr1a = lambda E : (self.dist_fd(E) - self.dist_fd(E + hf)) * funcg1(E)
        integr1b = lambda E : (1 - 2*self.dist_fd(E+hf)) * funcg1(E)
        integr2 = lambda E : (1 - 2*self.dist_fd(E+hf)) * funcg2(E)

        # integrate first part of sigma1
        intlim = [delta, 1e-20, 1e-18, 1e-16, 1e-12, 1e-10, 1e-5, 1e-2, 1e-1, 1e-0]  #
        sigma1aint = np.array([sp.integrate.quad(integr1a, intlim[i], intlim[i+1], epsabs = 1e-70)[0] for i in range(len(intlim)-1)])
        sigma1a = sigma1aint.sum()

        # integrate second part of sigma1
        if hf >= 2*delta :
            sigma1b = sp.integrate.quad(integr1b, delta - hf,  -delta, epsabs = 1e-70)[0]
        else:
            sigma1b = 0
            
        sigma1 = 2 / hf * sigma1a + 1 / hf * sigma1b

        # integrate sigma2
        if hf >= 2*delta :
            sigma2 = sp.integrate.quad(integr2, -delta, delta, epsabs = 1e-70)[0]
        else:
            sigma2 = sp.integrate.quad(integr2, delta - hf, delta, epsabs = 1e-70)[0]
        sigma2 = 1 / hf * sigma2

        self.sigma1 = sigma1*self.sigmaN
        self.sigma2 = sigma2*self.sigmaN
        return (self.sigma1, self.sigma2)

    def surfaceimpedance(self, **kwargs):
        """
        Pieter thesis eq 2.20
        """
        self.kwargs_in(**kwargs)
        if self.t == None:
            self.t = 10*self.lambda0
        sigma = self.sigma1 - 1j*self.sigma2
        self.Zs = np.sqrt(1j*spc.mu_0 * self.omega / sigma) / np.tanh(np.sqrt(1j*self.omega * spc.mu_0 * sigma)*self.t)
#        self.update_prop()
        return self.Zs

def staticLLS(t, f, Tc = 14.8, T = 0.12, rhoN = 100e-8):
    return Superconductor(Tc = Tc, rhoN =  rhoN, T = T, f = f, t = t)

def shuttledLLS(t, f, Tc = 15.4, rhoN = 190e-8, T = 0.12):
    return Superconductor(Tc = Tc, rhoN = rhoN, T = T, f = f, t = t)

def nordico(t, f, Tc = 15.0, T = 0.12, rhoN = 120e-8):
    return Superconductor(Tc = Tc, rhoN = rhoN, T = T, f = f, t = t)

def aluminum(t, f, Tc = 1.28, rhoN = 1.45e-8, T = 0.12, **kwargs):
    return Superconductor(Tc = Tc, rhoN = rhoN, T = T, f = f, t = t, **kwargs)




if __name__ == '__main__':
    f_range = np.arange(300, 1200, 20)*1e9
#    nbtin = nordico(100e-9, f_range[0])
    # output array
    Ls_range = np.zeros_like(f_range)

    Tc= 15
    foo = Superconductor(15.0, 104e-8, f = 350e9, T = 0.12, t = 150e-9)

    Rarr = []
    farr = []
    Larr = []
    sig1arr = []
    sigTarr = []
    al = aluminum(40e-9, 5e9, use_delta = True, Tdb = 420)
    
#    print al.delta, al.calc_deltaT(1.2, 420)
    
#    for fi in np.arange(300e9, 1000e9, 10e9):
#        foo.update(f = fi)
#        farr.append(fi)
#        Larr.append(foo.Ls)
#        Rarr.append(foo.Rs)
#        sig1arr.append(foo.sigma1)
#        sigTarr.append(foo.sigma1_thermal())
#    foo.update(f = 350e9)
#
#    fig,ax = plt.subplots()
#    ax.set_yscale('log')
#    ax.plot(farr, sig1arr)
#    ax.plot(farr, sigTarr)
#    ax.plot(farr, Rarr)

#    fig.savefig('C:\Users\sebastian\ownCloud\sigma1.png', dpi = 400)
#    for i, f in enumerate(f_range):
#        if f >= nbtin.deltaGHz*2e9: continue
#
#        Ls_range[i] = nbtin.update(f = f).Ls*1e12
#        print f*1e-9, Ls_range[i]
#    plt.plot(f_range, Ls_range)
#
#    # write to file:
#    filepath = r'C:\Users\sebastian\ownCloud\Chips\FP design\Wideband CPW\Coupler\nbtin_Ls_table.csv'
#    np.savetxt(filepath, np.array([f_range, Ls_range]).transpose(), delimiter = ',', header = 'Tc = 15.0K, rhoN = 120uOhmcm, t = 100nm', comments = '!')