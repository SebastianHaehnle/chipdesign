# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 11:41:20 2019

@author: sebastian
"""

import numpy as np
from .fabry_perot import interpolate_s, spc

class ABCDmatrix(object):
    def __init__(self, Fresolution = -1):
        self.matrix = np.zeros((2,2))
        self.F = np.array([0])
        self.Fresolution = Fresolution
        self.eps_eff = None

    def getS21(self, Z0char, otype = 'db'):
        Z0char = np.atleast_1d(Z0char)
        if len(Z0char) != len(self.F):
            Z0char = np.full(len(self.F), np.mean(Z0char))
        fct = {'db' : lambda x : 20*np.log10(np.abs(x)),
               'linear' : lambda x : np.abs(x),
               'imag' : lambda x : np.imag(x),
               'real' : lambda x : np.real(x),
               'complex' : lambda x : x,
               }
        return fct[otype](np.array([2/(self.matrix[0,0,i] + (self.matrix[0,1,i]/Z0char[i]) + self.matrix[1,0,i]*Z0char[i] + self.matrix[1,1,i]) for i in xrange(len(self.F))]))


    def getS11(self, Z0char, otype = 'db'):
        Z0char = np.atleast_1d(Z0char)
        if len(Z0char) != len(self.F):
            Z0char = np.full(len(self.F), np.mean(Z0char))
        fct = {'db' : lambda x : 20*np.log10(np.abs(x)),
               'linear' : lambda x : np.abs(x),
               'imag' : lambda x : np.imag(x),
               'real' : lambda x : np.real(x),
               'complex' : lambda x : x,
               }
        return fct[otype](np.array([(self.matrix[0,0,i] + (self.matrix[0,1,i]/Z0char[i]) - self.matrix[1,0,i]*Z0char[i] - self.matrix[1,1,i])/(self.matrix[0,0,i] + (self.matrix[0,1,i]/Z0char[i]) + self.matrix[1,0,i]*Z0char[i] + self.matrix[1,1,i]) for i in xrange(len(self.F))]))

    def setMatrix_C(self, C, F = []):
        self.F = np.atleast_1d(F)
        if self.F[1] - self.F[0] != self.Fresolution:
            self.F = np.arange(self.F[0], self.F[-1]+self.Fresolution, self.Fresolution)
        s11 = np.ones_like(self.F)
        s22 = np.ones_like(self.F)
        s12 = 1/(1j*2*np.pi*self.F*C)
        s21 = np.zeros_like(self.F)
        self.matrix = np.array([[s11, s12],[s21, s22]])


    def setMatrix_ts(self, ts,  ctype = 'sparam', inverse = False, length = 0, alpha = 0, eps_force = 0, F = [], Z0 = 0):
        """
        alpha in [dB/mm]
        """
        if len(F) != 0:
            s = interpolate_s(F, ts.f, ts.s)
            self.F = np.array(F)
        else:
            s = ts.s
            self.F = np.arange(ts.f[0], ts.f[-1] + self.Fresolution, self.Fresolution)
            s = interpolate_s(self.F, ts.f, ts.s)

        if ctype == 'sparam':
            self.Z0 = Z0
            self.setMatrix_Smatrix(np.moveaxis(s, 0, -1), np.mean(ts.z0), inverse)
        elif ctype == 'line':
            if eps_force > 0:
                eps_eff = eps_force
            elif eps_force < 0:
                eps_eff = np.interp(self.F, ts.f, ts.feps[:,1].real)
            else:
                eps_eff = np.mean(ts.feps[:,1].real)
            b = 2*np.pi*(self.F/spc.c)*np.sqrt(eps_eff) # propagation constant
            a = alpha*(self.F/(self.F[-1]+self.F[0])*2)**3 # alpha, set positive value in dB/mm for lossy case
            Z0 = np.full(self.F.shape, np.mean(ts.fz0[:,1].real))
            self.b = b
            self.eps_eff = eps_eff
            self.b = b
            self.a = a
            self.Z0 = Z0
            self.setMatrix_line(length, Z0, a, b)



    def setMatrix_Smatrix(self, Smatrix, Z0 = 50, inverse = False):
        """
        Sets self.matrix
        Input:
            Smatrix = 2x2xN np.array
            Z0 = normalization impedance (default = 50 ohm for sonnet files)
        """
        s11 = Smatrix[0,0]
        s21 = Smatrix[1,0]
        s12 = Smatrix[0,1]
        s22 = Smatrix[1,1]
        self.Z0 = Z0
        if inverse:
            s21 = Smatrix[0,1]
            s12 = Smatrix[1,0]
        self.matrix = np.array([[((1+s11)*(1-s22)+s12*s12)/(2*s12),
                         Z0*(((1+s11)*(1+s22)-(s12*s21))/(2*s21))],
                         [(1/Z0)*(((1-s11)*(1-s22)-(s12*s21))/(2*s21)),
                          (((1-s11)*(1+s22)+(s12*s21))/(2*s21))]])

    def setMatrix_line(self, length, Z0, alpha, beta, F = [], eps_eff = -1):
        if F != []:
            self.F = F
            if np.any(beta == -1):
                beta = 2*np.pi*(self.F/spc.c)*np.sqrt(eps_eff) # propagation constant
            alpha = alpha*(self.F/(self.F[-1]+self.F[0])*2)**3 # alpha, set positive value in dB/mm for lossy case
            Z0 = np.full(self.F.shape, Z0)
        self.eps_eff= np.full(self.F.shape, eps_eff)
        self.a = alpha
        self.b = beta
        self.Z0 = Z0
        self.length = length
        self.matrix = np.array([[np.cosh((alpha*length)+(1j*beta*length)),
                         Z0*np.sinh((alpha*length)+(1j*beta*length))],
        [(1/Z0)*np.sinh((alpha*length)+(1j*beta*length)),
         np.cosh((alpha*length)+(1j*beta*length))]])

    def multiply(self, *args):
        """
        Returns new matrix. Out = Self*args[0]*args[1]*...*args[N]
        """
        M = np.zeros_like(self.matrix)
        for i in xrange(len(self.F)):

            M[:,:,i] = self.matrix[:,:,i].dot(args[0].matrix[:,:,i])
            if len(args) > 1:
                for arg in args[1:]:
                    M[:,:,i] = M[:,:,i].dot(arg.matrix[:,:,i])
        out = ABCDmatrix()
        if np.any(self.eps_eff != None):
            out.eps_eff = self.eps_eff
        else:
            out.eps_eff = args[0].eps_eff
        out.F = self.F
        out.Fresolution = self.Fresolution
        out.matrix = M
        return out

