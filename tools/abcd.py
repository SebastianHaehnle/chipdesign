# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 11:41:20 2019

@author: sebastian
"""

import numpy as np
import scipy.constants as spc
import numexpr as ne
from multiprocessing import Pool
pool = Pool()

class ABCDmatrix(object):
    def __init__(self, Fresolution = -1):
        self.matrix = np.zeros((2,2))
        self.F = np.array([0])
        self.Fresolution = Fresolution
        self.eps_eff = None
        self.length = 0


    def getS21(self, Z0char, Z02 = None, otype = 'db'):
        """
        If Z02 is given, Z0char is set as port 1 impedance and Z02 as port 2 impedance. Otherwise Z0char is the impedance for both ports.
        
        Possible options for otype: 'db', 'linear', 'imag', 'real', 'complex'
        """
        Z0char = np.atleast_1d(Z0char)
        if len(Z0char) != len(self.F):
            Z0char = np.full(len(self.F), np.mean(Z0char))
        fct = {'db' : lambda x : 20*np.log10(np.abs(x)),
               'linear' : lambda x : np.abs(x),
               'imag' : lambda x : np.imag(x),
               'real' : lambda x : np.real(x),
               'complex' : lambda x : x,
               'power' : lambda x: np.abs(x)**2
               }
        Z02 = np.atleast_1d(Z02)
        if (Z02 != None).any():
            return fct[otype]((2*np.sqrt(Z0char*Z02))/(Z02*self.matrix[0,0,:] + (self.matrix[0,1,:]) + self.matrix[1,0,:]*Z0char*Z02 + self.matrix[1,1,:]*Z0char))
        else:
            return fct[otype](2/(self.matrix[0,0,:] + (self.matrix[0,1,:]/Z0char) + self.matrix[1,0,:]*Z0char + self.matrix[1,1,:]))


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

    def setMatrix_C(self, C, F):
        self.F = F
        A = np.ones_like(self.F)
        B = 1/(1j*2*np.pi*self.F*C)
        C = np.zeros_like(self.F)
        D = np.ones_like(self.F)
        self.matrix = np.array([[A, B],[C, D]])


    def setMatrix_ts(self, ts,  ctype = 'sparam', symmetric = True , inverse = False, length = 0, alpha = 0, eps_force = 0, F = [], Z0 = 0):
        """
        alpha in [dB/mm]
        """
        if len(F) != 0:
            self.F = np.array(F)
            s = self.interpolate_s(self.F, ts.f, ts.s, symmetric)
        else:
            s = ts.s
            self.F = np.arange(ts.f[0], ts.f[-1] + self.Fresolution, self.Fresolution)
            s = self.interpolate_s(self.F, ts.f, ts.s, symmetric)

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

    def setMatrix_line(self, length, Z0, alpha, beta, F = [], eps_eff = -1, n = 3):
        """
        alpha :: alpha given at the center frequencyNp/m
        n :: frequency dependence of alpha, default is 3 for radiation loss. 
        """
        if F != []:
            self.F = F
            if np.any(beta == -1):
                beta = 2*np.pi*(self.F/spc.c)*np.sqrt(eps_eff) # propagation constant
            alpha = alpha*(self.F/(self.F[-1]+self.F[0])*2)**n # alpha, set positive value in dB/mm for lossy case
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

    def multiply(self, input_ABCD, parallelize = True):
        """
        Returns new matrix. Out = Self*args[0]*args[1]*...*args[N]
        Do not use parallelize = False, only for legacy purposes. For Fres <0.001e9 parallelize = True speeds up about 10 times
        """
        M = np.zeros_like(self.matrix)
        if parallelize:
            M = np.einsum('abi,bci->aci', self.matrix, input_ABCD.matrix)
        else:
            for i in xrange(len(self.F)):
                M[:,:,i] = self.matrix[:,:,i].dot(input_ABCD.matrix[:,:,i])
        out = ABCDmatrix()
        if np.any(self.eps_eff != None):
            out.eps_eff = self.eps_eff
        else:
            out.eps_eff = input_ABCD.eps_eff
        out.F = self.F
        out.Fresolution = self.Fresolution
        out.matrix = M
        return out


    def interpolate_s(self, fo, fi, s, symmetric = False):
        """
        s input form: [freq, portA, portB], i.e. s21 = s[:,1, 0], s32 = s[:,2,1]
    
        returns interpolated complex s array
        """
        n_ports = s.shape[-1]
    
        output = np.zeros((len(fo), n_ports, n_ports), dtype = complex)
    
        if symmetric:
            mag = np.abs(s[:,0,1])
            angle = np.angle(s[:,0,1])
            new_mag = np.interp(fo, fi, mag)
            new_angle = np.interp(fo, fi, angle)
            output[:,0,1] = new_mag*np.exp(1j*new_angle)
            output[:,1,0] = output[:,0,1]
            
            mag = np.abs(s[:,0,0])
            angle = np.angle(s[:,0,0])
            new_mag = np.interp(fo, fi, mag)
            new_angle = np.interp(fo, fi, angle)
            output[:,0,0] = new_mag*np.exp(1j*new_angle)
            output[:,1,1] = output[:,0,0]
            
        else:
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