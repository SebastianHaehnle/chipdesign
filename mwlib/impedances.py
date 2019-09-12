# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 17:18:44 2016

This reads the appropriate touchstone file and writes to file:
(line geometry (from parameter))
(kinetic inductances (from parameter))
Z0
epseff

@author: sup-shahnle
"""

import chipdesign
import numpy as np
import scipy.interpolate as spi
from copy import *

def generate_Lk2z(pro, lkname, debug = False):    
    p0 = deepcopy(pro.params)
    p0.pop(lkname)
    Lkz0 = []
    for i, pri in enumerate(pro):
        if all([item in pri.params.iteritems() for item in p0.iteritems()]):
            findex = len(pri.ff)/2-1
            Lkz0.append([pri.params[lkname]*1e-12,pri.fz0[findex][0].real])
    Lkz0 = np.array(sorted(Lkz0, key = lambda x: x[0])).transpose()
    fct = spi.interp1d(Lkz0[0], Lkz0[1])
    if debug:
        print Lkz0
    return fct
    
def generate_Lk2eps(pro, lkname):
    p0 = deepcopy(pro.params)
    p0.pop(lkname)
    Lkeps = []
    for i, pri in enumerate(pro):
        if all([item in pri.params.iteritems() for item in p0.iteritems()]):
            findex = len(pri.ff)/2-1
            Lkeps.append([pri.params[lkname]*1e-12,pri.feps[findex][0].real])
    Lkeps = np.array(sorted(Lkeps, key = lambda x: x[0])).transpose()
    fct = spi.interp1d(Lkeps[0], Lkeps[1])
    return fct

def generate_f2Lk(sc, frange):
    Lkarr = np.zeros(len(frange))
    for i, f in enumerate(frange):
        sc.update(f = f)
        Lkarr[i] = sc.Ls
    fct = spi.interp1d(frange, Lkarr)
    return fct

def generate_f2zeps(f2lk, lk2z0, lk2eps):
    def f2zeps(f):
        return lk2z(f2lk(f)), lk2eps(f2lk(f)) 
    return f2zeps

def param_index(pro, params):
    for i, sweep in enumerate(pro.sweep):
        if all(item in params.iteritems() for item in sweep.params.iteritems()):
            return i

def print_fz0(pro, index = 0):    
    print pro.basename, '\tf = ', pro.ff[index], ';\tZ0 = ', \
    pro.fz0[index][0].real, ';\tEeff = ', pro.feps[index][0].real

def impedances_projects():
    '''
    this returns all Project instances in the impedances folder. Update this according to usage
    '''
    pathmsl = r'C:\Users\sebastian\ownCloud\Chip Design\msloc1\sonnet\msloc1_v1\impedances\msl\msl_Lkdependence\MSL'
    pathhybrid_nc = r'C:\Users\sebastian\ownCloud\Chip Design\msloc1\sonnet\msloc1_v1\impedances\hybrid\output_nc\impedance_hybrid_nc'
    pathhybrid_sc = r'C:\Users\sebastian\ownCloud\Chip Design\msloc1\sonnet\msloc1_v1\impedances\hybrid\output_sc\impedance_hybrid_sc'
    pathwide = r'C:\Users\sebastian\ownCloud\Chip Design\msloc1\sonnet\msloc1_v1\impedances\wide\output_cpw\CPW_NbTiN'
    pathro = r'C:\Users\sebastian\ownCloud\Chip Design\msloc1\sonnet\msloc1_v1\impedances\readout\output_ro\cpw_ro'
    wmsl = 1.4
    hmsl = 1.0

    whyb = 1.4
    shyb = 2.4

    wwide = 6.0
    swide = 16.0
    
    wro = 20
    sro = 10    
    
    Lk = 1.61
    Lkgnd = 0.71
    Lkal = 0.33
    Ral = 0.28

    msl = chipdesign.Project(pathmsl + '.son', {'Width' : wmsl, 'h' : hmsl, 'Lk' : Lk, 'Lk_gnd' : Lkgnd})
    hybridnc = chipdesign.Project(pathhybrid_nc + '.son', {'W_Al' : whyb, 'S_CPW' : shyb, 'Ral' : Ral, 'Lk_gnd' : Lkgnd})
    hybridsc = chipdesign.Project(pathhybrid_sc + '.son', {'W_Al' : whyb, 'S_CPW' : shyb, 'Lk_Al' : Lkal, 'Lk_gnd' : Lkgnd})     
    wide = chipdesign.Project(pathwide + '.son', {'W_line' : wwide, 'S_CPW' : swide, 'Lk_gnd' : Lkgnd})    
    ro = chipdesign.Project(pathro + '.son', {'s' : wro, 'w' : sro, 'Lk_gnd' : Lkgnd})
    allpros = [msl, hybridnc, hybridsc, wide, ro]      
    return allpros
    
  
if __name__ == '__main__':
    outputname = 'impedances.txt'
    
    pathmsl = r'impedances\msl\msl_Lkdependence\MSL'
    pathhybrid_nc = r'impedances\hybrid\output_nc\impedance_hybrid_nc'
    pathhybrid_sc = r'impedances\hybrid\output_sc\impedance_hybrid_sc'
    pathwide = r'impedances\wide\output_cpw\CPW_NbTiN'
    pathro = r'impedances\readout\output_ro\cpw_ro'

#==============================================================================
# line geometries and kinetic inductances
#==============================================================================
    wmsl = 1.4
    hmsl = 1.0

    whyb = 1.4
    shyb = 2.4

    wwide = 6.0
    swide = 16.0
    
    wro = 19
    sro = 8   
    
    Lk = 1.61
    Lkgnd = 0.71
    Lkal = 0.33
    Ral = 0.28

    msl = chipdesign.Project(pathmsl + '.son', {'Width' : wmsl, 'h' : hmsl, 'Lk' : Lk, 'Lk_gnd' : Lkgnd})
    hybridnc = chipdesign.Project(pathhybrid_nc + '.son', {'W_Al' : whyb, 'S_CPW' : shyb, 'Ral' : Ral, 'Lk_gnd' : Lkgnd})
    hybridsc = chipdesign.Project(pathhybrid_sc + '.son', {'W_Al' : whyb, 'S_CPW' : shyb, 'Lk_Al' : Lkal, 'Lk_gnd' : Lkgnd})     
    wide = chipdesign.Project(pathwide + '.son', {'W_line' : wwide, 'S_CPW' : swide, 'Lk_gnd' : Lkgnd})    
    ro = chipdesign.Project(pathro + '.son', {'s' : sro, 'w' : wro, 'Lk_gnd' : Lkgnd})    

#%%
#==============================================================================
# Read sonnet files
#==============================================================================
    
    allpros = np.array([msl, hybridnc, hybridsc, wide, ro]) 
    mask = np.where([1, 1, 0, 0, 1])
    maskedpros = np.take(allpros, mask)[0]
    for pro in maskedpros:
        pro.read_result(True)
    print '\n'
    for pro in maskedpros:
        print_fz0(pro)
    print '\n'

#%%    
#==============================================================================
# Do things with it
#==============================================================================
    msl.read_result(True)
    NbTiNthin = chipdesign.Superconductor(15.2, 190e-8, f = 600e9, t = 120e-9)
    frange = np.arange(100e9, 950e9, 25e9)
    f2lk = generate_f2Lk(NbTiNthin, frange)    
    lknew = f2lk(frange)
    lkname = 'Lk'
    lk2eps = generate_Lk2eps(msl, lkname)
    lk2z = generate_Lk2z(msl, lkname)
    f2zeps = generate_f2zeps(f2lk, lk2z, lk2eps)
    f = 860e9
    lamb2 = chipdesign.f_2_lambda(f, f2zeps(f)[1])*1e6/2.
    print lamb2
            