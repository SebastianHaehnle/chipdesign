# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 17:42:42 2017

@author: sebastian
"""

import chipdesign
import mwlib as mw
import matplotlib.pyplot as plt
import numpy as np
import os

sweeppath = r'C:\Users\sebastian\ownCloud\Chipdesign\General\Sonnet\impedance_ms\shuttled\\'
sweepfile = os.path.join(sweeppath, 'MSL_40nmShuttled.son')
hmsl = 2
Lk = 4.38
Lkgnd = 0.71
wmsl = 1.4

msl = chipdesign.Project(sweepfile, {'Width' : wmsl, 'h' : hmsl, 'Lk' : Lk, 'Lk_gnd' : Lkgnd})
msl.read_result(True)
print msl.params
#%%
listZ = []
listeeff = []
listW = []
for pro in msl:
    listW.append(pro.params['Width'])
    listZ.append(pro.fz0[0][0].real)
#    listeeff.append(pro.feps[0][0].real)
output = np.array([listW, listZ])
header = '#w_msl [um]   Impedance [Ohm]'
np.savetxt(os.path.join(sweeppath, 'Z_of_w_microstrip_4.38pHsq_2umthickness.dat'), output.transpose(), header = header)

plt.close(35836245)
fig, ax = plt.subplots(1,1, num = 35836245)
plt.plot(listW, listZ)
#plt.plot(listW, listeeff)