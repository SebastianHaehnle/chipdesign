# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 14:10:54 2016

@author: sebastian
"""

import numpy as np

def coth(x):
    return np.cosh(x)/np.sinh(x)
    
def csch(x):
    return 1/np.sinh(x)