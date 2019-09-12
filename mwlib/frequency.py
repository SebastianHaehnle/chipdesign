# -*- coding: utf-8 -*-
"""
Created on Wed Nov 02 13:52:12 2016

@author: sebastian
"""

import numpy as np
npa = np.array

class Frequency():
    # e.g. : f [unit:Hz] *unitScale[GHz] = f [unit:GHz]
    unitScale = {'Hz' : 1e0,
                 'kHz' : 1e-3,
                 'MHz' : 1e-6,
                 'GHz' : 1e-9,
                 'ghz' : 1e-9,
                 'THz' : 1e-12}
    def __init__(self, fStart, fStop = None, fStep = None, unit = 'GHz'):
        self.unit = unit
        if fStop and fStep:
            self._f = np.arange(fStart, fStop+fStep, fStep)
        elif isinstance(fStart, Frequency):
            self._f = fStart._f
            self.unit = fStart.unit
        else:
            self._f = npa(fStart)


    @property
    def fmin(self):
        return self._f[0]

    @property
    def fmax(self):
        return self._f[-1]

    @property
    def df(self):
        return self._f[1] - self._f[0]

    @property
    def N(self):
        return len(self)

    @property
    def fscaled(self):
        return self._f*Frequency.unitScale[self.unit]

    @property
    def f(self):
        return self._f

    def __len__(self):
        return len(self._f)
    def __getitem__(self, index):
        return self._f[index]
    def __iter__(self):
        for i in range(len(self)):
            yield self._f[i]
    def __contains__(self, item):
        return item in self._f
    def __eq__(self, other):
        if not isinstance(other, Frequency):
            return False
        elif (self.f == other.f).all():
            return True
        else:
            return False

    def __str__(self):
        return str(self._f)

if __name__ == '__main__':
    f = Frequency(300e9, 900e9, 100e9)
    f2 = Frequency(f.f)
    print f
    print f2
    print f2.f
    print f == f2
