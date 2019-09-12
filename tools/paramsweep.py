# -*- coding: utf-8 -*-
"""
Created on Fri Jul 08 15:33:27 2016

@author: sebastian
"""

from collections import Mapping
from itertools import *
import re
import os
import numpy as np
import mwavepy_mod as mw

class ParamSweep(Mapping):
    def __init__(self, filename, **kws):
        self.path = filename
        self.folder = os.path.dirname(self.path)
        self.ts_namelist, self.params = self.parse_summary()
        self._ts = [Touchstone(self, i, name, self.params[i]) for (i, name) in enumerate(self.ts_namelist)]

    def parse_summary(self, ):
        with open(self.path, 'r') as ff:
            fstr = ff.read()
            tsnames = []
            params = []
#            print 'reading ', self.path
            for i in (m.end() for m in re.finditer('OUTPUT_FILE ', fstr)):
                tsnames.append(''.join(takewhile(lambda x : x!='\n', fstr[i:])))
                temp = fstr[i:].split()
                temp = temp[2:temp.index('END')]
                params.append({temp[3*i] : float(temp[3*i+2]) for i in range(len(temp)/3)})
        return tsnames, params

    def printf(self, index = None):
        if index == None:
            filepath = self.path
        else:
            filepath = self._ts[index].path
        with open(filepath, 'r') as ff:
            print ff.read()

    ## methods to be like dict(), with inheritance of Mapping
    def __len__(self):
        return len(self._ts)
    def __getitem__(self, index):
        return self._ts[index]
    def __iter__(self):
        for i in range(len(self)):
            yield self._ts[i]
    def __contains__(self, item):
        return item in self._ts


class Touchstone(mw.Network):
    def __init__(self, parent, index, filename, params, **kws):
        self._parent = parent
        self._index = index
        self.path = os.path.join(parent.folder, filename)
        self.params = params
        mw.Network.__init__(self, self.path)
#        self._fmode = {'GHZ' : 1e9}


    def printf(self):
        with open(self.path, 'r') as ff:
            print ff.read()


#==============================================================================
# Incomplete touchstone parsing implementation. Using mwavelib_mod instead
#==============================================================================
#    def parse_data(self):
#        self._data = self.rawdata
#        self._mode = self.mode
#        self.f = data[0]*self._fmode[self.mode[1]]
#        print data
#
#
#    def _read_data(self):
#        with open(self.path, 'r') as ff:
#            lines = [line.strip('\n') for line in ff.readlines() if '!' not in line ]
#            lines = lines[1:]
##            data = np.array([np.fromstring(line) for line in lines])
#            data = np.array([np.fromstring(line, sep=' ') for line in lines])
#        return data
#
#    def _read_mode(self):
#        with open(self.path, 'r') as ff:
#            fstr = ff.read()
#            modestr = ''.join(takewhile(lambda x : x!='\n', fstr[fstr.index('#')+1:]))
#            modestr = modestr.split()
#        modestr = re.findall(r'\d+', self.path[self.path.index('.')+1:]) + modestr
#        return modestr
#
#    def __getattribute__(self, name):
#        if name == 'mode':
#            return self._read_mode()
#        elif name == 'rawdata':
#            return self._read_data()
#        else:
#            return object.__getattribute__(self, name)


if __name__ == '__main__':
    pass
#    path = 'parser_example'
#    txtname = 'CPW_output_files.txt'
#
#    filepath = path + '/' +txtname
#    tp = ParamSweep(filepath)


#    print tp[0].mode
#    data = tp[0]._read_data()
#    mode = tp[0]._read_mode()
#    print mode
#    print data[-1]

