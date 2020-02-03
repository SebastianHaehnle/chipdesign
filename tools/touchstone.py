# -*- coding: utf-8 -*-
"""
Created on Wed May 02 13:15:50 2018

@author: sebastian
"""

import chipdesign.mwavepy_mod as mwmod
import numpy as np

class Touchstone(mwmod.Network):
    def __init__(self, filepath, params = {}):
        mwmod.Network.__init__(self)
        self.read_touchstone(filepath)
        self.params = params
        if self.params == {}:
            with open(filepath, 'r') as ff:
                active = False
                lines = ff.readlines()
                for line in lines:
                    if active == True:
                        if line[3:-1] == 'END PARAMS':
                            active = False
                        else:
                            eq_id = line.index('=')
                            self.params[line[3:eq_id-1]] = float(line[eq_id+2:-1])
                    elif line[3:-1] == 'PARAMS':
                        active = True

    @property
    def s12(self):
        return self.s[:,0,1]

    @property
    def s21(self):
        return self.s[:,1,0]

    @property
    def s11(self):
        return self.s[:,0,0]

    @property
    def s22(self):
        return self.s[:,1,1]

    @property
    def s31(self):
        return self.s[:,2,0]

    @property
    def s32(self):
        return self.s[:,2,1]

    @property
    def s33(self):
        return self.s[:,2,2]

class CST_Plotdata(object):
    header_match = {'S1,1/abs,dB' : 's11',
                    'S2,1/abs,dB' : 's21',
                    'S1,2/abs,dB' : 's12',
                    'S2,2/abs,dB' : 's22',
                    'Frequency' : 'f'}
    setter_functions = {'S1,1/abs,dB' : lambda x :10**(x/20.),
                        'S2,1/abs,dB' : lambda x :10**(x/20.),
                        'S1,2/abs,dB' : lambda x :10**(x/20.),
                        'S2,2/abs,dB' : lambda x :10**(x/20.),
                        'Frequency' : lambda x : x*1e9}
    string_datastart = '----------------------------------------------------------------------\n'
    def __init__(self, filepath, *args):
        self.filepaths = [filepath] + [arg for arg in args]
#        print self.filepaths
        for filepath in self.filepaths:
            self.add_file(filepath)

    def add_file(self, filepath):
        if not filepath in self.filepaths:
            self.filepaths.append(filepath)
        with open(filepath, 'r') as ff:
            self.lines = ff.readlines()
            self.header_ids = [ii-1 for ii in xrange(len(self.lines)) if self.lines[ii] == self.string_datastart]
            self._data_dict = {}
            for jj, _id in enumerate(self.header_ids):
                # parse header to remove crap
                header_list = self.lines[_id].split()
                separator_ids = [ii for ii in xrange(len(header_list)) if header_list[ii] == '/']
                for sid in separator_ids[::-1]:
                    del header_list[sid:sid+2]
                assert len(header_list) == 2
                # create data slices
                if _id == self.header_ids[-1]:
                    slice_stop = -1
                else:
                    slice_stop = self.header_ids[jj+1]-1
                data_list = self.lines[_id+2: slice_stop]
                # split strings, convert to np.array, convert to float array, transpose array
                data_list = np.array(map(str.split, data_list)).astype(float).transpose()
                self._data_dict[header_list[0]] = data_list[0]
                self._data_dict[header_list[1]] = data_list[1]

        for k, v in self._data_dict.iteritems():
            if k in self.header_match:
                setattr(self, self.header_match[k], self.setter_functions.get(k, lambda x: x)(v))
            else:
                print "WARNING: header does not have a proper header_match"
                setattr(self, k, v)

    @property
    def s(self):
        try:
            return np.array([[self.s11, self.s12],[self.s21, self.s22]])
        except:
            sys.exit('S-parameter not fully initialized in CST touchstone object')

    @property
    def s_db(self):
        try:
            return 20*np.log10(self.s)
        except:
            sys.exit('S-parameter not fully initialized in CST touchstone object')





if __name__ == '__main__':
    filepath = r'C:/Users/sebastian/ownCloud/Chips/FP design/M400x/wideband_chip/cst_shorted/s_parameters.txt'
    foo = CST_Plotdata(filepath)