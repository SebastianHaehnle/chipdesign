import numpy as np
from .touchstone import Touchstone
from collections import Mapping
import os
from itertools import takewhile
import re


class SonnetSweep(Mapping):
    def __init__(self, path, **kwargs):
        self.path_param = path
        self.folder = os.path.dirname(self.path_param)
        self.sweep_params = np.atleast_1d(kwargs.pop('sweep_params', []))
        self.paths_ts, self.params = self.parse_summary(self.path_param)
        self.paths_ts = [os.path.join(self.folder, path) for path in self.paths_ts]
        self.ts = []
        for ii, path in enumerate(self.paths_ts):
            self.ts.append(Touchstone(self.paths_ts[ii], self.params[ii]))

    @staticmethod
    def parse_summary(path):
        with open(path, 'r') as ff:
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

    def sort_by_param(self, parameter, reverse = False):
        self.ts.sort(key = lambda cc: cc.params[parameter], reverse = reverse)
        self.params.sort(key = lambda cc: cc[parameter], reverse = reverse)

    def sweep_values(self, parameter):
        return [param[parameter] for param in self.params]

    def get_ts(self, parameter, value):
        return self.ts[self.sweep_values(parameter).index(value)]


    def __len__(self):
        return len(self.ts)
    def __getitem__(self, index):
        return self.ts[index]
    def __iter__(self):
        for i in range(len(self)):
            yield self.ts[i]
    def __contains__(self, item):
        return item in self.ts
