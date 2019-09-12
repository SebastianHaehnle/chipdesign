# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 16:16:22 2016

cd E:\data\Sebastian\Sonnet\MSLOC1\pysontest

@author: sup-shahnle
"""



import os
from glob import glob
import mwavepy_mod as mw
import paramsweep as ps

def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return float('.'.join([i, (d+'0'*n)[:n]]))
       
class Project(mw.Network):
    def __init__(self, fullpath, variables, sweep = False, **kwargs):
        self.path = fullpath
        self.outputfolder = ''
        self.params = variables
        for k, v in self.params.iteritems():
            self.params[k] =  truncate(v, 14)
        self.sweepbool = sweep
        self.sweep = None
        for k,v in kwargs.items():
            setattr(self, k, v)
        mw.Network.__init__(self)
        
    @property
    def filename(self):
        return os.path.basename(self.path)
    
    @property
    def basename(self):
        return self.filename[:-4]

    @property
    def folder(self):
        return os.path.split(self.path)[0]
        
    def param_index(self, debug = False):
        '''
        only use with sweep initialized
        '''
        for i, sweep in enumerate(self.sweep):
            if debug:
                pass
#                print [item in sweep.params.iteritems() for item in self.params.iteritems()]
            if all(item in sweep.params.iteritems() for item in self.params.iteritems()):
                return i        
        
    def read_result(self, sweep = None, debug = False):
        if sweep == None:
            pass
        else:
            self.sweepbool = sweep
        if self.sweepbool:
            if self.outputfolder == '':
                self.sweep = ps.ParamSweep(self.path[:-4]+'_output_files.txt')
                i = self.param_index(debug)
                if debug:
                    print self.basename, '; default output == index ', i
                if i:                
                    self.read_touchstone(self.sweep[i].path)
            else:
                self.sweep = ps.ParamSweep(self.outputfolder+self.basename+'_output_files.txt')
                i = self.param_index(debug)
                self.read_touchstone(self.sweep[i].path)
        else:
            if self.outputfolder == '':
                tsname = glob(self.path[:-4]+'.s*p')
            else:
                tsname = glob(self.outputfolder + self.basename+'.s*p')                
            if tsname != []:
                self.read_touchstone(tsname[0])
            else:
                print 'No touchstone file found for '+self.filename
    
            ## methods to be like dict(), with inheritance of Mapping
    def __len__(self):
        return len(self.sweep)
    def __getitem__(self, index):
        return self.sweep[index]
    def __iter__(self):
        for i in range(len(self)):
            yield self.sweep[i]
    def __contains__(self, item):
        return item in self.sweep


    def run_project(self, eng):       
        pro = eng.SonnetProject(self.path)
            
        newsim = glob(self.folder + '\\' + self.basename+'.s*p') == []
        for k,v in self.params.iteritems():
            if eng.getVariableValue(pro, k) != v:
                eng.modifyVariableValue(pro, k, v, nargout = 0)
                newsim = True
                print k, eng.getVariableValue(pro, k)
        eng.save(pro, nargout = 0)
        if newsim:
            print 'Running new simulation for ' + self.filename
            eng.simulate(pro, nargout = 0)
            print 'Simulation finished, close the Sonnet analysis window to continue. \n'
        else:
            print 'Correct response data already exists for ' + self.filename + '\n'
            
    def apply_on_all(self, method, *args, **kwargs):
        for obj in self:
             getattr(obj, method)(*args, **kwargs)

if __name__ == '__main__':
    pass
#    pass
#    prolib = {}
#    whyb = 1.4
#    shyb = 2.4
#    Ral = 0.28
#    Lkgnd = 0.71
#        
#    
#    cpwnb = Project(r'E:\data\Sebastian\Sonnet\msloc1_v1\impedances\hybrid\output_nc\impedance_hybrid_nc.son', {'W_Al' : whyb, 'S_CPW' : shyb, 'Ral' : Ral, 'Lk_gnd' : Lkgnd}, sweep = True)
#    cpwnb.read_result(True)
#    print cpwnb[0]
#    
#    cpwalsc = Project(r'KID\CPW_Al_sc.son',{})
##    cpwal.read_result()
#    
#    mslthz = Project(r'MSL\MSL.son', {'Width' : 1.5, 'h' : 0.5})
#    eng = matlab.engine.start_matlab()
#    eng.cd(r'E:\\data\\Sebastian\\Sonnet\\MSLOC1')
###    eng.addpath(opdir)
#    eng.addpath('E:\data\Sebastian\pySon\SonnetLab\Scripts')
#    pro = eng.SonnetProject(mslthz.path)
#    print eng.returnSelectedFrequencySweep(pro)
#    eng.quit()
#    def run_project(project, eng):       
#        pro = eng.SonnetProject(project.path)
#            
##        newsim = glob(project.folder + '\\' + project.basename+'.s*p') == []
##        for k,v in project.vars.items():
##            if eng.getVariableValue(pro, k) != v:
##                eng.modifyVariableValue(pro, k, v, nargout = 0)
##                newsim = True
##                print k, eng.getVariableValue(pro, k)
#    #    eng.saveAs(pro, project.path, nargout = 0)
#        
#        eng.save(pro, nargout = 0)
#        newsim = True
#        if newsim:
#            print 'Running new simulation for ' + project.name
#        #    eng.addFileOutput(pro, )
#            eng.simulate(pro, nargout = 0)
#            print 'Simulation finished, close the Sonnet analysis window to continue. \n'
#        else:
#            print 'Correct response data already exists for ' + project.name + '\n'
#    run_project(mslson, eng)
#
#    #pro = eng.SonnetProject(mslson.path)
#    #eng.save(pro, nargout = 0)
#    
#    #pro = eng.SonnetProject('MSL\\MSL.son')
#    eng.quit()
#    print 'Exiting matlab engine'
        
    
