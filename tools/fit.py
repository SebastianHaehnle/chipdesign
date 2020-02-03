import numpy as np
import scipy as sp
import scipy.optimize as spopt

def lorentzian_log(x, f0, fwhm, I):
    hwhm = fwhm /2.
    return 10*np.log10(I * (hwhm**2 / ((x-f0)**2 + hwhm**2) ))

def lorentzian(x, f0, fwhm, I):
    hwhm = fwhm /2.
    return I * (hwhm**2 / ((x-f0)**2 + hwhm**2) )

def fit_lorentzian(f, s21, pguess, dtype = 'mag', fitmethod = 'log'):
    func = {'mag' : lambda x:20*np.log10(x),
            'db' : lambda x:x if fitmethod == 'log' else 10**(x/10.)}
    fitfunc = {'log' : lorentzian_log,
               'linear' : lorentzian}
    popt, pcov = spopt.curve_fit(fitfunc[fitmethod], f, func[dtype](s21), p0 =pguess)
    return popt, pcov

if __name__ == '__main__':
    f0 = 350e9
    fwhm = 1e9
    f_sample = np.arange(340e9, 360e9, 0.01e9)
    s_sample = lorentzian_log(f_sample, f0, fwhm, 1)
    
    
    import matplotlib.pyplot as plt 
    fig,ax = plt.subplots(1,1)
    ax.plot(f_sample, s_sample)
    
    