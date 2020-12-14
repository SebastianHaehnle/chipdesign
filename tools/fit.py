import numpy as np
import scipy as sp
import scipy.optimize as spopt
import lmfit

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

def fit_lorentzian_curve(x,y, error = False):
    pinit = lorentzian_initial_params(x,y)
    p0 = [pinit['hwhm'], pinit['peak_center'], pinit['intensity'], pinit['offset']]
    popt, pocv = spopt.curve_fit(lorentzian_curvefit, x, y, p0)
    if error:
        return popt, pocv
    return popt

def lorentzian_curvefit(x, p1, p2, p3, p4):
    numerator = p1**2
    denominator = (x-p2)**2 + p1**2
    y = p3*(numerator/denominator) + p4
    return y

def lorentzian_initial_params(x, y):
    params = lmfit.Parameters()
    idx1 = np.where(y >= 0.5*np.amax(y))
    fwhm_hat = (x[1] - x[0]) * len(idx1[0])
    hwhm_hat = 0.5*fwhm_hat
    #params.add('hwhm', value=0.001, min=0)
    params.add('hwhm', value=hwhm_hat, min=0)
    params.add('peak_center', value=x[np.argmax(y)])
    params.add('intensity', value=np.amax(y), min=0)
    #params.add('offset', value=x[np.argmin(y)], min=0)
    params.add('offset', value=y[np.argmin(y)])
    return params


if __name__ == '__main__':
    f0 = 350e9
    fwhm = 1e9
    I = 0.83
    f_sample = np.arange(340e9, 360e9, 0.01e9)
    s_sample = lorentzian_log(f_sample, f0, fwhm, I)
    popt, pocv = fit_lorentzian(f_sample, s_sample, [f0, fwhm, I], dtype = 'db', fitmethod = 'linear')
    s_fit = lorentzian_log(f_sample, *popt)
    print popt
    import matplotlib.pyplot as plt 
    fig,ax = plt.subplots(1,1)
    ax.plot(f_sample, s_sample)
    ax.plot(f_sample, s_fit)
    
    