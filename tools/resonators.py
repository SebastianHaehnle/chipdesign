import numpy as np
import scipy.constants as spc


class Resonator_Series(object):
    def __init__(self):
        pass

    @property
    def Qi(self):
        return self._Qi

    @Qi.setter
    def Qi(self, Qi):
        self._Qi = Qi

    @classmethod
    def Qarr_from_peak(cls, Q_loaded, s21_max):
        Qc = Q_loaded/s21_max
        Qi = Q_loaded*Qc/(np.abs(Qc - Q_loaded))
        Ql = Q_loaded
        return [Ql, Qi, Qc]
