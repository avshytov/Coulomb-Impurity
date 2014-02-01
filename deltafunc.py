import math
import numpy as np
from scipy import special

class DeltaCubic:
    def __init__ (self, r_0):
        self.r_0 = r_0
        
    def rho(self, r):
        return 1.0 / 16.0 * self.r_0 / np.sqrt(r**2 + self.r_0**2)**3
    def U(self, r):
        return 1.0 / np.sqrt(r**2 + self.r_0**2)
    
class DeltaGauss:
    def __init__ (self, r_0):
        self.r_0 = r_0
    def rho(self, r):
        return 1.0 / np.pi / 2.0 / self.r_0**2 * np.exp(-r**2 / 2.0 / self.r_0**2)
    def U(self, r):
        C = np.sqrt(2.0 * np.pi) / self.r_0 / 2.0;
        x = r**2 / 4.0 / self.r_0**2
        # ive is bessel function of
        # imaginary argument * decaying exponent
        return special.ive(0, x) * C 

