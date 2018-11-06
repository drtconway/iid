import math
import pystat.basic as basic
from pystat.dist import dist
import pystat.special as special

class norm(dist):
    def __init__(self, mu = 0.0, sig = 1.0):
        self.mu = mu
        self.sig = sig
        self.sr2 = sig * math.sqrt(2)
        self.C = 1.0/math.sqrt(2*basic.pi*sig*sig)
        self.lC = math.log(self.C)

    def mean(self):
        return self.mu

    def median(self):
        return self.mu

    def var(self):
        return self.sig**2

    def pdf(self, x, **args):
        z = (x - self.mu)**2 / (2*self.sig*self.sig)
        if 'log' in args and args['log']:
            return self.lC - z
        else:
            return self.C * math.exp(-z)

    def cdf(self, x, **args):
        log = ('log' in args and args['log'])
        upper = ('upper' in args and args['upper'])

        z = (x - self.mu)/self.sr2
        if log:
            if upper:
                return 0.5 * special.logErfc(z)
            else:
                return 0.5 * special.logErfc(-z)
        else:
            if upper:
                return 0.5 * special.erfc(z)
            else:
                return 0.5 * special.erfc(-z)

