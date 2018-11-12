import math
import iid.basic as basic
from iid.dist import dist
import iid.special as special

class norm(dist):
    def __init__(self, lam)
        self.lam = lam
        self.llam = math.log(lam)

    def mean(self):
        return self.lam

    def median(self):
        return int(math.floor(self.lam + 1.0/3.0 - 0.02/self.lam))

    def var(self):
        return self.lam

    def pmf(self, k, **args):
        lr = k*self.llam - self.lam - basic.logFac(k)
        if 'log' in args and args['log']:
            return lr
        else:
            return math.exp(lr)

    def cdf(self, x, **args):
        log = ('log' in args and args['log'])
        upper = ('upper' in args and args['upper'])

        if log:
            if upper:
                return special.logGammaP(k+1, self.lam)
            else:
                return special.logGammaQ(k+1, self.lam)
        else:
            if upper:
                return special.gammaP(k+1, self.lam)
            else:
                return special.gammaQ(k+1, self.lam)

