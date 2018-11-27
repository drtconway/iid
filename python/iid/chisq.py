import math
import iid.basic as basic
from iid.dist import dist
import iid.special as special

class chisq(dist):
    def __init__(self, k):
        assert(k >= 1)
        self.k = k

    def mean(self):
        return self.k

    def median(self):
        return self.k*math.pow(1 - 2.0/(9*k), 3)

    def var(self):
        return 2*self.k

    def pdf(self, x, **args):
        assert(self.k != 1 or x > 0)
        assert(self.k <= 1 or x >= 0)
        k = self.k
        lx = math.log(x)
        l2 = math.log(2)
        lr = (0.5*k - 1) * lx  - 0.5*x - 0.5*k*l2 - basic.logGamma(0.5*k)
        if 'log' in args and args['log']:
            return lr
        else:
            return math.exp(lr)

    def cdf(self, x, **args):
        assert(self.k != 1 or x > 0)
        assert(self.k <= 1 or x >= 0)
        k = self.k
        if 'upper' in args and args['upper']:
            lr = special.logGammaQ(0.5*k, 0.5*x)
        else:
            lr = special.logGammaP(0.5*k, 0.5*x)

        if 'log' in args and args['log']:
            return lr
        else:
            return math.exp(lr)

