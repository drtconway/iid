import math
import iid.basic as basic
from iid.dist import dist
import iid.special as special

class geom(dist):
    def __init__(self, p):
        self.p = p
        self.q = 1 - p
        self.lp = math.log(p)
        self.lq = basic.log1p(-p)

    def mean(self):
        return self.q / self.p

    def median(self):
        return int(math.ceil(-1/self.lq)) - 1

    def var(self):
        return self.q/(self.p*self.p)

    def pmf(self, k, **args):
        assert(k >= 0)
        if 'log' in args and args['log']:
            return k*self.lq + self.lp
        else:
            return math.pow(self.q, k)*self.p

    def cdf(self, k, **args):
        lr = (k+1)*self.lq
        if not ('upper' in args and args['upper']):
            lr = basic.log1mexp(lr)

        if 'log' in args and args['log']:
            return lr
        else:
            return math.exp(lr)

