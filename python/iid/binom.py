import math
import iid.basic as basic
from iid.dist import dist
import iid.special as special

class binom(dist):
    def __init__(self, p, n):
        self.p = p
        self.q = 1 - p
        self.lp = math.log(p)
        self.lq = basic.log1p(-p)
        self.n = n

    def mean(self):
        return self.p * self.n

    def median(self):
        return int(math.floor(self.mean()))

    def var(self):
        return self.n * self.p * self.q

    def pmf(self, k, **args):
        assert(k >= 0)
        assert(k <= self.n)
        if 'log' in args and args['log']:
            lr = basic.logChoose(self.n, k)
            return lr + k*self.lp + (self.n - k)*self.lq
        else:
            v = k*self.lp + (self.n - k)*self.lq
            if v > -700:
                r = basic.choose(self.n, k)
                return r * (self.p**k) * (self.q**(self.n-k))
            else:
                lr = basic.logChoose(self.n, k)
                lr += k*self.lp + (self.n - k)*self.lq
                return math.exp(lr)

    def cdf(self, k, **args):
        ls = self.n*self.lq
        i = 1
        while i <= k:
            lt = basic.logChoose(self.n, i) + i*self.lp + (self.n - i)*self.lq
            ls = basic.logAdd(ls, lt)
            i += 1

        if 'upper' in args and args['upper']:
            ls = basic.log1mexp(ls)

        if 'log' in args and args['log']:
            return ls
        else:
            return math.exp(ls)

    def quant(self, q):
        return self.quantK(q, 0, self.n)

    def rnd(self):
        if self.n < 25:
            k = 0
            for i in xrange(self.n):
                if self.rndX() < self.p:
                    k += 1
            return k
        
        u = self.rndX()
        return self.quant(u)
