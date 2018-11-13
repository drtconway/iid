import math
import iid.basic as basic
from iid.dist import dist

class hyper(dist):
    def __init__(self, N, K, n):
        self.N = N
        self.K = K
        self.n = n

    def mean(self):
        N = float(self.N)
        K = float(self.K)
        n = float(self.n)
        return n*K/N

    def median(self):
        return None

    def var(self):
        N = float(self.N)
        K = float(self.K)
        n = float(self.n)
        return  n*(K/N)*((N-K)/N)*((N-n)/(N - 1))

    def pmf(self, k, **args):
        assert(k >= 0)
        assert(k <= self.n)
        assert(k <= self.K)
        N = self.N
        K = self.K
        n = self.n
        if 'log' in args and args['log']:
            return basic.logChoose(K, k) + basic.logChoose(N - n, N - k - 1) - basic.logChoose(N, n)
        else:
            return float(basic.choose(K, k)) * float(basic.choose(N - n, N - k - 1)) /  float(basic.choose(N, n))

    def cdf(self, k, **args):
        N = self.N
        K = self.K
        n = self.n

        lr0 = basic.logChoose(n, k+1) + basic.logChoose(N-n, K-k-1) - basic.logChoose(N, K)
        lr1 = basic.logHyper([1, k+1-K, k+1-n], [k+2, N+k+2-K-n], 1)
        lr = lr0 + lr1

        if not ('upper' in args and args['upper']):
            lr = basic.log1mexp(ls)

        if 'log' in args and args['log']:
            return lr
        else:
            return math.exp(lr)

