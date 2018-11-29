import math
import iid.basic as basic
from iid.dist import dist

class hyper(dist):
    def __init__(self, N, K, n):
        self.N = N
        self.K = K
        self.n = n
        self.kMin = max(0, n + K - N)
        self.kMax = min(n, K)

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
        N = self.N
        K = self.K
        n = self.n

        assert self.kMin <= k
        assert k <= self.kMax

        if 'log' in args and args['log']:
            return basic.logChoose(K, k) + basic.logChoose(N - K, n - k) - basic.logChoose(N, n)
        else:
            return float(basic.choose(K, k)) * float(basic.choose(N - K, n - k)) /  float(basic.choose(N, n))

    def cdf(self, k, **args):
        N = self.N
        K = self.K
        n = self.n

        assert self.kMin <= k
        assert k <= self.kMax

        # Find the mode, to decide which direction to do the summation
        #
        w = (n + 1.0)*(K + 1.0)/(N + 2.0)

        if k == self.kMax:
            lls = 0.0
            lus = None
        elif k < w:
            s = 0
            for j in range(k, self.kMin-1, -1):
                t = self.pmf(j)
                s += t
            lls = math.log(s)
            lus = None
        else:
            s = 0
            for j in range(k+1, self.kMax+1):
                t = self.pmf(j)
                s += t
            lls = None
            lus = math.log(s)

        if 'upper' in args and args['upper']:
            if lus is None:
                lr = basic.log1mexp(lls)
            else:
                lr = lus
        else:
            if lls is None:
                lr = basic.log1mexp(lus)
            else:
                lr = lls

        if 'log' in args and args['log']:
            return lr
        return math.exp(lr)


