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
        N = self.N
        K = self.K
        n = self.n

        assert k >= max(0, n + K - N)
        assert k <= min(n, K)

        if 'log' in args and args['log']:
            return basic.logChoose(K, k) + basic.logChoose(N - K, n - k) - basic.logChoose(N, n)
        else:
            return float(basic.choose(K, k)) * float(basic.choose(N - K, n - k)) /  float(basic.choose(N, n))

    def cdf(self, k, **args):
        N = self.N
        K = self.K
        n = self.n

        assert k >= max(0, n + K - N)
        assert k <= min(n, K)

        if k == min(n, K):
            return 1.0

        w = (n + 1.0)*(K + 1.0)/(N + 2.0)

        if k < w:
            s = 0
            for j in range(max(0, n + K - N), k+1):
                t = self.pmf(j)
                s += t
            lls = math.log(s)
            lus = None
        else:
            s = 0
            for j in range(k+1, min(n,K)+1):
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

        #lr0 = basic.logChoose(n, k+1) + basic.logChoose(N-n, K-k-1) - basic.logChoose(N, K)
        #lr1 = basic.logHyper([1, k+1-K, k+1-n], [k+2, N+k+2-K-n], 1)
        #lr = lr0 + lr1
        #if not ('upper' in args and args['upper']):
        #    lr = basic.log1mexp(lr)

        #if 'log' in args and args['log']:
        #    return lr
        #else:
        #    return math.exp(lr)

