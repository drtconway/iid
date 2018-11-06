import random

class dist(object):
    '''base class for distributions with some helper functions'''

    def quantK(self, q):
        '''compute quantiles of discrete distribution by bisection search on the CDF'''
        lo = 0
        hi = self.n

        if self.cdf(lo) >= q:
            return lo
        if self.cdf(hi) <= q:
            return hi

        idx = {}
        while lo < hi:
            mid = (lo + hi) // 2
            idx[mid] = self.cdf(mid)
            if idx[mid] >= q:
                hi = mid - 1
            else:
                lo = mid + 1

        if lo in idx:
            q0 = idx[lo]
        else:
            q0 = self.cdf(lo)

        while q0 < q:
            lo += 1
            q0 += self.pmf(lo)
        assert self.cdf(lo) >= q
        return lo

    def rndK(self, lo, hi):
        return random.randint(lo, hi)

    def rndX(self):
        return random.random()
