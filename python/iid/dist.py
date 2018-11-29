import random

class dist(object):
    '''base class for distributions with some helper functions'''

    def quantK(self, q, lo, hi):
        '''compute quantiles of discrete distribution by bisection search on the CDF'''

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
        return lo

    def quantX(self, q, lo, hi):
        if self.cdf(lo) >= q:
            return lo
        if self.cdf(hi) <= q:
            return hi

        eps = 1e-9
        while (hi - lo)/hi > eps:
            mid = (lo + hi) / 2.0
            q0 = self.cdf(mid)
            if q0 < q:
                lo = mid
            else:
                hi = mid
        return lo

    def rndK(self, lo, hi):
        return random.randint(lo, hi)

    def rndX(self):
        return random.random()
