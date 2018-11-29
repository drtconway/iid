import os
import sys
import yaml
from iid.binom import binom

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def addTests():
    thisModule = sys.modules[__name__]
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'dist', 'binom.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        n = itm[0]
        p = itm[1]
        k = itm[2]
        m = itm[3]
        v = itm[4]
        rP = itm[5]
        rL = itm[6]
        rU = itm[7]

        testName = 'test_binom_%d_%f_%d' % (n, p, k)
        def thisTest():
            dst = binom(p, n)
            assert same(m, dst.mean(), 1e-14)
            assert same(v, dst.var(), 1e-14)
            rP0 = dst.pmf(k)
            assert same(rP, rP0, 1e-12)
            rL0 = dst.cdf(k)
            assert same(rL, rL0, 1e-12)
        setattr(thisModule, testName, thisTest)

addTests()
