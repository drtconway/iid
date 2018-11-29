import os
import sys
import yaml
from iid.beta import beta

def sameAbs(a, b, eps):
    return abs(a - b) < eps

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def addTests():
    thisModule = sys.modules[__name__]
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'dist', 'beta.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        a = itm[0]
        b = itm[1]
        x = itm[2]
        m = itm[3]
        v = itm[4]
        rP = itm[5]
        rL = itm[6]
        testName = 'test_beta_%d_%d_%f' % (a, b, x)
        def thisTest():
            dst = beta(a, b)
            assert same(m, dst.mean(), 1e-14)
            assert same(v, dst.var(), 1e-14)
            rP0 = dst.pdf(x)
            assert same(rP, rP0, 1e-12)
            rL0 = dst.cdf(x)
            assert same(rL, rL0, 1e-12)
        setattr(thisModule, testName, thisTest)

addTests()
