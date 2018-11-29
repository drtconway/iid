import os
import sys
import yaml
from iid.chisq import chisq

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def addTests():
    thisModule = sys.modules[__name__]
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'dist', 'chisq.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        n = itm[0]
        x = itm[1]
        m = itm[2]
        v = itm[3]
        rP = itm[4]
        rL = itm[5]
        rU = itm[6]
        testName = 'test_chisq_%f_%f' % (n, x)
        def thisTest():
            dst = chisq(n)
            assert same(m, dst.mean(), 1e-14)
            assert same(v, dst.var(), 1e-14)
            rP0 = dst.pdf(x)
            if not same(rP, rP0, 1e-13):
                print(itm, abs((rP-rP0)/rP))
            assert same(rP, rP0, 1e-13)
            rL0 = dst.cdf(x)
            assert same(rL, rL0, 1e-13)
        setattr(thisModule, testName, thisTest)

addTests()
