import os
import sys
import yaml
from iid.geom import geom

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def addTests():
    thisModule = sys.modules[__name__]
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'dist', 'geom.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        p = itm[0]
        k = itm[1]
        m = itm[2]
        v = itm[3]
        rP = itm[4]
        rL = itm[5]
        rU = itm[6]
        dst = geom(p)

        testName = 'test_geom_%f_%d' % (p, k)
        def thisTest():
            assert same(m, dst.mean(), 1e-14)
            assert same(v, dst.var(), 1e-14)
            rP0 = dst.pmf(k)
            assert same(rP, rP0, 1e-14)
            rL0 = dst.cdf(k)
            assert same(rL, rL0, 1e-14)
        setattr(thisModule, testName, thisTest)

addTests()
