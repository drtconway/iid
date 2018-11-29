import os
import sys
import yaml
from iid.special import gammaP, gammaQ, gammaPQ, logGammaP, logGammaQ, logGammaPQ

def sameAbs(a, b, eps):
    return abs(a - b) < eps

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def addTests():
    thisModule = sys.modules[__name__]
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'special', 'igamma.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        a = itm[0]
        x = itm[1]
        rL = itm[2]
        lrL = itm[3]
        rU = itm[4]
        lrU = itm[5]
        testName = 'test_igamma_%f_%f' % (a, x)
        def thisTest():
            rL0 = gammaP(a, x)
            assert same(rL, rL0, 1e-12)
            rU0 = gammaQ(a, x)
            assert same(rU, rU0, 1e-12)
            (rL0, rU0) = gammaPQ(a, x)
            assert same(rL, rL0, 1e-12)
            assert same(rU, rU0, 1e-12)

            lrL0 = logGammaP(a, x)
            assert same(lrL, lrL0, 1e-12)
            lrU0 = logGammaQ(a, x)
            assert same(lrU, lrU0, 1e-12)
            (lrL0, lrU0) = logGammaPQ(a, x)
            assert same(lrL, lrL0, 1e-12)
            assert same(lrU, lrU0, 1e-12)
        setattr(thisModule, testName, thisTest)

addTests()
