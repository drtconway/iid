import os
import sys
import yaml

from iid.special import erf

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def addTests():
    thisModule = sys.modules[__name__]
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'special', 'erf.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        x = itm[0]
        r = itm[1]
        testName = 'test_erf_%f' % (x, )
        def thisTest():
            r0 = erf(x)
            assert same(r, r0, 1e-12)
        setattr(thisModule, testName, thisTest)

addTests()
