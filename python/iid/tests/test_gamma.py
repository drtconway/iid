import os
import sys
import yaml
from iid.basic import gamma, logGamma

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def addTests():
    thisModule = sys.modules[__name__]
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'special', 'gamma.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        x = itm[0]
        r = itm[1]
        lr = itm[2]
        testName = 'test_gamma_%f' % (x, )
        def thisTest():
            r0 = gamma(x)
            assert same(r, r0, 1e-12)
            lr0 = logGamma(x)
            assert same(lr, lr0, 1e-12)
        setattr(thisModule, testName, thisTest)

addTests()
