import os
import sys
import yaml

from iid.special import betaInt, logBetaInt

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def addTests():
    thisModule = sys.modules[__name__]
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'special', 'betaInt.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        a = itm[0]
        b = itm[1]
        r = itm[2]
        lr = itm[3]
        testName = 'test_betaInt_%f_%f' % (a, b)
        def thisTest():
            r0 = betaInt(a, b)
            assert same(r, r0, 1e-12)
            lr0 = logBetaInt(a, b)
            assert same(lr, lr0, 1e-12)
        setattr(thisModule, testName, thisTest)

addTests()
