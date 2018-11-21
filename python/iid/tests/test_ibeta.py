import os
import yaml
from iid.special import lowerBeta
import iid.tests.data_ibetaInt as ibI

def sameAbs(a, b, eps):
    return abs(a - b) < eps

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_ibetaInt():
    for itm in ibI.data:
        a = itm[0]
        b = itm[1]
        x = itm[2]
        rL = itm[3]
        rL0 = lowerBeta(a, b, x)
        assert same(rL, rL0, 1e-12)

def test_ibetaGood():
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data-ibeta-good.yaml')) as f:
        data = yaml.load(f)
    for itm in data:
        a = itm[0]
        b = itm[1]
        x = itm[2]
        rL = itm[3]
        rL0 = lowerBeta(a, b, x)
        assert same(rL, rL0, 1e-12)

def test_ibetaBad():
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data-ibeta-bad.yaml')) as f:
        data = yaml.load(f)
    for itm in data:
        a = itm[0]
        b = itm[1]
        x = itm[2]
        rL = itm[3]
        rL0 = lowerBeta(a, b, x)
        assert not same(rL, rL0, 1e-12)

