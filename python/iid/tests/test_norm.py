import math
import os
import yaml
from iid.norm import norm

def sameAbs(a, b, eps):
    return abs(a - b) < eps

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_norm():
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'dist', 'norm.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        mu = itm[0]
        sig = itm[1]
        x = itm[2]
        m = itm[3]
        v = itm[4]
        rP = itm[5]
        rL = itm[6]
        dst = norm(mu, sig)
        assert same(m, dst.mean(), 1e-14)
        assert same(v, dst.var(), 1e-14)

        if rP == 0 or rP > 1e-250:
            rP0 = dst.pdf(x)
            if not same(rP, rP0, 1e-12):
                print(mu, sig, x, rP, rP0, abs((rP - rP0)/rP))
            assert same(rP, rP0, 1e-12)
        else:
            lrP = math.log(rP)
            eps = 1e-12
            if lrP < -720:
                eps = 1e-4
            lrP0 = dst.pdf(x, log=True)
            if not same(lrP, lrP0, eps):
                print(mu, sig, x, lrP, lrP0, abs((lrP - lrP0)/lrP))
            assert same(lrP, lrP0, eps)
        if rL == 0 or rL > 1e-250:
            rL0 = dst.cdf(x)
            if not same(rL, rL0, 1e-12):
                print(mu, sig, x, rL, rL0, abs((rL - rL0)/rL))
            assert same(rL, rL0, 1e-12)
        else:
            lrL = math.log(rL)
            eps = 1e-12
            if lrL < -720:
                eps = 1e-4
            lrL0 = dst.cdf(x, log=True)
            if not same(lrL, lrL0, eps):
                print(mu, sig, x, lrL, lrL0, abs((lrL - lrL0)/lrL))
            assert same(lrL, lrL0, eps)

