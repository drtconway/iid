import os
import yaml
from iid.hyper import hyper

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_binom():
    wd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(wd, 'data', 'dist', 'hyper.yaml')) as f:
        data = yaml.load(f)['data']
    for itm in data:
        N = itm[0]
        K = itm[1]
        n = itm[2]
        k = itm[3]
        m = itm[4]
        v = itm[5]
        rP = itm[6]
        rL = itm[7]
        rU = itm[8]
        dst = hyper(N, K, n)
        assert same(m, dst.mean(), 1e-14)
        assert same(v, dst.var(), 1e-14)
        rP0 = dst.pmf(k)
        if not same(rP, rP0, 1e-14):
            print(itm, rP0, abs((rP-rP0)/rP))
        assert same(rP, rP0, 1e-14)
        rL0 = dst.cdf(k)
        if not same(rL, rL0, 1e-14):
            print(itm, rL0, abs((rL-rL0)/rL))
        assert same(rL, rL0, 1e-14)

