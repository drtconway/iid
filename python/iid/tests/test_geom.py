import os
import yaml
from iid.geom import geom
import iid.tests.data_ibetaInt as ibI

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_binom():
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
        assert same(m, dst.mean(), 1e-14)
        assert same(v, dst.var(), 1e-14)
        rP0 = dst.pmf(k)
        if not same(rP, rP0, 1e-14):
            print(itm)
        assert same(rP, rP0, 1e-14)
        rL0 = dst.cdf(k)
        assert same(rL, rL0, 1e-14)


