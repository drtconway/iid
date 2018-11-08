from iid.basic import gamma, logGamma
from iid.tests.data_gamma import data

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_gamma():
    for (x, r, lr) in data:
        r0 = gamma(x)
        assert same(r, r0, 1e-13)

def test_logGamma():
    for (x, r, lr) in data:
        lr0 = logGamma(x)
        assert same(lr, lr0, 1e-13)

