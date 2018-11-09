from iid.special import gammaP, gammaQ, gammaPQ
from iid.tests.data_igamma import data

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_gammaP():
    for (a, x, rL, lrL, rU, lrU) in data:
        rL0 = gammaP(a, x)
        assert same(rL, rL0, 2e-14)

def test_gammaQ():
    for (a, x, rL, lrL, rU, lrU) in data:
        rU0 = gammaQ(a, x)
        assert same(rU, rU0, 2e-14)

def test_gammaPQ():
    for (a, x, rL, lrL, rU, lrU) in data:
        (rL0, rU0) = gammaPQ(a, x)
        assert same(rL, rL0, 2e-14)
        assert same(rU, rU0, 2e-14)

