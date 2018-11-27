from iid.special import erfc, logErfc
from iid.tests.data_erfc import data

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_erfc():
    for (x, r, lr) in data:
        r0 = erfc(x)
        assert same(r, r0, 1e-13)

def test_logErfc():
    for (x, r, lr) in data:
        lr0 = logErfc(x)
        assert same(lr, lr0, 5e-13)

