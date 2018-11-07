from iid.special import erf
from iid.tests.data_erf import data

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_erf():
    for (x, r) in data:
        r0 = erf(x)
        print x, r, r0, abs((r - r0)/r)
        assert same(r, r0, 5e-10)
