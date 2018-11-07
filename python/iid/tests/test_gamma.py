from iid.basic import gamma, logGamma
from iid.tests.data_gamma import data

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_gamma():
    for (x, r, lr) in data:
        r0 = gamma(x)
        if r == 0:
            print x, r, r0, abs((r - r0))
        else:
            print x, r, r0, abs((r - r0)/r)
        assert same(r, r0, 1e-13)

def test_logGamma():
    for (x, r, lr) in data:
        lr0 = logGamma(x)
        if lr0 == 0:
            print x, lr, lr0, abs((lr - lr0))
        else:
            print x, lr, lr0, abs((lr - lr0)/lr)
        assert same(lr, lr0, 5e-14)

