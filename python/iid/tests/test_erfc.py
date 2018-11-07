from iid.special import erfc, logErfc
from iid.tests.data_erfc import data

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_erfc():
    for (x, r, lr) in data:
        r0 = erfc(x)
        if r == 0:
            print x, r, r0, abs((r - r0))
        else:
            print x, r, r0, abs((r - r0)/r)
        assert same(r, r0, 1e-13)

def test_logErfc():
    for (x, r, lr) in data:
        lr0 = logErfc(x)
        if lr0 == 0:
            print x, lr, lr0, abs((lr - lr0))
        else:
            print x, lr, lr0, abs((lr - lr0)/lr)
        assert same(lr, lr0, 5e-14)

