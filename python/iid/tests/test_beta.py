from iid.special import beta, logBeta
from iid.tests.data_beta import data

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_beta():
    for (a, b, r, lr) in data:
        r0 = beta(a, b)
        assert same(r, r0, 1e-12)

def test_logBeta():
    for (a, b, r, lr) in data:
        lr0 = logBeta(a, b)
        assert same(lr, lr0, 1e-12)

