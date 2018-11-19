from iid.special import betaInt, logBetaInt
from iid.tests.data_betaInt import data

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_betaInt():
    for (a, b, r, lr) in data:
        r0 = betaInt(a, b)
        assert same(r, r0, 1e-14)

def test_logBetaInt():
    for (a, b, r, lr) in data:
        lr0 = logBetaInt(a, b)
        assert same(lr, lr0, 1e-14)

