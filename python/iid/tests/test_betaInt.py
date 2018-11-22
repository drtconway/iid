import os
import yaml

from iid.special import betaInt, logBetaInt

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

wd = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(wd, 'data', 'special', 'betaInt.yaml')) as f:
    data = yaml.load(f)['data']

def test_betaInt():
    for (a, b, r, lr) in data:
        r0 = betaInt(a, b)
        assert same(r, r0, 1e-14)

def test_logBetaInt():
    for (a, b, r, lr) in data:
        lr0 = logBetaInt(a, b)
        assert same(lr, lr0, 1e-14)

