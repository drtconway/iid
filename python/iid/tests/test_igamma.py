from iid.special import gammaP, gammaQ, gammaPQ, logGammaP, logGammaQ, logGammaPQ
from iid.special import pickGammaMethod
from iid.tests.data_igamma import data

def sameAbs(a, b, eps):
    return abs(a - b) < eps

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

def test_logGammaP():
    for (a, x, rL, lrL, rU, lrU) in data:
        lrL0 = logGammaP(a, x)
        assert same(lrL, lrL0, 4e-14)

def test_logGammaQ():
    for (a, x, rL, lrL, rU, lrU) in data:
        lrU0 = logGammaQ(a, x)
        if lrL < -85:
            assert sameAbs(lrU, lrU0, 4e-16)
            continue
        assert same(lrU, lrU0, 1e-13)

def test_logGammaPQ():
    for (a, x, rL, lrL, rU, lrU) in data:
        (lrL0, lrU0) = logGammaPQ(a, x)
        assert same(lrL, lrL0, 4e-14)
        if lrL < -85:
            assert sameAbs(lrU, lrU0, 4e-16)
            continue
        assert same(lrU, lrU0, 1e-13)

