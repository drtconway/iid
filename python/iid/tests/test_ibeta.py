from iid.special import lowerBeta
import iid.tests.data_ibeta as ib
import iid.tests.data_ibetaInt as ibI

def sameAbs(a, b, eps):
    return abs(a - b) < eps

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_ibetaInt():
    for itm in ibI.data:
        a = itm[0]
        b = itm[1]
        x = itm[2]
        rL = itm[3]
        rL0 = lowerBeta(a, b, x)
        if x > 0.001 and x < 0.99:
            if not same(rL, rL0, 5e-11):
                print a, b, x, rL, rL0, abs((rL - rL0)/rL)
            assert same(rL, rL0, 5e-11)
        else:
            if not same(rL, rL0, 5e-9):
                print a, b, x, rL, rL0, abs((rL - rL0)/rL)
            assert same(rL, rL0, 5e-9)
