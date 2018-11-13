import math
from iid.basic import logGamma, log1p
from iid.special import lowerBeta
from iid.tests.data_ibeta import data

def sameAbs(a, b, eps):
    return abs(a - b) < eps

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def logLowerBetaSeries(a, b, x):
    lx = math.log(x)
    lv = logGamma(a+b) - logGamma(a + 1) - logGamma(b) + a * lx

    s = 0
    mbPoc = 1
    nFac = 1
    xN = 1
    n = 1
    while True:
        mbPoc *= n - b
        nFac *= n
        xN *= x
        t = (mbPoc / (nFac * (a + n))) * xN
        s += t
        if abs(t/s) < 1e-14:
            break
        n += 1
    lu = log1p(a*s)
    return lv + lu

def lowerBetaSeries(a, b, x):
    return math.exp(logLowerBetaSeries(a, b, x))

def lowerBetaClass(a, b, x):
    y = 1 - x
    p = a/(a+b)
    q = b/(a+b)
    if min(a, b) <= 1:
        if any([max(a, b) <= 1 and a > min(0.2, b),
                max(a, b) < 1 and a < min(0.2, b) and math.pow(x, a) < 0.9,
                max(a, b) > 1 and b < 1,
                max(a, b) > 1 and b > 1 and x < 0.1, math.pow(b*x, a) < 0.7]):
            yield 'ser(a, b, x)'
        if any([max(a, b) > 1 and b > 1 and x < 0.1, math.pow(b*x, a) < 0.7,
                max(a, b) > 1 and b > 1 and x > 0.3]):
            yield 'ser(b, a, y)'
        if any([max(a, b) > 1 and b > 15 and 0.1 < x and x < 0.3,
                max(a, b) > 1 and b > 15 and x < 0.1 and math.pow(b*x, a) > 0.7]):
            yield 'bgrat(b, a, y, x, w0=0)'
        if any([max(a, b) > 1 and b > 1 and 0.1 < x and x < 0.3 and b < 15,
                max(a, b) > 1 and b > 1 and x < 0.1 and math.pow(b*x, a) > 0.7 and b <= 15,
                max(a, b) < 1 and a < min(0.2, b), math.pow(x, a) > 0.9 and x <0.3]):
            yield 'bgrat(b, a, y, x, w0=bup(b, a, y, x, n = 20))'
    else:
        if any([b < 40 and b*x < 0.7]):
            yield 'ser(a, b, x)'
        if any([b < 40 and b*x > 0.7 and x <= 0.7]):
            yield 'bup(b\', a, y, x, n) + ser(a, b, x)'
        if any([b < 40 and x > 0.7 and a > 15]):
            yield 'bgrat(a+m, b, x, y, w0 = bup(b\', a, y, x, n) + bup(a, b\', x, y, m)), m = 20'
        if any([b >= 40 and a <= b and a <= 15,
                b >= 40 and 100 < a and a <= b and x < 0.97*p,
                b >= 40 and a > b and b <= 100,
                b >= 40 and 100 < a and a < b and y > 1.03*q]):
            yield 'bfrac(a, b, x, y)'
        if any([b >= 40 and 100 < a and a <= b and x >= 0.97*p,
                b >= 40 and 100 < b and b < a and y <= 1.03*q]):
            yield 'basym(a, b, x, y)'


for itm in data:
    a = itm[0]
    b = itm[1]
    x = itm[2]
    cls = list(lowerBetaClass(a, b, x))
    if len(cls) == 0:
        cls = ['!' + cx for cx in lowerBetaClass(b, a, 1 - x)]
    print a, b, x, len(cls), cls

