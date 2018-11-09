import math
import random
from iid.basic import minMax, logAdd, logSub

def same(a, b, eps):
    if a != 0:
        return abs((a - b)/a) < eps
    return abs(a - b) < eps

def test_logAdd_1():
    random.seed(17)
    N = 1000
    for i in range(N):
        a = random.expovariate(0.2)
        la = math.log(a)
        b = random.expovariate(0.2)
        lb = math.log(b)
        c = a + b
        lc = logAdd(la, lb)
        c0 = math.exp(lc)
        assert same(c, c0, 1e-14)

def test_logAdd_2():
    random.seed(18)
    N = 1000
    for i in range(N):
        a = random.expovariate(0.0002)
        la = math.log(a)
        b = random.expovariate(0.2)
        lb = math.log(b)
        c = a + b
        lc = logAdd(la, lb)
        c0 = math.exp(lc)
        assert same(c, c0, 1e-14)

def test_logSub_1():
    random.seed(18)
    N = 1000
    for i in range(N):
        a0 = random.expovariate(0.2)
        b0 = random.expovariate(0.2)
        (b, a) = minMax(a0, b0)
        la = math.log(a)
        lb = math.log(b)
        c = a - b
        lc = logSub(la, lb)
        c0 = math.exp(lc)
        assert same(c, c0, 1e-13)

def test_logSub_2():
    random.seed(18)
    N = 1000
    for i in range(N):
        a0 = random.expovariate(0.0002)
        b0 = random.expovariate(0.2)
        (b, a) = minMax(a0, b0)
        la = math.log(a)
        lb = math.log(b)
        c = a - b
        lc = logSub(la, lb)
        c0 = math.exp(lc)
        assert same(c, c0, 1e-13)

