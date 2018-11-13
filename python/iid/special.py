import math
import iid.basic as basic

def pickGammaMethod(a, x):
    '''choose between different methods for computing gamma.'''
    if x > 0.25:
        alpha = x + 0.25
    else:
        alpha = math.log(0.5) / math.log(x)

    if a > alpha:
        return 'g'

    x0 = 1.5
    if x < x0 and -0.5 <= a:
        return 'G1'
    if x < x0 and a < -0.5:
        assert False, 'incomplete gamma for a < -0.5 not implemented'
        #return 'G2'
    return 'G3'

def upperGammaKummer(a, x):
    '''compute upper gamma using Kummer's series.'''
    u = (basic.gamma(1 + a) - 1) / a - (math.pow(x, a) - 1) / a

    k = 1
    tk = 1
    ss = tk
    while True:
        tk *= - (a + k)*x / ((a + k + 1)*(k + 1))
        ss += tk
        if abs(tk/ss) < 1e-16:
            break
        k += 1
    v = (math.pow(x, a + 1) / (a + 1)) * ss

    return u + v

def logUpperGammaKummer(a, x):
    '''compute log of upper gamma using Kummer's series.'''
    la = math.log(a)
    lx = math.log(x)

    lga = basic.logGamma(a)
    lxaoa = a * lx - la
    if lga > lxaoa:
        flip = False
        lu = basic.logSub(lga, lxaoa)
    else:
        flip = True
        lu = basic.logSub(lxaoa, lga)

    ltk = 0
    tkPos = True
    lss = ltk
    k = 1
    while True:
        tkPos = not tkPos
        ltk += math.log(a + k) + lx - (math.log(a + k + 1) + math.log(k + 1))
        if tkPos:
            lss = basic.logAdd(lss, ltk)
        else:
            lss = basic.logSub(lss, ltk)
        if ltk - lss < -45:
            break
        k += 1
    lv = (a + 1)*lx - math.log(a + 1) + lss

    if not flip:
        return basic.logAdd(lu, lv)
    else:
        return basic.logSub(lv, lu)

def upperGammaLegendre(a, x):
    '''compute upper gamma using a series derived from Legendre's continued fraction.'''
    pk = 0
    qk = (x - 1 - a)*(x + 1 - a)
    rk = 4*(x + 1 - a)
    sk = 1 - a
    rhok = 0
    tk = 1
    k = 1
    ss = tk
    while True:
        pk += sk
        qk += rk
        rk += 8
        sk += 2
        tauk = pk*(1 + rhok)
        rhok = tauk / (qk - tauk)
        tk *= rhok
        ss += tk
        if abs(tk/ss) < 1e-16:
            break
        k += 1
    return math.exp(-x)*math.pow(x, a)*ss / (x + 1 - a)

def logUpperGammaLegendre(a, x):
    '''compute log of upper gamma using a series derived from Legendre's continued fraction.'''
    pk = 0
    qk = (x - 1 - a)*(x + 1 - a)
    rk = 4*(x + 1 - a)
    sk = 1 - a
    rhok = 0
    tk = 1
    k = 1
    ss = tk
    while True:
        pk += sk
        qk += rk
        rk += 8
        sk += 2
        tauk = pk*(1 + rhok)
        rhok = tauk / (qk - tauk)
        tk *= rhok
        ss += tk
        if abs(tk/ss) < 1e-16:
            break
        k += 1
    return a * math.log(x) - x + math.log(ss) - math.log(x + 1 - a)

def lowerRegularizedGammaSeries(a, x):
    '''compute lower regularized gamma with a Taylor series.'''
    ss =  1.0 / basic.gamma(a + 1)
    xn = 1
    n = 1
    while True:
        xn *= x
        tn = xn / basic.gamma(a + n + 1)
        ss += tn
        if abs(tn/ss) < 1e-16:
            break
        n += 1
    return math.exp(-x)*math.pow(x, a)*ss

def logLowerRegularizedGammaSeries(a, x):
    '''compute lower regularized gamma with a Taylor series.'''
    lx = math.log(x)
    lss =  - basic.logGamma(a + 1)
    lxn = 0
    n = 1
    while True:
        lxn += lx
        ltn = lxn - basic.logGamma(a + n + 1)
        lss = basic.logAdd(lss, ltn)
        if ltn - lss < -45:
            break
        n += 1
    return a*lx - x + lss

def gammaP(a, x):
    '''compute lower regularised gamma.'''
    alg = pickGammaMethod(a, x)
    if alg == 'g':
        return lowerRegularizedGammaSeries(a, x)

    if alg == 'G1':
        Gam = upperGammaKummer(a, x)
    if alg == 'G3':
        Gam = upperGammaLegendre(a, x)

    ga = basic.gamma(a)
    bigG = Gam / ga
    gs = 1 - bigG
    return gs

def logGammaP(a, x):
    '''compute the log of lower regularised gamma.'''
    alg = pickGammaMethod(a, x)
    if alg == 'g':
        return logLowerRegularizedGammaSeries(a, x)

    if alg == 'G1':
        Gam = upperGammaKummer(a, x)
        lGam = logUpperGammaKummer(a, x)
    if alg == 'G3':
        lGam = logUpperGammaLegendre(a, x)

    lga = basic.logGamma(a)
    lBigG = lGam - lga
    lgs = basic.log1mexp(lBigG)
    return lgs

def gammaQ(a, x):
    '''compute upper regularised gamma.'''
    alg = pickGammaMethod(a, x)

    if alg == 'g':
        gs = lowerRegularizedGammaSeries(a, x)
        return 1 - gs

    if alg == 'G1':
        Gam = upperGammaKummer(a, x)
    if alg == 'G3':
        Gam = upperGammaLegendre(a, x)

    return Gam / basic.gamma(a)

def logGammaQ(a, x):
    '''compute log upper regularised gamma.'''
    alg = pickGammaMethod(a, x)

    if alg == 'g':
        lgs = logLowerRegularizedGammaSeries(a, x)
        return basic.log1mexp(lgs)

    if alg == 'G1':
        lGam = logUpperGammaKummer(a, x)
    if alg == 'G3':
        lGam = logUpperGammaLegendre(a, x)

    return lGam - basic.logGamma(a)

def gammaPQ(a, x):
    '''compute both lower and upper regularised gamma.'''
    alg = pickGammaMethod(a, x)

    if alg == 'g':
        gs = lowerRegularizedGammaSeries(a, x)
        return (gs, 1 - gs)

    if alg == 'G1':
        Gam = upperGammaKummer(a, x)
    if alg == 'G3':
        Gam = upperGammaLegendre(a, x)

    bigG = Gam / basic.gamma(a)
    return (1 - bigG, bigG)

def logGammaPQ(a, x):
    '''compute both log lower and log upper regularised gamma.'''
    alg = pickGammaMethod(a, x)

    if alg == 'g':
        lgs = logLowerRegularizedGammaSeries(a, x)
        return (lgs, basic.log1mexp(lgs))

    if alg == 'G1':
        lGam = logUpperGammaKummer(a, x)
    if alg == 'G3':
        lGam = logUpperGammaLegendre(a, x)

    lBigG = lGam - basic.logGamma(a)
    return (basic.log1mexp(lBigG), lBigG)

def logBeta(a, b):
    '''compute log beta'''
    return basic.logGamma(a) + basic.logGamma(b) - basic.logGamma(a + b)

def beta(a, b):
    '''compute beta'''
    return math.exp(logBeta(a, b))

def lowerBetaLin(a, b, x):
    '''compute lower incomplete beta using a series approximation for integer a & b'''
    #assert type(a) == int
    #assert type(b) == int
    a = int(a)
    b = int(b)
    y = 1 - x
    m = a
    n = a + b - 1
    s = 0
    for j in xrange(m, n+1):
        c = basic.choose(n, j)
        t = c * (x**j) * (y**(n-j))
        s += t
    return s

def logLowerBetaLog(a, b, x):
    '''compute log lower incomplete beta using a log series approximation for integer a & b'''
    #assert type(a) == int
    #assert type(b) == int
    a = int(a)
    b = int(b)
    y = 1 - x
    lx = math.log(x)
    ly = basic.log1p(-x)
    m = a
    n = a + b - 1
    s = None
    for j in xrange(m, n+1):
        lc = basic.logChoose(n, j)
        t = lc  + j*lx + (n-j)*ly
        if s is None:
            s = t
        else:
            s = basic.logAdd(s, t)
    if s is None:
        # ~log(1/inf)
        return -999
    return s

def lowerBetaLog(a, b, x):
    '''compute lower incomplete beta using a log series approximation for integer a & b'''
    return math.exp(logLowerBetaLog(a, b, x))

def lowerBetaCont(a, b, x):
    '''compute lower incomplete beta using a continued fraction approximation'''
    #assert a > 1
    #assert b > 1
    y = 1 - x
    lx = math.log(x)
    ly = basic.log1p(-x)

    def aa(j):
        if j == 1:
            return 1

        m = j - 1
        return ((a + m - 1)*(a + b + m - 1)*m*(b - m)*(x**2)) / ((a + 2*m - 1)**2)

    def bb(j):
        m = j - 1
        return m + (m*(b-m)*x)/(a + 2*m - 1) + ((a + m)*(a - (a + b)*x + 1 + m*(2 - x)))/(a + 2*m + 1)

    return (x**a * y**b / beta(a, b)) * basic.contFrac(aa, bb)

def txLowerBeta(f, a, b, x, n):
    '''apply a duplication formula to compute Ix(a, b) in terms of Ix(a + n, b)'''
    y = 1 - x
    lx = math.log(x)
    ly = basic.log1p(-x)
    lgb = basic.logGamma(b)
    p = f(a + n, b, x)
    ls = None
    for j in xrange(1, n + 1):
        lt = basic.logGamma(a + b + j - 1) + (j - 1)*lx - (lgb + basic.logGamma(a + j))
        if ls is None:
            ls = lt
        else:
            ls = basic.logAdd(ls, lt)
    s = math.exp(a*lx + b*ly + ls)
    return p + s

def lowerBeta(a, b, x):
    '''compute lower incomplete beta Ix(a, b) heuristics to choose between different approximations'''
    if a < b:
        return 1.0 - lowerBeta(b, a, 1 - x)

    if math.floor(a)/a > 0.95 and math.floor(b)/b > 0.95:
        v = basic.logChoose(a + b, b)
        if v < 100:
            return lowerBetaLin(a, b, x)
        else:
            return math.exp(logLowerBetaLog(a, b, x))
    else:
        return lowerBetaCont(a, b, x)

def erf(x):
    '''compute the error function erf(x)'''
    if x < 0:
        return - erf(-x)

    if x > 5:
        return 1 - math.exp(-x*x)/basic.sqrt_pi

    if x > 4:
        return 1 - math.exp(-x*x)/(basic.sqrt_pi*x)

    x = float(x)
    s = x
    n = 1
    p = 1
    while True:
       p *= -x*x / n
       t = x / (2*n + 1) * p
       s += t
       if s == 0 or abs(t/s) < 1e-12:
           break
       n += 1
    return 2.0/basic.sqrt_pi * s

def erfc(x):
    '''compute the complement error function erfc(x) == 1 - erf(x)'''
    if x < - 0.5:
        return 2.0 - erfc(-x)
    if x < 0.5:
        return 1 - erf(x)
    aa = lambda j : 0.5*j
    bb = lambda j : x
    v = basic.contFrac(aa, bb)
    return math.exp(-x*x)/(basic.sqrt_pi*(x + v))

def logErfc(x):
    '''compute log(erfc(x))'''
    if x < - 0.5:
        return math.log(2.0 - erfc(-x))
    if x < 0.5:
        return basic.log1p(-erf(x))
    aa = lambda j : 0.5*j
    bb = lambda j : x
    v = basic.contFrac(aa, bb)
    return -x*x - basic.log_sqrt_pi - math.log(x + v)

