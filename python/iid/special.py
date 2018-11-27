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

def upperGammaSeries(a, x):
    s = 0
    k = 1
    m1k = 1
    xk = 1
    kfac = 1
    ts = []
    while True:
        m1k *= -1
        xk *= x
        kfac *= k
        t = m1k*xk/((a+k)*kfac)
        s += t
        ts.append(t)
        if abs(t/s) < 1e-14:
            break
        k += 1

    s1 = basic.kahanSum(ts)
    gap1m1 = basic.gamma(a + 1) - 1
    return (gap1m1 - basic.powm1(x, a)) / a + math.pow(x, a)*s

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

def lowerGammaSeries(a, x):
    t = 1.0
    s = 0.0
    while True:
        a += 1
        t *= x/a
        s += t
        if abs(t/s) < 1e-14:
            break
    return s

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

def digamma(x):
    '''compute the digamma function of x.'''
    g = 0.5772156649015328606065120900824024310421
    n = 1
    s = 0
    while True:
        t = 1.0 / (n + x) - 1.0/n
        s += t
        if abs(t/s) < 1e-14:
            break
        n += 1
    return s

def logBeta(a, b):
    '''compute log beta'''
    return basic.logGamma(a) + basic.logGamma(b) - basic.logGamma(a + b)

def beta(a, b):
    '''compute beta'''
    return math.exp(logBeta(a, b))

def betaIntImpl(a, b):
    assert type(a) == int and a > 0
    assert type(b) == int and b > 0

    if a < b:
        return betaIntImpl(b, a)

    n = 1
    d = a
    for i in range(1,b):
        n *= i
        d *= (a + i)
    f = basic.gcd(n, d)
    n /= f
    d /= f
    ln = math.log(n)
    ld = math.log(d)
    if ln > 700 or ld > 700:
        return math.exp(ln - ld)
    return float(n)/float(d)

def logBetaInt(a, b):
    lga = basic.logGamma(a)
    lgb = basic.logGamma(b)
    lgab = basic.logGamma(a+b)
    lb = lga + lgb - lgab
    if lb > -700:
        return math.log(betaIntImpl(a, b))
    return lb

def betaInt(a, b):
    lga = basic.logGamma(a)
    lgb = basic.logGamma(b)
    lgab = basic.logGamma(a+b)
    lb = lga + lgb - lgab
    if lb > -700:
        return betaIntImpl(a, b)
    return math.exp(lb)

def lowerBetaInt(a, b, x):
    y = 1 - x
    if x > float(a)/float(a+b):
        return 1 - lowerBetaInt(b, a, y)

    lx = math.log(x)
    ly = basic.log1p(-x)
    lpfx = a*lx + b*ly

    if a < 350 and b < 350 and lpfx > -700:
        xa = math.pow(x, a)
        yb = math.pow(y, b)
        bx = xa*yb/a * basic.hyper([a+b, 1], [a+1], x)
        bab = betaInt(a, b)
        return bx/bab
    else:
        # Try log space
        lbx = lpfx - math.log(a) + basic.logHyper([a+b, 1], [a+1], x)
        lbab = logBetaInt(a, b)
        if lbx - lbab > 0:
            # avoid rounding errors
            lbab = lbx
        return math.exp(lbx - lbab)

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

def logLowerBetaSeries(a, b, x):
    lx = math.log(x)
    lv = basic.logGamma(a+b) - basic.logGamma(a + 1) - basic.logGamma(b) + a * lx

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
    lu = basic.log1p(a*s)
    return lv + lu

def lowerBetaSeries(a, b, x):
    return math.exp(logLowerBetaSeries(a, b, x))

def lowerBetaOffset(a, b, x, y, n):
    lx = math.log(x)
    ly = math.log(y)
    lv = a*lx + b*ly

    ls = None
    ldi = 0
    for j in range(1, n+1):
        lt = basic.logGamma(a + b + j - 1) - basic.logGamma(b) - basic.logGamma(a + j) + (j-1)*lx
        if ls is None:
            ls = lt
        else:
            ls = basic.logAdd(ls, lt)
    return math.exp(lv+ls)

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

def lowerBetaClass(a, b, x):
    y = 1 - x
    p = a/(a+b)
    q = b/(a+b)
    if type(a) == int and type(b) == int:
        return lowerBetaInt(a, b, x)
    if min(a, b) <= 1:
        if x > 0.5:
            r = lowerBetaClass(b, a, 1 - x)
            if type(r) == str:
                return '!'+r
            return 1 - r

        if any([max(a, b) <= 1 and a > min(0.2, b),
                max(a, b) < 1 and a < min(0.2, b) and math.pow(x, a) < 0.9,
                max(a, b) > 1 and b < 1,
                max(a, b) > 1 and b > 1 and x < 0.1, math.pow(b*x, a) < 0.7]):
            return lowerBetaSeries(a, b, x)
        if any([max(a, b) > 1 and b > 1 and x < 0.1, math.pow(b*x, a) < 0.7,
                max(a, b) > 1 and b > 1 and x > 0.3]):
            return 1 - lowerBetaSeries(b, a, y)
        if any([max(a, b) > 1 and b > 15 and 0.1 < x and x < 0.3,
                max(a, b) > 1 and b > 15 and x < 0.1 and math.pow(b*x, a) > 0.7]):
            return  1 - lowerBetaSeries(b, a, y)
        if any([max(a, b) > 1 and b > 1 and 0.1 < x and x < 0.3 and b < 15,
                max(a, b) > 1 and b > 1 and x < 0.1 and math.pow(b*x, a) > 0.7 and b <= 15,
                max(a, b) < 1 and a < min(0.2, b), math.pow(x, a) > 0.9 and x <0.3]):
            n = 20
            u = lowerBetaOffset(b + n, a, y, x, n)
            v = lowerBetaSeries(b + n, a, y)
            return 1 - u + v
    else:
        if x > p:
            r = lowerBetaClass(b, a, 1 - x)
            if type(r) == str:
                return '!'+r
            return 1 - r

        if any([b < 40 and b*x < 0.7]):
            return lowerBetaSeries(a, b, x)
        if any([b < 40 and b*x > 0.7 and x <= 0.7]):
            n = int(math.floor(b))
            if n == b:
                n -= 1
            u = lowerBetaOffset(b - n, a, y, x, n)
            v = lowerBetaSeries(a, b - n, x)
            return u + v
        if any([b < 40 and x > 0.7 and a > 15]):
            #return 'bgrat(a+m, b, x, y, w0 = bup(b\', a, y, x, n) + bup(a, b\', x, y, m)), m = 20'
            return lowerBetaSeries(a, b, x)
        if any([b < 40 and x > 0.7 and a <= 15]):
            return lowerBetaSeries(a, b, x)
        if any([b >= 40 and a <= b and a <= 100,
                b >= 40 and 100 < a and a <= b and x < 0.97*p,
                b >= 40 and a > b and b <= 100,
                b >= 40 and 100 < a and b < a and y > 1.03*q]):
            return lowerBetaCont(a, b, x)
        if any([b >= 40 and 100 < a and a <= b and x >= 0.97*p,
                b >= 40 and 100 < b and b < a and y <= 1.03*q]):
            #return 'basym(a, b, x, y)'
            return lowerBetaSeries(a, b, x)

    return 'none'

def lowerBeta(a, b, x):
    '''compute lower incomplete beta Ix(a, b) heuristics to choose between different approximations'''
    return lowerBetaClass(a, b, x)

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
    if x < 50:
        aa = lambda j : 0.5*j
        bb = lambda j : x
        v = basic.contFrac(aa, bb)
        return math.exp(-x*x)/(basic.sqrt_pi*(x + v))
    return gammaQ(0.5, x*x)/basic.sqrt_pi

def logErfc(x):
    '''compute log(erfc(x))'''
    if x < - 0.5:
        return math.log(2.0 - erfc(-x))
    if x < 0.5:
        return basic.log1p(-erf(x))
    if x < 50:
        aa = lambda j : 0.5*j
        bb = lambda j : x
        v = basic.contFrac(aa, bb)
        return -x*x - basic.log_sqrt_pi - math.log(x + v)
    return logGammaQ(0.5, x*x) - basic.log_sqrt_pi

