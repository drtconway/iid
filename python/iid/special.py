import math
import iid.basic as basic

def logLowerGamma(s, x):
    '''compute log lower incomplete gamma'''
    lx = math.log(x)
    lgs = basic.logGamma(s)
    g = lgs + math.log(s)
    w = -g
    k = 0
    while True:
        k += 1
        g += math.log(s + k)
        t = k * lx - g
        w = basic.logAdd(w, t)
        if t - w < -45:
            break
    return s*lx + lgs - x + w

def lowerGamma(s, x):
    '''compute lower incomplete gamma'''
    return math.exp(logLowerGamma(s, x))

def logUpperGammaInner(s, x):
    if x == 0:
        M, N = 1, basic.gamma(s)

    def aa1(n):
        if n == 1:
            return 1
        return -(n - 1)*(n - s - 1)

    aa2 = lambda n: aa1(n + 1)

    bb1 = lambda n: x + 2*n - 1 - s
    bb2 = lambda n: bb1(n + 1)

    if x == s - 1:
        M = aa1(1) / basic.contFrac(aa2, bb2)
    else:
        M = basic.contFrac(aa1, bb1)
    N = s*math.log(x) - x
    return M, N

def logUpperGamma(s, x):
    '''compute log upper incomplete gamma'''
    M, N = logUpperGammaInner(s, x)
    return N + math.log(M)

def upperGamma(s, x):
    '''compute upper incomplete gamma'''
    M, N = logUpperGammaInner(s, x)
    return M * math.exp(N)

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
        return 1 - math.exp(-x*x)/math.sqrt(math.pi)

    if x > 4:
        return 1 - math.exp(-x*x)/(math.sqrt(math.pi)*x)

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
    return math.exp(-x*x)/math.sqrt(pi)/(x + v)

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

