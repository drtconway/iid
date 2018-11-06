import math

'''

Sources for the methods include:
- https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
- https://en.wikipedia.org/wiki/Taylor_series

'''

eps = 1e-16
'''float: machine precision such that (1 + eps) == 1.0.

The value 1e-16 is computed with Python 2.7 on x86_64.
'''

pi = 3.14159265358979323846264338327950288
'''float: high precision value for pi'''

log_pi = math.log(pi)
'''float: log(pi)'''

log_2 = math.log(2)
'''float: log(2)'''

sqrt_pi = math.sqrt(pi)
'''float: sqrt(pi)'''

log_sqrt_pi = math.log(sqrt_pi)
'''float: log(sqrt(pi))'''

def minMax(a, b):
    '''select the minimum and maximum of two values.
    
    Args:
        a: a comparable value
        b: a comparable value
        
    Returns:
        tuple: (min(a, b), max(a, b))
    '''
    if a >= b:
        return (b, a)
    else:
        return (a, b)

def contFrac(aa, bb):
    '''Compute a continued fraction.

    The literature on continued fractions is a bit inconsistent:
    in this implementation the a_i terms are used in the numerator
    and the b_i terms in the denominator.

    Args:
        aa (function): numerator terms represented by a function (int -> float) that returns a_i for i >= 1
        bb (function): denominator terms represented by a function (int -> float) that returns b_i for i >= 1

    Returns:
        The value of the continued fraction to the machine precision.
    '''
    dm = 1e-300
    a1 = aa(1)
    b1 = bb(1)
    f = aa(1)/bb(1)
    C = a1/dm
    D = 1/b1
    n = 2
    while True:
        an = aa(n)
        bn = bb(n)

        D = D * an + bn
        if D == 0:
            D = dm

        C = bn + an/C
        if C == 0:
            C = dm

        D = 1/D
        delta = C*D
        f = f * delta
        n = n + 1
        if abs(delta - 1) < eps:
            return f

def log1p(x):
    '''compute log(1+x) accurately for cases where |x| is small.

    Implemented using the Maclaurin series for log(1 + x).

    Args:
        x (float): a floating point number > -1.

    Returns:
        an accurate evaluation of log(1 + x).
    '''

    if abs(x) >= 1:
        return math.log(1 + x)

    mx = -x
    mxn = 1
    s = 0
    n = 0
    while True:
        n = n + 1
        mxn *= mx
        s -= mxn/n
        if abs(mxn/s) < eps:
            return s

def expm1(x):
    '''compute exp(x) - 1 accurately for cases where x < 1.

    Args:
        x (float): a floating point value.

    Returns:
        an accurate evaluation of exp(x) - 1.
    '''
    if abs(x) >= 1:
        return math.exp(x) - 1

    s = 0
    n = 0
    xn = 1
    nfac = 1
    while True:
        n = n + 1
        xn *= x
        nfac = nfac * n
        t = xn / nfac
        s += t
        if abs(t/s) < eps:
            return s

def log1mexp(x):
    '''compute log(1 - exp(x) accurately.'''
    a = -x
    if a < log_2:
        return math.log(-expm1(x))
    else:
        return log1p(-math.exp(x))

def log1pexp(x):
    '''compute log(1 + exp(x) accurately.'''
    if x <= -37:
        return math.exp(x)
    elif x < 18:
        return log1p(math.exp(x))
    elif x < 33.3:
        return x + math.exp(-x)
    else:
        return x

def logAdd(a, b):
    '''compute log(exp(a) + exp(b))'''
    y, x = minMax(a, b)
    w = y - x
    return x + log1pexp(w)

def logGamma(z):
    '''compute log(gamma(x))'''
    x0 = 9

    if z < x0:
        n = int(math.floor(x0) - math.floor(z))
        p = 1.0
        for k in xrange(0, n):
            p *= z+k
        return logGamma(z + n) - math.log(p)
    else:
        z2 = z*z
        z3 = z2*z
        z5 = z3*z2
        return z*math.log(z) - z - 0.5*math.log(z/(2*math.pi)) + 1.0/(12*z) + 1.0/(360*z3) + 1.0/(1260*z5)

def gamma(z):
    '''compute gamma(x)'''
    return math.exp(logGamma(z))

def logFac(n):
    '''compute log(n!)'''
    return logGamma(n + 1)

def fac(n):
    '''compute n!'''
    r = 1
    for i in xrange(1,n+1):
        r *= i
    return r

def logChoose(n, k):
    '''compute log of binomial choice log(n k)'''
    return logFac(n) - (logFac(n - k) + logFac(k))

def choose(N, K):
    '''compute binomial choice (n k)'''
    (J,I) = minMax(N - K, K)
    num = 1
    n = I + 1
    while n <= N:
        num *= n
        n += 1
    den = 1
    j = 1
    while j <= J:
        den *= j
        j += 1
    return num // den

