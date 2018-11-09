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

def logSub(a, b):
    w = b - a
    return a + log1mexp(w)

def fac(n):
    '''compute n!'''
    r = 1
    for i in range(1,n+1):
        r *= i
    return r

def doubleFac(n):
    '''compute n!! = n*(n - 2)*(n - 4)...'''
    assert type(n) == int
    k = n
    r = 1
    while k > 0:
        r *= k
        k -= 2
    return r

def gammanph(n):
    '''compute gamma(n + 1/2)'''
    return sqrt_pi * doubleFac(2*n - 1) / math.pow(2, n)

def logGammaP1(a):
    '''
    log(gamma(1+a)) for -0.2 <= a <= 1.25
    Derived from AS63.
    '''
    p0 = 0.577215664901533e+00
    p1 = 0.844203922187225e+00
    p2 = -0.168860593646662e+00
    p3 = -0.780427615533591e+00
    p4 = -0.402055799310489e+00
    p5 = -0.673562214325671e-01
    p6 = -0.271935708322958e-02

    q1 = 0.288743195473681e+01
    q2  = 0.312755088914843e+01
    q3 = 0.156875193295039e+01
    q4  = 0.361951990101499e+00
    q5 = 0.325038868253937e-01
    q6  = 0.667465618796164e-03

    r0 = 0.422784335098467e+00
    r1  = 0.848044614534529e+00
    r2 = 0.565221050691933e+00
    r3  = 0.156513060486551e+00
    r4 = 0.170502484022650e-01
    r5  = 0.497958207639485e-03

    s1 = 0.124313399877507e+01
    s2  = 0.548042109832463e+00
    s3 = 0.101552187439830e+00
    s4  = 0.713309612391000e-02
    s5 = 0.116165475989616e-03
    if a < 0.6:
        w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0) / ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0)
        return -a * w

    x = (a - 0.5) - 0.5
    w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0) / (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0)
    return x * w

def logGamma(z):
    '''
    compute log(gamma(x))
    Implementation due to AS63
    '''
    x0 = 25

    if z < 0.8:
        return logGammaP1(z) - math.log(z)

    if z < 2.25:
      return logGammaP1(z - 1.0)

    if z < 10:
        return logGamma(z - 1.0) + math.log(z - 1.0)

    d = 0.418938533204673

    c0 = 0.833333333333333e-01
    c1 = -0.277777777760991e-02
    c2 = 0.793650666825390e-03
    c3 = -0.595202931351870e-03
    c4 = 0.837308034031215e-03
    c5 = -0.165322962780713e-02

    t = (1.0/z) ** 2
    w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0) / z
    return (d+w) + (z-0.5) * (math.log(z)-1.0)

def gamma(z):
    '''compute gamma(x)'''
    return math.exp(logGamma(z))

def logFac(n):
    '''compute log(n!)'''
    return logGamma(n + 1)

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

