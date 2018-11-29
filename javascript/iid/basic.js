var pi = 3.14159265358979323846264338327950288

var log_pi = Math.log(pi)

var log_2 = Math.log(2)

var sqrt_pi = Math.sqrt(pi)

var log_sqrt_pi = Math.log(sqrt_pi)

function abs(x) {
    if (x >= 0) {
        return x;
    } else {
        return -x;
    }
}

function sgn(x) {
    if (x > 0) {
        return 1;
    }
    if (x < 0) {
        return -1;
    }
    return 0;
}

function min(a, b) {
    if (a <= b) {
        return a;
    } else {
        return b;
    }
}

function max(a, b) {
    if (a >= b) {
        return a;
    } else {
        return b;
    }
}

function minMax(a, b) {
    if (a <= b) {
        return [a, b];
    } else {
        return [b, a];
    }
}

function contFrac(aa, bb) {
    var dm = 1e-300;
    var a1 = aa(1);
    var b1 = bb(1);
    var f = a1/b1;
    var C = a1/dm;
    var D = 1/b1;
    var n = 2;
    while (true) {
        var an = aa(n);
        var bn = bb(n);
        D = D * an + bn;
        if (D == 0) {
            D = dm;
        }
        C = bn + an/C;
        if (C == 0) {
            C = dm;
        }
        var delta = C*D;
        f *= delta;
        n += 1;
        if (abs(delta - 1) < 1e-14) {
            return f;
        }
    }
}

function log1p(x) {
    if (abs(x) > 0.8) {
        return Math.log(1 + x);
    }
    var mx = -x;
    var mxn = 1;
    var s = 0;
    var n = 0;
    while (true) {
        n += 1;
        mxn *= mx;
        s -= mxn / n;
        if (abs(mxn/s) < 1e-16) {
            return s;
        }
    }
}

function expm1(x) {
    if (abs(x) >= 1) {
        return Math.exp(x) - 1;
    }

    var s = 0;
    var n = 0;
    var xn = 1;
    var nfac = 1;
    while (true) {
        n += 1;
        xn *= x;
        nfac *= n;
        var t = xn / nfac;
        s += t;
        if (abs(t/s) < 1e-16) {
            return s;
        }
    }
}

function log1mexp(x) {
    var a = -x;
    if (a < log_2) {
        return Math.log(-expm1(x));
    } else {
        return log1p(-Math.exp(x));
    }
}

function powm1(x, y) {
    var l = y * Math.log(x);
    if (l < 0.5) {
        return expm1(l);
    }
    return Math.pow(x, y) - 1;
}

function logAdd(a, b) {
    var t = minMax(a, b);
    var y = t[0];
    var x = t[1];
    var w = y - x;
    return x + log1pexp(w);
}

function logSub(a, b) {
    w = b - a;
    return a + log1mexp(w);
}

function fac(n) {
    var r = 1;
    for (var i = 1; i <= n; i++) {
        r *= i;
    }
    return r;
}

function doubleFac(n) {
    var k = n;
    var r = 1;
    while (k > 0) {
        r *= k;
        k -= 2;
    }
    return r;
}

function kahanSum(xs) {
    var r = 0;
    var c = 0;
    for (var i = 0; i < xs.length; i++) {
        var x = xs[i];
        var y = x - c;
        var t = r + y;
        c = (t - r) - y;
        r = t;
    }
    return r;
}

function gammanph(n) {
    return sqrt_pi * doubleFac(2*n - 1) / Math.pow(2, n);
}

function logGammaP1(a) {
    var p0 = 0.577215664901533e+00;
    var p1 = 0.844203922187225e+00;
    var p2 = -0.168860593646662e+00;
    var p3 = -0.780427615533591e+00;
    var p4 = -0.402055799310489e+00;
    var p5 = -0.673562214325671e-01;
    var p6 = -0.271935708322958e-02;

    var q1 = 0.288743195473681e+01;
    var q2  = 0.312755088914843e+01;
    var q3 = 0.156875193295039e+01;
    var q4  = 0.361951990101499e+00;
    var q5 = 0.325038868253937e-01;
    var q6  = 0.667465618796164e-03;

    var r0 = 0.422784335098467e+00;
    var r1  = 0.848044614534529e+00;
    var r2 = 0.565221050691933e+00;
    var r3  = 0.156513060486551e+00;
    var r4 = 0.170502484022650e-01;
    var r5  = 0.497958207639485e-03;

    var s1 = 0.124313399877507e+01;
    var s2  = 0.548042109832463e+00;
    var s3 = 0.101552187439830e+00;
    var s4  = 0.713309612391000e-02;
    var s5 = 0.116165475989616e-03;
    if (a < 0.6) {
        var w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0) / ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0);
        return -a * w;
    }

    x = a - 1.0;
    var w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0) / (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0);
    return x * w;
}

function logGamma(z) {
    if (z < 0.8) {
        return logGammaP1(z) - Math.log(z);
    }

    if (z < 2.25) {
      return logGammaP1(z - 1.0);
    }

    if (z < 10) {
        return logGamma(z - 1.0) + Math.log(z - 1.0);
    }

    var d = 0.418938533204673;

    var c0 = 0.833333333333333e-01;
    var c1 = -0.277777777760991e-02;
    var c2 = 0.793650666825390e-03;
    var c3 = -0.595202931351870e-03;
    var c4 = 0.837308034031215e-03;
    var c5 = -0.165322962780713e-02;

    var t = (1.0/z) ** 2;
    var w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0) / z;
    return (d+w) + (z-0.5) * (Math.log(z)-1.0);
}

function logGammaSum(a, b) {
    var x = a + b - 2.0;
    if (x <= 0.25) {
        return logGammaP1(1.0 + x);
    }
    if (x <= 1.25) {
        return logGammaP1(x) + log1p(x);
    }
    return logGammaP1(x - 1.0) + Math.log(x*(1.0 + x));
}

function gamma(z) {
    return Math.exp(logGamma(z));
}

function logFac(n) {
    return logGamma(n + 1);
}

function logChoose(n, k) {
    if (k == 0 || k == n) {
        return 0.0;
    }
    if (k == 1 || n - k == 1) {
        return Math.log(n);
    }
    return logFac(n) - (logFac(n - k) + logFac(k));
}

function choose(N, K) {
    var t = minMax(N - K, K);
    var J = t[0];
    var I = t[1];
    var r = 1;
    var n = I + 1;
    var j = 1;
    while (n <= N && j <= J) {
        r *= n/j;
        n += 1;
        j += 1;
    }
    while (n <= N) {
        r *= n;
        n += 1;
    }
    while (j <= J) {
        r /= j;
        j += 1;
    }
    return r;
}

function hyper(num, den, x) {
    var N = max(num.length, den.length);
    var xN = 1;
    var n = 1;
    var s = 1;
    var u = 1.0;
    while (true) {
        for (var i = 0; i < N; i++) {
            if (i < num.length) {
                ai = num[i];
                v = ai + n - 1;
                if (v == 0) {
                    return s;
                }
                u *= v;
            }
            if (i < den.length) {
                bi = den[i];
                v = (bi + n - 1);
                u /= v;
            }
        }
        u /= n;
        xN *= x
        var t = u * xN
        if (!isFinite(t)) {
            return s;
        }
        s += t;
        if (abs(t/s) < 1e-14) {
            return s;
        }
        n += 1;
    }
}

function logHyper(num, den, x) {
    var lx = Math.log(x);
    var sgnNumPoc = [];
    var lnumPoc = [];
    for (var i = 0; i < num.length; i++) {
        sgnNumPoc.push(1);
        lnumPoc.push(0);
    }
    var sgnDenPoc = [];
    var ldenPoc = [];
    for (var i = 0; i < den.length; i++) {
        sgnDenPoc.push(1);
        ldenPoc.push(0);
    }
    var n = 1;
    var ls = 0;
    while (true) {
        var sgnT = 1;
        var lt = 0;
        for (var i = 0; i < num.length; i++) {
            var ai = num[i];
            var w = ai + n - 1;
            if (w == 0) {
                return ls;
            }
            sgnNumPoc[i] *= sgn(w);
            lnumPoc[i] += Math.log(abs(w));
            sgnT *= sgnNumPoc[i];
            lt += lnumPoc[i];
        }
        for (var i = 0; i < den.length; i++) {
            var bi = den[i];
            var w = bi + n - 1;
            sgnDenPoc[i] *= sgn(w);
            ldenPoc[i] += Math.log(abs(w));
            sgnT *= sgnDenPoc[i];
            lt -= ldenPoc[i];
        };
        lt += n * lx - logFac(n)
        if (sgnT > 0) {
            ls = logAdd(ls, lt);
        } else {
            ls = logSub(ls, lt);
        }
        if (lt - ls < -45) {
            return ls;
        }
        n += 1;
    }
}
