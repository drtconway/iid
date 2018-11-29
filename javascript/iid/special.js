var basic = require('./basic');

function pickGammaMethod(a, x) {
    var alpha;
    if (x > 0.25) {
        alpha = x + 0.25;
    } else {
        alpha = Math.log(0.5) / Math.log(x);
    }

    if (a > alpha) {
        return 'g';
    }

    var x0 = 1.5;
    if (x < x0 && -0.5 <= a) {
        return 'G1';
    }
    if (x < x0 && a < -0.5) {
        assert(false, 'incomplete gamma for a < -0.5 not implemented');
        //return 'G2'
    }
    return 'G3';
}

function upperGammaKummer(a, x) {
    var u = (basic.gamma(1 + a) - 1) / a - (Math.pow(x, a) - 1) / a;

    var k = 1;
    var tk = 1;
    var ss = tk;
    while (true) {
        tk *= - (a + k)*x / ((a + k + 1)*(k + 1));
        ss += tk;
        if (abs(tk/ss) < 1e-16) {
            break
        }
        k += 1;
    }
    var v = (Math.pow(x, a + 1) / (a + 1)) * ss;
    return u + v;
}

function logUpperGammaKummer(a, x) {
    var la = Math.log(a);
    var lx = Math.log(x);

    var lga = basic.logGamma(a);
    var lxaoa = a * lx - la;
    var flip;
    var lu;
    if (lga > lxaoa) {
        flip = false;
        lu = basic.logSub(lga, lxaoa);
    } else {
        flip = true;
        lu = basic.logSub(lxaoa, lga);
    }

    var ltk = 0;
    var tkPos = true;
    var lss = ltk;
    var k = 1;
    while (true) {
        tkPos = !tkPos;
        ltk += Math.log(a + k) + lx - (Math.log(a + k + 1) + Math.log(k + 1));
        if (tkPos) {
            lss = basic.logAdd(lss, ltk);
        } else {
            lss = basic.logSub(lss, ltk);
        }
        if (ltk - lss < -45) {
            break;
        }
        k += 1;
    }

    var lv = (a + 1)*lx - Math.log(a + 1) + lss;

    if (!flip) {
        return basic.logAdd(lu, lv);
    } else {
        return basic.logSub(lv, lu);
    }
}

function upperGammaSeries(a, x) {
    var s = 0;
    var k = 1;
    var m1k = 1;
    var xk = 1;
    var kfac = 1;
    var ts = [];
    while (true) {
        m1k *= -1;
        xk *= x;
        kfac *= k;
        t = m1k*xk/((a+k)*kfac);
        s += t;
        ts.append(t);
        if (abs(t/s) < 1e-14) {
            break;
        }
        k += 1;
    }

    var gap1m1 = basic.gamma(a + 1) - 1;
    return (gap1m1 - basic.powm1(x, a)) / a + Math.pow(x, a)*s;
}

function upperGammaLegendre(a, x) {
    var pk = 0;
    var qk = (x - 1 - a)*(x + 1 - a);
    var rk = 4*(x + 1 - a);
    var sk = 1 - a;
    var rhok = 0;
    var tk = 1;
    var k = 1;
    var ss = tk;
    while (true) {
        pk += sk;
        qk += rk;
        rk += 8;
        sk += 2;
        tauk = pk*(1 + rhok);
        rhok = tauk / (qk - tauk);
        tk *= rhok;
        ss += tk;
        if (abs(tk/ss) < 1e-16) {
            break;
        }
        k += 1;
    }
    return Math.exp(-x)*Math.pow(x, a)*ss / (x + 1 - a);
}

function logUpperGammaLegendre(a, x) {
    var pk = 0;
    var qk = (x - 1 - a)*(x + 1 - a);
    var rk = 4*(x + 1 - a);
    var sk = 1 - a;
    var rhok = 0;
    var tk = 1;
    var k = 1;
    var ss = tk;
    while (true) {
        pk += sk;
        qk += rk;
        rk += 8;
        sk += 2;
        tauk = pk*(1 + rhok);
        rhok = tauk / (qk - tauk);
        tk *= rhok;
        ss += tk;
        if (abs(tk/ss) < 1e-16) {
            break;
        }
        k += 1;
    }
    return a * Math.log(x) - x + Math.log(ss) - basic.log1p(x - a);
}

function lowerGammaSeries(a, x) {
    var t = 1.0;
    var s = 0.0;
    while (true) {
        a += 1;
        t *= x/a;
        s += t;
        if (abs(t/s) < 1e-14) {
            break;
        }
    }
    return s;
}

function lowerRegularizedGammaSeries(a, x) {
    var ss =  1.0 / basic.gamma(a + 1);
    var xn = 1
    var n = 1;
    while (true) {
        xn *= x;
        var tn = xn / basic.gamma(a + n + 1);
        ss += tn;
        if (abs(tn/ss) < 1e-16) {
            break;
        }
        n += 1;
    }
    return Math.exp(-x)*Math.pow(x, a)*ss;
}

function logLowerRegularizedGammaSeries(a, x) {
    var lx = Math.log(x);
    var lss =  - basic.logGamma(a + 1);
    var lxn = 0;
    var n = 1;
    while (true) {
        lxn += lx;
        var ltn = lxn - basic.logGamma(a + n + 1);
        lss = basic.logAdd(lss, ltn);
        if (ltn - lss < -45) {
            break;
        }
        n += 1;
    }
    return a*lx - x + lss;
}

function gammaP(a, x) {
    var alg = pickGammaMethod(a, x);
    if (alg == 'g') {
        return lowerRegularizedGammaSeries(a, x);
    }

    var Gam;
    if (alg == 'G1') {
        Gam = upperGammaKummer(a, x);
    }
    if (alg == 'G3') {
        Gam = upperGammaLegendre(a, x);
    }

    var ga = basic.gamma(a);
    var bigG = Gam / ga;
    var gs = 1 - bigG;
    return gs;
}

function logGammaP(a, x) {
    var alg = pickGammaMethod(a, x);
    if (alg == 'g') {
        return logLowerRegularizedGammaSeries(a, x);
    }

    var lGam;
    if (alg == 'G1') {
        lGam = logUpperGammaKummer(a, x);
    }
    if (alg == 'G3') {
        lGam = logUpperGammaLegendre(a, x);
    }

    var lga = basic.logGamma(a);
    var lBigG = lGam - lga;
    var lgs = basic.log1mexp(lBigG);
    return lgs;
}

function gammaQ(a, x) {
    var alg = pickGammaMethod(a, x);

    if (alg == 'g') {
        var gs = lowerRegularizedGammaSeries(a, x);
        return 1 - gs;
    }

    var Gam;
    if (alg == 'G1') {
        Gam = upperGammaKummer(a, x)
    }
    if (alg == 'G3') {
        Gam = upperGammaLegendre(a, x)
    }

    return Gam / basic.gamma(a);
}

function logGammaQ(a, x) {
    var alg = pickGammaMethod(a, x);

    if (alg == 'g') {
        var lgs = logLowerRegularizedGammaSeries(a, x);
        return basic.log1mexp(lgs);
    }

    var lGam;
    if (alg == 'G1') {
        lGam = logUpperGammaKummer(a, x)
    }
    if (alg == 'G3') {
        lGam = logUpperGammaLegendre(a, x)
    }

    return lGam - basic.logGamma(a);
}

function gammaPQ(a, x) {
    var alg = pickGammaMethod(a, x);

    if (alg == 'g') {
        var gs = lowerRegularizedGammaSeries(a, x);
        return [gs, 1 - gs];
    }

    var Gam;
    if (alg == 'G1') {
        Gam = upperGammaKummer(a, x)
    }
    if (alg == 'G3') {
        Gam = upperGammaLegendre(a, x)
    }

    var bigG = Gam / basic.gamma(a);
    return [1 - bigG, bigG];
}

function logGammaPQ(a, x) {
    var alg = pickGammaMethod(a, x);

    if (alg == 'g') {
        var lgs = logLowerRegularizedGammaSeries(a, x);
        return [lgs, basic.log1mexp(lgs)];

    var lGam;
    if (alg == 'G1') {
        lGam = logUpperGammaKummer(a, x);
    }
    if (alg == 'G3') {
        lGam = logUpperGammaLegendre(a, x);
    }

    var lBigG = lGam - basic.logGamma(a);
    return [basic.log1mexp(lBigG), lBigG];
}

function logBeta(a, b) {
    return basic.logGamma(a) + basic.logGamma(b) - basic.logGamma(a + b);
}

function beta(a, b) {
    return Math.exp(logBeta(a, b));
}

function betaIntImpl(a, b) {
    if (a < b) {
        return betaIntImpl(b, a);
    }

    var r = 1;
    for (var i = 1; i < b; i++) {
        r *= i / (a + i);
    }
    return r;
}

function logBetaInt(a, b) {
    var lga = basic.logGamma(a);
    var lgb = basic.logGamma(b);
    var lgab = basic.logGamma(a+b);
    var lb = lga + lgb - lgab;
    if (lb > -700) {
        return Math.log(betaIntImpl(a, b));
    }
    return lb;
}

function betaInt(a, b) {
    var lga = basic.logGamma(a);
    var lgb = basic.logGamma(b);
    var lgab = basic.logGamma(a+b);
    var lb = lga + lgb - lgab;
    if (lb > -700) {
        return betaIntImpl(a, b);
    }
    return Math.exp(lb);
}

function lowerBetaInt(a, b, x) {
    var y = 1 - x;
    if (x > a/(a+b)) {
        return 1 - lowerBetaInt(b, a, y);
    }

    var lx = Math.log(x);
    var ly = basic.log1p(-x);
    var lpfx = a*lx + b*ly;

    if (a < 350 && b < 350 && lpfx > -700) {
        var xa = Math.pow(x, a);
        var yb = Math.pow(y, b);
        var bx = xa*yb/a * basic.hyper([a+b, 1], [a+1], x);
        var bab = betaInt(a, b);
        return bx/bab;
    } else {
        var lbx = lpfx - Math.log(a) + basic.logHyper([a+b, 1], [a+1], x);
        var lbab = logBetaInt(a, b);
        if (lbx - lbab > 0) {
            lbab = lbx;
        }
        return Math.exp(lbx - lbab);
    }
}

function erf(x) {
    if (x < 0) {
        return - erf(-x);
    }

    if (x > 5) {
        return 1 - Math.exp(-x*x)/basic.sqrt_pi;
    }

    if (x > 4) {
        return 1 - Math.exp(-x*x)/(basic.sqrt_pi*x);
    }

    var s = x;
    var n = 1;
    var p = 1;
    while (true) {
       p *= -x*x / n;
       var t = x / (2*n + 1) * p;
       s += t;
       if (s == 0 || abs(t/s) < 1e-12) {
           break;
        }
       n += 1;
    }
    return 2.0/basic.sqrt_pi * s;
}

function erfc(x) {
    if (x < - 0.5) {
        return 2.0 - erfc(-x);
    }
    if (x < 0.5) {
        return 1 - erf(x);
    }
    if (x < 50) {
        var aa = function(j) { return 0.5*j };
        var bb = function(j) : { return x };
        var v = basic.contFrac(aa, bb);
        return Math.exp(-x*x)/(basic.sqrt_pi*(x + v));
    }
    return gammaQ(0.5, x*x)/basic.sqrt_pi;
}

function logErfc(x) {
    if (x < - 0.5 ) {
        return Math.log(2.0 - erfc(-x));
    }
    if (x < 0.5) {
        return basic.log1p(-erf(x));
    }
    if (x < 50) {
        var aa = function(j) { return 0.5*j };
        var bb = function(j) : { return x };
        var v = basic.contFrac(aa, bb);
        return -x*x - basic.log_sqrt_pi - Math.log(x + v);
    }
    return logGammaQ(0.5, x*x) - basic.log_sqrt_pi;
}
