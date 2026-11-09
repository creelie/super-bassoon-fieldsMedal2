"""Certified constants for the paper "Split Primes and the Exponent of Imaginary Quadratic Class Groups".

Ball arithmetic of Arb through python-flint (pip install python-flint).  Every assertion below is a
statement about intervals, so a passing run proves the inequality it checks.

Usage: python3 -I certify.py
"""
import math
from flint import arb, arb_series, ctx

ctx.prec = 256


def mzeta(s):
    """-zeta'/zeta at the real point s > 1."""
    c = arb_series([arb(s), 1], prec=2).zeta().coeffs()
    return -c[1] / c[0]


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]


# ---------------------------------------------------------------- Theorem 1: the least split prime
def split_margin(B, L, s0, lq0, c):
    """Proposition 4.1 at log q0 = lq0, with y chosen so that Y = e^L y = (c log q0)^k."""
    B, L, s0, lq0, c = arb(B), arb(L), arb(s0), arb(lq0), arb(c)
    k = 1 / (1 - B)
    logy = k * (c * lq0).log() - L
    y = logy.exp()
    Y = (L + logy).exp()
    assert y > 41 and y > L.exp()
    Psi = ((1 - 1 / logy) * (L.exp() - 1) - arb('1.03883') * (1 - (-L).exp())) * y / L
    R = lq0 + 2 * L * lq0 / arb(2).log()
    E = arb('1.03883') * (Y.sqrt() + Y.root(3) * Y.log() / arb(2).log())
    T = (L.exp() / y) * arb.pi() ** 2 / 8 / L
    h0 = -arb.pi().log() / 2 + ((s0 + 1) / 2).digamma() / 2 + mzeta(s0)
    assert h0 > 0
    H = lq0 / 2 + h0
    m1, m2 = 4 * (s0 - arb('0.5')), 1 / (s0 - B)
    assert m1 > m2 or m2 > m1
    M = m1 if m1 > m2 else m2
    assert s0 / (1 - B) ** 2 > 1 / (s0 - arb('0.5'))
    Z = 4 / L * H * ((B * logy).exp() * (B * L / 2).cosh() ** 2 * M
                     + y.sqrt() * (L / 4).cosh() ** 2 * s0 / (1 - B) ** 2)
    margin = Psi - R - 2 * E - T - Z
    assert margin > 0
    return dict(logy0=logy, h0=h0, M=M, Psi_y=Psi / y, R_y=R / y, E_y=E / y, T_y=T / y,
                Z_y=Z / y, margin_y=margin / y)


def small_range(c, k, lq0):
    """(|d|+4)/4 < (c log|d|)^k for 7 <= |d| <= e^lq0, checked on intervals [a, 1.05a]."""
    a = 7.0
    n = 0
    while a < math.exp(lq0):
        b = a * 1.05
        assert (arb(b) + 4) / 4 < (arb(c) * arb(a).log()) ** k, (a, b)
        a = b
        n += 1
    assert 7 < (arb(c) * arb(3).log()) ** k and 5 < (arb(c) * arb(4).log()) ** k
    return n


# ---------------------------------------------------------------- Theorem 2: exponent bounds
def exponent_bound(k, e):
    """Largest x with log(x/4) <= k e log(2 log x); returned as log10 x (upper end of the ball)."""
    f = lambda t: (t - arb(4).log()) - k * e * (2 * t).log()   # t = log x
    lo, hi = arb(k * e), arb(10000)                            # f increasing for t > k e
    assert f(hi) > 0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return hi / arb(10).log()


# ---------------------------------------------------------------- Theorem 3: L(1, chi) for q >= 2^40
B1 = arb('0.26149721284764278375542683860869585905156664826120')   # Mertens' constant
G6 = arb(0)
for p in primes_upto(10 ** 6):
    G6 += 1 / arb(p) - (1 + 1 / arb(p)).log()


def lbound(B, y0, ell, u, q0=arb(2) ** 40):
    B, y0, ell, u = arb(B), arb(y0), arb(ell), arb(u)
    k = 1 / (1 - B)
    x0 = ell.exp() * y0
    c0 = ((1 + u / 2).digamma() - arb.pi().log()) / 2 + mzeta(1 + u)
    assert c0 > 0 and y0 >= 10 ** 6
    S0 = (1 + u - B) / (1 - B) ** 2 * (((q0 / arb.pi()).log() + (1 + u / 2).digamma()) / 2 + mzeta(1 + u))
    E0 = ((x0 ** (B - 1) / x0.log() + y0 ** (B - 1) / y0.log()) * S0
          + arb.pi() ** 2 / 8 * (1 / (x0 * x0.log()) + 1 / (y0 * y0.log()))) / ell
    kappa0 = y0.log() + ell / 2 - k * q0.log().log()
    Cm = (G6 - B1 - 1 / (2 * y0.log() ** 2) - 3 / y0.sqrt() - E0).exp()
    hq = lambda q: Cm * q.sqrt() / (arb.pi() * (k * q.log().log() + kappa0))
    return dict(k=k, c0=c0, S0=S0, E0=E0, kappa0=kappa0, Cm=Cm, h40=hq(q0), hq=hq)


def exponent_two_tail(r):
    """For fundamental d with |d| > 2^40 and Cl(d) of exponent <= 2: contradiction.

    h(d) = 2^(omega-1).  If 2^(omega-1) < h40 the lower bound h(d) >= h40 is violated.  Otherwise
    |d| >= p_t# for t = omega, and h_lower(p_t#) > 2^(t-1) at the first such t; the ratio
    h_lower(qp)/h_lower(q) >= sqrt(p) A/(A + k log 2) >= 2 propagates it to all larger t."""
    k, kappa0, hq = r['k'], r['kappa0'], r['hq']
    t = 1
    while arb(2) ** (t - 1) < r['h40']:
        t += 1
    ps = primes_upto(1000)
    prim = arb(1)
    for p in ps[:t]:
        prim *= p
    assert prim > arb(2) ** 40
    q = prim
    assert hq(q) > arb(2) ** (t - 1)
    A = k * (arb(2) ** 40).log().log() + kappa0
    assert arb(ps[t]).sqrt() * A / (A + k * arb(2).log()) >= 2
    return t, q, hq(q)


if __name__ == '__main__':
    print('== Theorem 1 ==')
    for B, L, s0 in [('0.875', 8, '1.225'), (arb(11) / 12, 8, '1.27')]:
        r = split_margin(B, L, s0, 22, 2)
        print('B =', B, {kk: v.str(8) for kk, v in r.items()})
    print('small range k=8 intervals:', small_range(2, 8, 22))
    print('small range k=12 intervals:', small_range(2, 12, 22))
    print('== Theorem 2: log10 of the largest |d| with E(d) <= e ==')
    for e in range(2, 11):
        print(e, exponent_bound(8, e).str(5), exponent_bound(12, e).str(5))
    print('== Theorem 3 ==')
    print('G(10^6) =', G6.str(15))
    for B, y0, ell, u in [(arb(7) / 8, 5 * 10 ** 6, 14, arb(1) / 10),
                          (arb(11) / 12, arb('1.8e9'), 20, arb('0.08'))]:
        r = lbound(B, y0, ell, u)
        print('B =', B.str(6), {kk: r[kk].str(10) for kk in ['c0', 'S0', 'E0', 'kappa0', 'Cm', 'h40']})
        t, q, hq = exponent_two_tail(r)
        print('   exponent two: omega >=', t, 'needs |d| >=', q.str(5), 'where h >=', hq.str(6),
              '>', 2 ** (t - 1))
    print('all assertions passed')
