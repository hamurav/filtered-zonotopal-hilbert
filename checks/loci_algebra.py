#!/usr/bin/env python3
"""Symbolic identities of Section 6.10 (Remarks 6.21 and 6.22); needs SymPy.

1. (6.19): v_f = f'(g) = 1/g' has v_1 = 2 alpha, v_2 = 3 beta - 2 alpha^2,
   v_3 = 4(alpha^3 - 2 alpha beta + gamma).
2. Remark 6.21: the loci A, B, C, Q, L, W, E, K in the chart
   f = u + u^2 + lambda u^3 + mu u^4; [u^3]g = 0 on A and [u^4]g = 0 on B;
   the incidences (C, Q, L, W through E; A, B, Q, L, W through K); A, L, Q
   are the 2 x 2 minors of (1 v1 v2; v1 v2 v3), and the same minors for the
   coefficients of 1/v_f = g' give C, L, Q.
3. Remark 6.22: on Q (v1 != 0), v_f = (1 + sigma u)/(1 + tau u) mod u^4 and
   the formulas for g mod u^5; the (a, b) form of W, its discriminant
   4(a-1)^2(5a^2-4a+1) in b, and the affine chart change (W in the chart is
   -64 times the (a, b) form).

Usage:  python3 loci_algebra.py
"""

import sys

import sympy as sp

u, al, be, ga, de, lam, mu = sp.symbols('u alpha beta gamma delta lambda mu')
v1, v2, v3, a, b, s, t = sp.symbols('v1 v2 v3 a b sigma tau')
OK = []


def check(cond, msg):
    OK.append(bool(cond))
    print(f'{"ok  " if cond else "FAIL"} {msg}')


def zero(e):
    return sp.simplify(sp.expand(e)) == 0


def series_inverse(f, order):
    """Compositional inverse g of f = u + ..., modulo u^order."""
    g = u
    for _ in range(order):
        g = sp.expand(g - (f.subs(u, g) - u))
        g = sum(g.coeff(u, k) * u ** k for k in range(1, order))
    return g


def coeffs(expr, n):
    e = sp.expand(sp.series(expr, u, 0, n + 1).removeO())
    return [e.coeff(u, k) for k in range(n + 1)]


def main():
    f = u + al * u ** 2 + be * u ** 3 + ga * u ** 4 + de * u ** 5
    g = series_inverse(f, 6)
    vf = coeffs(1 / sp.diff(g, u), 3)
    check(zero(vf[1] - 2 * al) and zero(vf[2] - (3 * be - 2 * al ** 2)) and
          zero(vf[3] - 4 * (al ** 3 - 2 * al * be + ga)),
          '(6.19): v1 = 2 alpha, v2 = 3 beta - 2 alpha^2, '
          'v3 = 4(alpha^3 - 2 alpha beta + gamma)')
    vf2 = coeffs(sp.diff(f, u).subs(u, g), 3)
    check(all(zero(x - y) for x, y in zip(vf, vf2)), "v_f = f'(g) = 1/g'")

    # the chart alpha = 1, beta = lambda, gamma = mu
    ch = {al: 1, be: lam, ga: mu}
    V = {v1: vf[1].subs(ch), v2: vf[2].subs(ch), v3: vf[3].subs(ch)}
    gc = [sp.expand(g.coeff(u, k).subs(ch)) for k in range(6)]
    W = (v1**2 * v3**2 + 2 * v1 * v2**2 * v3 - 4 * v2**4 + 8 * v1**2 * v2**3
         - 6 * v1**3 * v2 * v3 - 3 * v1**4 * v2**2 + 2 * v1**5 * v3)
    loci = {
        'A': (v2 - v1**2, lam - 2),
        'B': (v3 - (2 * v1 * v2 - v1**3), mu - (5 * lam - 5)),
        'C': (v2, lam - sp.Rational(2, 3)),
        'Q': (v1 * v3 - v2**2, 8 * mu - (9 * lam**2 + 4 * lam - 4)),
        'L': (v3 - v1 * v2, 2 * mu - (7 * lam - 4)),
        'W': (W, 81 * lam**4 - 360 * lam**3 - 36 * lam**2 * mu + 272 * lam**2
              + 256 * lam * mu + 80 * lam - 16 * mu**2 - 208 * mu - 64),
    }
    for name, (veq, cheq) in loci.items():
        e = sp.expand(veq.subs(V))
        q = sp.cancel(e / cheq)
        check(q.free_symbols == set() and q != 0,
              f'{name}: chart equation {sp.sstr(cheq)} = 0 '
              f'(v-equation / chart equation = {q})')
    check(sp.cancel(gc[3] / (lam - 2)).free_symbols == set(),
          '[u^3]g vanishes exactly on A (lambda = 2)')
    check(sp.cancel(gc[4] / (mu - 5 * lam + 5)).free_symbols == set(),
          '[u^4]g vanishes exactly on B (mu = 5 lambda - 5)')

    E = {lam: sp.Rational(2, 3), mu: sp.Rational(1, 3)}
    K = {lam: 2, mu: 5}
    fE = sp.expand(sp.series((sp.exp(2 * u) - 1) / 2, u, 0, 5).removeO())
    check(fE == (u + u**2 + lam * u**3 + mu * u**4).subs(E),
          'E = (2/3, 1/3) is the exponential (e^(2u) - 1)/2')
    gK = series_inverse(u + u**2 + 2 * u**3 + 5 * u**4, 5)
    check(sp.expand(gK - (u - u**2)) == 0,
          'K = (2, 5): f^(-1) = u - u^2 mod u^5 (inverse quadratic)')
    vE = {k: x.subs(E) for k, x in V.items()}
    vK = {k: x.subs(K) for k, x in V.items()}
    check(vE[v2] == 0 and vE[v3] == 0, 'E: v2 = v3 = 0')
    check(vK[v2] == vK[v1]**2 and vK[v3] == vK[v1]**3,
          'K: v2 = v1^2, v3 = v1^3')
    for name in 'CQLW':
        check(loci[name][1].subs(E) == 0, f'{name} passes through E')
    for name in 'ABQLW':
        check(loci[name][1].subs(K) == 0, f'{name} passes through K')
    check(loci['A'][1].subs(E) != 0 and loci['B'][1].subs(E) != 0 and
          loci['C'][1].subs(K) != 0,
          'A, B do not pass through E; C does not pass through K')

    Mv = sp.Matrix([[1, v1, v2], [v1, v2, v3]])
    minors = [Mv[:, [0, 1]].det(), Mv[:, [1, 2]].det(), Mv[:, [0, 2]].det()]
    check(zero(minors[0] - loci['A'][0]) and zero(minors[1] - loci['Q'][0])
          and zero(minors[2] - loci['L'][0]),
          'A, Q, L are the 2 x 2 minors of (1 v1 v2; v1 v2 v3)')
    w = coeffs(1 / (1 + v1 * u + v2 * u**2 + v3 * u**3), 3)
    Mw = sp.Matrix([[1, w[1], w[2]], [w[1], w[2], w[3]]])
    mw = [sp.expand(Mw[:, [0, 1]].det()), sp.expand(Mw[:, [1, 2]].det()),
          sp.expand(Mw[:, [0, 2]].det())]
    check(zero(mw[0] + v2) and zero(mw[1] - (v1 * v3 - v2**2)) and
          zero(mw[2] + (v3 - v1 * v2)),
          'the same minors for g\' = 1/v_f are -v2 (C), Q and -(L)')

    # Remark 6.22
    sig, tau = v1 - v2 / v1, -v2 / v1
    fr = coeffs((1 + s * u) / (1 + t * u), 3)
    check(zero(fr[1] - (s - t)) and zero(fr[2] + t * (s - t)) and
          zero(fr[3] - t**2 * (s - t)),
          '(1+sigma u)/(1+tau u) = 1 + (sigma-tau)u - tau(sigma-tau)u^2 '
          '+ tau^2(sigma-tau)u^3 + ...')
    sub = [x.subs({s: sig, t: tau}) for x in fr]
    check(zero(sub[1] - v1) and zero(sub[2] - v2) and
          zero((sub[3] - v3) * v1 - (v2**2 - v1 * v3)),
          'with sigma = v1 - v2/v1, tau = -v2/v1 this is v_f mod u^4 '
          'exactly on Q')
    glog = t / s * u + (s - t) / s**2 * sp.log(1 + s * u)
    target = sp.integrate(sp.series((1 + t * u) / (1 + s * u), u, 0, 4)
                          .removeO(), u)
    check(all(zero(x - y) for x, y in zip(coeffs(glog, 4), coeffs(target, 4))),
          "g = (tau/sigma) u + ((sigma-tau)/sigma^2) log(1+sigma u) mod u^5 "
          "integrates g' = (1+tau u)/(1+sigma u)")
    target0 = sp.integrate((1 + t * u), u)
    check(zero(target0 - (u + t / 2 * u**2)), 'sigma = 0: g = u + tau/2 u^2')

    Wab = b**2 + (2 * a**2 - 6 * a + 2) * b - 4 * a**4 + 8 * a**3 - 3 * a**2
    check(zero(W.subs({v2: a * v1**2, v3: b * v1**3}) - v1**8 * Wab),
          'W / v1^8 = b^2 + (2a^2 - 6a + 2)b - 4a^4 + 8a^3 - 3a^2')
    disc = sp.discriminant(Wab, b)
    check(zero(disc - 4 * (a - 1)**2 * (5 * a**2 - 4 * a + 1)),
          f'discriminant in b = {sp.factor(disc)}')
    check(sp.discriminant(5 * a**2 - 4 * a + 1, a) != 0,
          '5a^2 - 4a + 1 has distinct roots (discriminant -4)')
    A_, B_ = (3 * lam - 2) / 4, (1 - 2 * lam + mu) / 2
    check(zero(V[v2] / V[v1]**2 - A_) and zero(V[v3] / V[v1]**3 - B_),
          'chart: a = (3 lambda - 2)/4, b = (1 - 2 lambda + mu)/2')
    check(zero(loci['W'][1] + 64 * Wab.subs({a: A_, b: B_})),
          'W in the chart = -64 (a, b) form')

    ok = all(OK)
    print('all checks passed' if ok else 'SOME CHECK FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
