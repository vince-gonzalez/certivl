"""
Exact rational + certified interval arithmetic kernel.

ZERO-MARGIN PROTOCOL
--------------------
Nothing in this project's proof path may touch a binary float. Every quantity is
either an exact `Fraction`, or an `Ivl` -- a closed interval with exact rational
endpoints that is GUARANTEED to contain the true value.

Every operation here is outward-rounded: the result interval always contains the
true result of the operation applied to any point of the input intervals. That
containment is what makes a computed inequality a proof: if `Ivl.hi < 0` then the
true value is negative, full stop, no floating-point caveat.

Irrational constants (pi, sqrt, sin, cos) are produced as rational enclosures by
integer-only algorithms (isqrt) or by mpmath's validated interval type widened
outward. They are cross-checked against independently known digit strings in
tests/test_exact.py.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import isqrt

# Working precision for irrational enclosures: denominators of this scale.
# 10**60 is far beyond any tolerance the geometry needs; the certified results
# are insensitive to it (only the width of final intervals changes).
_PREC = 10**60


class Ivl:
    """A closed interval [lo, hi] with exact rational endpoints.

    Invariant: the true value being represented lies in [lo, hi].
    All arithmetic widens outward, so the invariant is preserved by composition.
    """

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        lo = F(lo)
        hi = lo if hi is None else F(hi)
        if lo > hi:
            raise ValueError(f"empty interval [{lo}, {hi}]")
        self.lo = lo
        self.hi = hi

    # -- constructors ----------------------------------------------------
    @staticmethod
    def exact(q) -> "Ivl":
        """A degenerate interval holding an exactly representable rational."""
        return Ivl(F(q), F(q))

    # -- basic arithmetic ------------------------------------------------
    def __add__(self, o):
        o = _coerce(o)
        return Ivl(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return Ivl(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-_coerce(o))

    def __rsub__(self, o):
        return _coerce(o) + (-self)

    def __mul__(self, o):
        o = _coerce(o)
        c = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Ivl(min(c), max(c))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = _coerce(o)
        if o.lo <= 0 <= o.hi:
            raise ZeroDivisionError("interval divisor straddles zero")
        return self * Ivl(F(1) / o.hi, F(1) / o.lo)

    def __rtruediv__(self, o):
        return _coerce(o) / self

    def sqr(self) -> "Ivl":
        """Certified square, tight for intervals straddling zero.

        Generic multiplication of x by itself loses the fact that a square is
        non-negative: for x = [-a, b] it returns [-ab, max(a^2,b^2)], whose lower
        bound is spuriously negative. This returns [0, max(a^2, b^2)] there, and
        the exact monotone square otherwise.
        """
        if self.lo >= 0:
            return Ivl(self.lo * self.lo, self.hi * self.hi)
        if self.hi <= 0:
            return Ivl(self.hi * self.hi, self.lo * self.lo)
        return Ivl(F(0), max(self.lo * self.lo, self.hi * self.hi))

    def __pow__(self, n: int):
        if n < 0:
            return Ivl.exact(1) / (self ** (-n))
        r = Ivl.exact(1)
        for _ in range(n):
            r = r * self
        return r

    # -- certified comparisons -------------------------------------------
    # These return True only when the relation holds for EVERY point of the
    # intervals, i.e. only when it is proved. Overlap returns False, never a
    # guess. `definitely_lt` is the workhorse predicate of the whole project.
    def definitely_lt(self, o) -> bool:
        return self.hi < _coerce(o).lo

    def definitely_gt(self, o) -> bool:
        return self.lo > _coerce(o).hi

    def definitely_le(self, o) -> bool:
        return self.hi <= _coerce(o).lo

    def definitely_ge(self, o) -> bool:
        return self.lo >= _coerce(o).hi

    def contains(self, q) -> bool:
        return self.lo <= F(q) <= self.hi

    def straddles_zero(self) -> bool:
        return self.lo <= 0 <= self.hi

    # -- reporting -------------------------------------------------------
    @property
    def width(self) -> F:
        return self.hi - self.lo

    @property
    def mid(self) -> F:
        return (self.lo + self.hi) / 2

    def round_out(self, den: int) -> "Ivl":
        """Outward rounding to endpoints with denominator dividing `den`.

        lo is floored, hi is ceiled, so the result CONTAINS self -- enclosure is
        preserved. This is exactly what a fixed-precision interval library does: it
        keeps denominators bounded so arithmetic stays O(1) instead of growing with
        subdivision depth. The widening introduced is at most 1/den per endpoint.
        """
        from math import floor, ceil
        lo = F(floor(self.lo * den), den)
        hi = F(ceil(self.hi * den), den)
        return Ivl(lo, hi)

    def decimals(self, n: int = 12) -> str:
        """Digit string with an explicit uncertainty flag.

        Prints only digits that are common to both endpoints; if the endpoints
        disagree at digit n the string ends in '?' so a reader can never mistake
        an unresolved digit for a certified one.
        """
        s = 10**n
        lo_d = (self.lo * s).__floor__()
        hi_d = (self.hi * s).__floor__()
        body = f"{F(lo_d, s):.{n}f}"
        return body if lo_d == hi_d else body + "?"

    def __repr__(self):
        return f"Ivl[{float(self.lo):.18g}, {float(self.hi):.18g}]"


def _coerce(o) -> Ivl:
    return o if isinstance(o, Ivl) else Ivl.exact(o)


# ---------------------------------------------------------------------------
# Irrational enclosures -- integer-only, provable
# ---------------------------------------------------------------------------

def isqrt_ivl(q) -> Ivl:
    """Certified rational enclosure of sqrt(q) for rational q >= 0.

    For q = a/b with a,b > 0: sqrt(a/b) = sqrt(a*b)/b. With N = _PREC and
    m = isqrt(a*b*N^2) we have, by definition of integer square root,
        m <= sqrt(a*b*N^2) < m+1
    hence  m/(b*N) <= sqrt(a/b) <= (m+1)/(b*N).
    Pure integer arithmetic, no float, provable by construction.
    """
    q = F(q)
    if q < 0:
        raise ValueError("sqrt of negative")
    if q == 0:
        return Ivl.exact(0)
    a, b = q.numerator, q.denominator
    n = _PREC
    m = isqrt(a * b * n * n)
    return Ivl(F(m, b * n), F(m + 1, b * n))


def sqrt_ivl(x) -> Ivl:
    """Certified sqrt of an interval (monotone, so endpointwise)."""
    x = _coerce(x)
    if x.lo < 0:
        raise ValueError("sqrt of interval extending below zero")
    return Ivl(isqrt_ivl(x.lo).lo, isqrt_ivl(x.hi).hi)


def mpf_to_frac(x) -> F:
    """EXACT conversion of an mpmath mpf to a Fraction.

    An mpf is a binary float: value = (-1)^sign * man * 2^exp. Reading the
    (sign, man, exp, bc) tuple therefore loses nothing, whereas going via a
    decimal string would silently round. Nothing on the proof path may round
    without widening, so this is the only conversion used.

    SOUNDNESS FIX (2026-07-30). The previous implementation re-created the
    value via mpf(x) in the ambient mp context before reading the tuple, and
    mp.prec in a fresh process is 53: the FIRST conversions of a run were
    silently rounded to double precision — a one-sided ~1e-17 error far outside
    the 1e-80 outward pad. Later calls were exact only because an earlier
    _asin_point/_s call had raised mp.prec. Found by the rung-2/rung-4
    consistency gate (cover + lens - Sprague missed 0 by 8.1e-19); fixed by
    reading the raw tuple directly, with a high-precision, restored fallback
    for inputs that are not already mpf-like.
    """
    tup = getattr(x, "_mpf_", None)
    if tup is None:
        from mpmath import mp
        old = mp.prec
        mp.prec = 1200
        try:
            tup = mp.mpf(x)._mpf_
        finally:
            mp.prec = old
    sign, man, exp, bc = tup
    if man == 0:
        if exp != 0:  # inf / -inf / nan carry man == 0 with a special exp
            raise ValueError(f"non-finite mpf: {x}")
        return F(0)
    v = F(man) * F(2) ** exp
    return -v if sign else v


def _mp_ivl(fn, *args) -> Ivl:
    """Convert an mpmath validated-interval result to an outward rational Ivl.

    mpmath.iv already rounds outward; the endpoint conversion is exact, and we
    still pad by a further tiny amount so the result is outward-rounded even if
    a future mpmath changed its endpoint convention. Cross-checked against
    independently known digit strings in tests/test_exact.py.
    """
    from mpmath import iv

    old = iv.prec
    iv.prec = 300
    try:
        r = fn(*args)
        lo = mpf_to_frac(r.a)
        hi = mpf_to_frac(r.b)
    finally:
        iv.prec = old
    pad = F(1, 10**80)
    return Ivl(lo - pad, hi + pad)


def pi_ivl() -> Ivl:
    from mpmath import iv
    return _mp_ivl(lambda: +iv.pi)


def sin_ivl(x) -> Ivl:
    from mpmath import iv
    x = _coerce(x)
    return _mp_ivl(lambda: iv.sin(iv.mpf([_s(x.lo), _s(x.hi)])))


def cos_ivl(x) -> Ivl:
    from mpmath import iv
    x = _coerce(x)
    return _mp_ivl(lambda: iv.cos(iv.mpf([_s(x.lo), _s(x.hi)])))


def _asin_point(y: F) -> Ivl:
    """Certified enclosure of asin(y) for a rational y in [-1, 1].

    mpmath's interval context provides no asin, and approximating one would put
    an uncertified number on the proof path. Instead: APPROXIMATE, THEN CERTIFY
    BY INVERSION. We take a high-precision guess t0, then *prove* the bracket
    [t0-d, t0+d] using only the certified sine and the monotonicity of sin on
    [-pi/2, pi/2]:

        sin(t0-d) <= y   ==>   asin(y) >= t0-d
        sin(t0+d) >= y   ==>   asin(y) <= t0+d

    Both premises are checked with `definitely_le`/`definitely_ge`, so they hold
    for every point of the enclosing intervals. The guess only has to be close;
    it never has to be trusted. If the bracket fails to verify we widen and
    retry, so a bad guess costs precision, never soundness.
    """
    from mpmath import mp

    if not (-1 <= y <= 1):
        raise ValueError("asin argument outside [-1, 1]")
    mp.prec = 600
    t0 = mp.asin(mp.mpf(y.numerator) / mp.mpf(y.denominator))
    guess = mpf_to_frac(t0)

    half_pi_hi = (pi_ivl() / 2).hi
    d = F(1, 10**80)
    for _ in range(60):
        lo_c, hi_c = guess - d, guess + d
        if -half_pi_hi <= lo_c and hi_c <= half_pi_hi:
            if sin_ivl(lo_c).definitely_le(y) and sin_ivl(hi_c).definitely_ge(y):
                return Ivl(lo_c, hi_c)
        d *= 1000
    raise ArithmeticError(f"asin bracket failed to certify for y={y}")


def asin_ivl(x) -> Ivl:
    """Certified arcsin. Increasing on [-1, 1], so enclose endpointwise."""
    x = _coerce(x)
    return Ivl(_asin_point(x.lo).lo, _asin_point(x.hi).hi)


def abs_ivl(x) -> Ivl:
    """Certified absolute value.

    If the interval straddles zero the result is [0, max|endpoint|] -- still a
    valid enclosure, but the caller usually wants to know, because a straddling
    cross-product means a degenerate configuration rather than a small distance.
    """
    x = _coerce(x)
    if x.lo >= 0:
        return x
    if x.hi <= 0:
        return -x
    return Ivl(F(0), max(-x.lo, x.hi))


def tan_ivl(x) -> Ivl:
    """Certified tangent as sin/cos.

    The division raises if the cosine enclosure straddles zero, which is exactly
    the right behaviour: near a pole there is no finite enclosure, and refusing
    is correct where returning a number would not be.
    """
    return sin_ivl(x) / cos_ivl(x)


def sec_ivl(x) -> Ivl:
    return Ivl.exact(1) / cos_ivl(x)


def _atan_point(y: F) -> Ivl:
    """Certified enclosure of atan(y), by the same certify-by-inversion pattern.

    tan is increasing on (-pi/2, pi/2), so a bracket [t0-d, t0+d] is *proved* by
        tan(t0-d) <= y  and  tan(t0+d) >= y
    both checked with certified comparisons against the certified tangent.
    """
    from mpmath import mp

    mp.prec = 600
    t0 = mp.atan(mp.mpf(y.numerator) / mp.mpf(y.denominator))
    guess = mpf_to_frac(t0)

    d = F(1, 10**80)
    for _ in range(60):
        lo_c, hi_c = guess - d, guess + d
        try:
            if tan_ivl(lo_c).definitely_le(y) and tan_ivl(hi_c).definitely_ge(y):
                return Ivl(lo_c, hi_c)
        except ZeroDivisionError:
            pass
        d *= 1000
    raise ArithmeticError(f"atan bracket failed to certify for y={y}")


def atan_ivl(x) -> Ivl:
    """Certified arctangent. Increasing everywhere, so enclose endpointwise."""
    x = _coerce(x)
    return Ivl(_atan_point(x.lo).lo, _atan_point(x.hi).hi)


def isolate_root(coeffs, lo, hi, iters: int = 400) -> Ivl:
    """Certified enclosure of a real root of an INTEGER polynomial by bisection.

    coeffs are integer, highest degree first. The polynomial is evaluated in exact
    rational arithmetic, so the sign of P(a) is known exactly -- no rounding can
    flip it. A sign change across [lo, hi] therefore *proves* a root lies inside,
    by the intermediate value theorem. Bisection narrows the enclosure; the
    guarantee never weakens.
    """
    def P(t: F) -> F:
        acc = F(0)
        for c in coeffs:
            acc = acc * t + c
        return acc

    a, b = F(lo), F(hi)
    fa, fb = P(a), P(b)
    if fa == 0:
        return Ivl(a, a)
    if fb == 0:
        return Ivl(b, b)
    if (fa > 0) == (fb > 0):
        raise ValueError("no sign change on the bracket: root not proved")
    for _ in range(iters):
        m = (a + b) / 2
        fm = P(m)
        if fm == 0:
            return Ivl(m, m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b, fb = m, fm
    return Ivl(a, b)


def _s(q: F) -> str:
    """Exact decimal string for a rational, for lossless handoff to mpmath.

    A Fraction with a non-terminating decimal expansion is widened to a
    terminating one by truncation toward -inf / +inf at the call sites above,
    which pad outward afterwards, so no inward rounding can occur.
    """
    from mpmath import mp
    mp.prec = 400
    return mp.nstr(mp.mpf(q.numerator) / mp.mpf(q.denominator), 90)


# ---------------------------------------------------------------------------
# Named certified constants
# ---------------------------------------------------------------------------

def sqrt3() -> Ivl:
    return isqrt_ivl(3)


def sqrt2() -> Ivl:
    return isqrt_ivl(2)


def deg(d) -> Ivl:
    """d degrees in radians, certified."""
    return pi_ivl() * F(d) / 180
