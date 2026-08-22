"""Exact rational and certified interval arithmetic.

Every quantity is an exact `Fraction` or an `Ivl` -- a closed interval with
exact rational endpoints guaranteed to contain the true value. Every operation
rounds outward, so containment survives composition. That is what makes a
computed inequality a proof: if `x.hi < 0` then the true value is negative,
with no floating-point caveat.

    from certivl import Ivl, pi_ivl, sqrt_ivl

    p = pi_ivl()
    assert p.lo < 314159265 / 100000000 < p.hi
"""
from .exact import (  # noqa: F401
    Ivl,
    abs_ivl,
    asin_ivl,
    atan_ivl,
    cos_ivl,
    deg,
    isolate_root,
    isqrt_ivl,
    mpf_to_frac,
    pi_ivl,
    sec_ivl,
    sin_ivl,
    sqrt2,
    sqrt3,
    sqrt_ivl,
    tan_ivl,
)

__version__ = "0.1.2"
__all__ = [
    "Ivl", "abs_ivl", "asin_ivl", "atan_ivl", "cos_ivl", "deg", "isolate_root",
    "isqrt_ivl", "mpf_to_frac", "pi_ivl", "sec_ivl", "sin_ivl", "sqrt2",
    "sqrt3", "sqrt_ivl", "tan_ivl",
]
