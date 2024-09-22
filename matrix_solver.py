import warnings

import numpy as np
import scipy
import sympy as sp

from utils import *

warnings.filterwarnings("ignore")


def matrix_solver(fx, fx_der, xn, interval, n):
    """
    Solve the system of linear equations: pn(xi) + (-1)^i * E = f(xi), where i = 0, 1,..., n+1.

    The goal is to minimize the absolute error of the polynomial approximation for the current alternation point set (xn) and by calculating the residual between the function (fx) and the approximation polynomial (px), identify a new alternation point set (xn_new), thereby further optimizing the approximation in subsequent steps.

    Parameters
    ----------
    fx : string-like
        The function to be approximated by polynomials.
    fx_der : string-like
        The derivative of function fx.
    xn : array-like
        Array of `n+2` points where the function and polynomial will be evaluated.
    n : int
        The degree of the approximation polynomial, must be positive.
    interval : list of two floats
        A 2-element list representing the interval over which the polynomial approximation will be evaluated.

    Returns
    -------
    an : list
        List of `n+1` coeffients of polynomial expression.
    xn_new : array-like
        The new set of `n+2` extremal points where the maximum absolute error occurs.
    ind : int
        The index of the alternation point set with the largest absolute error.
    e : float
        The value of system error term `E`.
    """

    x = sp.Symbol("x")
    f = lambda x: eval(fx)

    P = np.array([[xi**i for i in range(n + 1)] for xi in xn])
    E = np.array([[(-1) ** i] for i in range(n + 2)])
    A = np.concatenate((P, E), axis=1)
    F = [[f(xi)] for xi in xn]
    an = scipy.linalg.solve(A, F)
    e, an = an[-1, 0], an[:-1, 0]

    # form the polynomial approximation function px
    # px = f"{an[0]}" + "".join([f" + {an[i]} * x**{i}" for i in range(1, n + 1)])
    # or more elegant
    px = ensemble_polynomial(an, use_smart_round=False, keep_first_zeros=True)

    # error function and its derivative, do not touch this line
    err = lambda x: eval(fx) - eval(px)
    neg_px_der = str(sp.diff(-eval(px), x))
    err_der = lambda x: eval(fx_der) + eval(neg_px_der)

    # find minimum of abs(error) / roots of the error function
    xn_root = [
        scipy.optimize.minimize_scalar(
            lambda x: abs(err(x)),
            bounds=(xn[i], xn[i + 1]),
            method="bounded",
            options={"xatol": 1e-52},
        ).x
        for i in range(n + 1)
    ]
    xn_root = np.array(interval[:1] + xn_root + interval[1:])

    # find xn_new if different sign of error derivative appear / choose the root with greater absolute error
    xn_new = np.zeros(n + 2)
    for i in range(n + 2):
        if err_der(xn_root[i]) * err_der(xn_root[i + 1]) < 0:
            xn_new[i] = scipy.optimize.minimize_scalar(
                lambda x: abs(err_der(x)),
                bounds=(xn_root[i], xn_root[i + 1]),
                method="bounded",
                options={"xatol": 1e-52},
            ).x
        else:
            xn_new[i] = (
                xn_root[i] if abs(err(xn_root[i])) > abs(err(xn_root[i + 1])) else xn_root[i + 1]
            )
    ind = np.argmax([abs(err(xn_new[i])) for i in range(n + 2)])
    return an, xn_new, ind, e
