from matrix_solver import *


def remez(fx, fx_der, interval, n, max_iters=50, tol=1e-30):
    """
    The iterative approach is adopted, progressively adjusting the extremal points to find the alternation point set (xn) and the approximation polynomial of given degree (px) that minimize the max absolute error (at given interval).

    Parameters
    ----------
    fx : string-like
        The explicit fomula of the function to be approximated (limited to Numpy math grammar).
    fx_der : string-like
        The derivative of function fx (limited to Numpy math grammar).
    n : int
        The degree of the approximation polynomial, must be positive.
    interval : list of two floats
        A 2-element list representing the left and right end of interval over which the polynomial approximation will be evaluated.

    Returns
    -------
    px : str
        A string representing the optimal polynomial approximation in the form:
        "a0 + a1*x + a2*x^2 + ... + an*x^n", where n is the input degree.
    xn : array-like
        The final alternation points set where the error is maximized.
    history : dict
        A dictionary containing the history of the iterative process with keys:
        - 'px': list of polynomial approximations at each iteration.
        - 'xn': list of alternation points at each iteration.
        - 'ind': list of indices of alternation point corresponding to the largest error at each iteration.
        - 'e': list of the absolute errors at each iteration.

    Examples
    --------
    >>> fx = "np.sin(x)"
    >>> fx_der = "np.cos(x)"
    >>> n = 2
    >>> interval = [-1, 1]
    >>> px, xn, history = remez(fx, fx_der, interval, n)
    """

    # Initial inference of xn (n + 2 points):
    a, b = interval
    Chebyshev_points = [
        (a + b + (b - a) * np.cos(((2 * k + 1) * np.pi) / (2 * n))) / 2 for k in range(n)
    ]
    xn = np.array(interval[:1] + Chebyshev_points + interval[1:])
    xn.sort()

    history = {"an": [], "xn": [], "ind": [], "e": []}
    for _ in range(max_iters):
        try:
            an, xn_new, ind, e = matrix_solver(fx, fx_der, xn, interval, n)
            for key, value in zip(history.keys(), [an, xn, ind, abs(e)]):
                history[key].append(value)
            if abs(xn_new[ind] - xn[ind]) < tol:
                break
            if ind < (len(xn) - 1) and abs(xn_new[ind] - xn[ind + 1]) < tol:
                break
            xn = xn_new
        except np.linalg.LinAlgError:
            # In case matrix become singular after many iterations and matrix_solver fails
            xn = history["xn"][-1]
            break
    px = ensemble_polynomial(an, smart_round=False, keep_first_zeros=True)
    return px, xn, history
