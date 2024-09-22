import numpy as np
import scipy


def compare_approximation_interpolation(fx, px, xn, interval, n):
    def f(x):
        return eval(fx)

    # interval points
    a, b = interval
    x = np.linspace(a, b, 200)

    # approximation in n + 2 points
    def Approx_func(x):
        return eval(px)

    Approx_err_xn = [f(i) - Approx_func(i) for i in xn]
    Approx_err_interval = [f(j) - Approx_func(j) for j in x]

    # Interpolation in n + 1 points
    Interp_x = np.array(
        [
            (a + b + (b - a) * np.cos(((2 * k + 1) * np.pi) / (2 * (n + 1)))) / 2
            for k in range(n + 1)
        ]
    )
    Interp_y = np.array([f(i) for i in Interp_x])
    interpolate = scipy.interpolate.BarycentricInterpolator(Interp_x, Interp_y)

    def Interp_func(x):
        return interpolate(x)

    Interp_err_xn = [f(i) - Interp_func(i) for i in Interp_x]
    Interp_err_interval = [f(j) - Interp_func(j) for j in x]

    return x, xn, Interp_x, Approx_err_xn, Approx_err_interval, Interp_err_xn, Interp_err_interval
