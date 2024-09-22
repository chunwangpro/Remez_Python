# Numerical Analysis: Best Uniform Approximation


"""
This is the original version of the Remez algorithm, which is deprecated.
Please refer to the Remez.py for the latest version.
"""


"""
Tuesday Dec/17/2019 12:41:25

@author: Chun Wang, School of Mathematics & Statistics, Wu Han University, China.
"""
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import sympy as sy
from scipy.interpolate import BarycentricInterpolator
from scipy.optimize import minimize_scalar


def matrix_solver(fx, fx_der, xn, n, interval):
    """
    Solve the linear system of equations: pn(xi) + (-1)^i*E = f(xi), i = 0, 1,..., n+1.
    """
    # symbolize x and create the function f(x)
    x = sy.Symbol("x")

    def f(x):
        return eval(fx)

    # create the linear system of equations and get the solution
    P = np.array([[xi**i for i in range(n + 1)] for xi in xn])
    E = np.array([[(-1) ** i] for i in range(n + 2)])
    A = np.concatenate((P, E), axis=1)
    F = [[f(xi)] for xi in xn]
    an = sp.linalg.solve(A, F)
    # get the string-like polynomial function px of the Best Uniform Approximation
    # e is the value of error and an is the coefficients only
    e, an = an[-1, 0], an[:-1, 0]
    px = f"{an[0]}"
    for i in range(1, n + 1):
        px += f"+{an[i]}*x**{i}"

    # define the error function
    def err(x):
        return eval(fx) - eval(px)

    # the first and last element of xn_root is the start point and end point of the interval
    # in between we fill with the roots of the error function
    xn_root = [
        minimize_scalar(
            lambda x: abs(err(x)),
            bounds=(xn[i], xn[i + 1]),
            method="bounded",
            options={"xatol": 1e-52},
        ).x
        for i in range(n + 1)
    ]
    xn_root = np.concatenate(([interval[0]], xn_root, [interval[1]]), axis=0)
    # define the derivative of error
    minus_px_der = str(sy.diff(-eval(px), x))

    def err_der(x):
        return eval(fx_der) + eval(minus_px_der)

    # find the extrem points according to the sign of derivative of error in adjacent points
    extrem_point, value = [-np.inf] * (n + 2), [np.inf] * (n + 2)
    for i in range(n + 2):
        # if they have different sign, find the extrem point in between
        if err_der(xn_root[i]) * err_der(xn_root[i + 1]) < 0:
            extrem_point[i] = minimize_scalar(
                lambda x: abs(err_der(x)),
                bounds=(xn_root[i], xn_root[i + 1]),
                method="bounded",
                options={"xatol": 1e-52},
            ).x
        else:  # choose the root with greater absolute value of error
            extrem_point[i] = (
                xn_root[i] if abs(err(xn_root[i])) > abs(err(xn_root[i + 1])) else xn_root[i + 1]
            )
        value[i] = abs(err(extrem_point[i]))
    return px, xn, extrem_point, np.argmax(value), e


def remez(fx, fx_der, n, interval, raise_exception=True):
    """
    fx: (string-like) the function that will be given the polynomials of best uniform approximation of.
    fx_der: (string-like) the derivative of function fx.
    n: (positive int) the degree of the polynomial.
    interval: (list) a 2-list
    so pn(x) is like (a0 + a1*x + a2*x^2 + ...+ an*x^n).
    """
    # the first choice of the (n+2) crossing-points
    # created by Chebyshev roots (n) plus the start and end 2 points of interval
    a, b = interval[0], interval[1]
    xn = np.array(
        [a]
        + [(a + b + (b - a) * np.cos(((2 * k + 1) * np.pi) / (2 * n))) / 2 for k in range(n)]
        + [b]
    )
    # the following points perform bad:
    # (bad) equally distance points: xn = np.linspace(-1, 1, n+2)
    xn.sort()  # make the points array sorted ascendingly
    epoches = 1  # count times of iteration
    record = {"epoch": [], "px": [], "xn": [], "extrem_point": [], "ind": [], "e": []}
    while epoches < 11:  # max times of iteration
        try:
            px, xn, extrem_point, ind, e = matrix_solver(fx, fx_der, xn, n, interval)
            # update the record dictionary
            record["epoch"].append(epoches)
            record["px"].append(px)
            record["xn"].append(xn)
            record["extrem_point"].append(extrem_point)
            record["ind"].append(ind)
            record["e"].append(abs(e))
        except np.linalg.LinAlgError:
            # sometimes may fall in a singular matrix after many iterations
            # then just return the params in last iteration
            if raise_exception:
                print(
                    "\n"
                    + "#" * 10
                    + f" Next {fx} Matrix is singular! (n = {n}, epoch = {epoches}) "
                    + "#" * 10
                )
            return (
                record["px"][-1],
                record["xn"][-1],
                record["epoch"][-1],
                record["epoch"],
                record["e"],
            )
        else:
            # if difference between new point and the corressponding one in old array
            # is less than a certain threshold, then quit
            # also compare it with the following point if it is not the last point
            if abs(xn[ind] - extrem_point[ind]) < 1e-30:
                break
            elif ind < (len(xn) - 1) and abs(xn[ind + 1] - extrem_point[ind]) < 1e-30:
                break
            else:
                xn = extrem_point  # replace the old points with the new points
                epoches += 1  # update the time count of iteration
    return record["px"][-1], record["xn"][-1], record["epoch"][-1], record["epoch"], record["e"]


if __name__ == "__main__":

    def err_analysis(fx, fx_der, n, interval, history_error, plot_all=True):
        message = False  # # True if plot_all else False
        tic = time.time()
        px, xn, epoches, epoch_list, err_list = remez(
            fx, fx_der, n, interval, raise_exception=message
        )
        toc = time.time()

        def f(x):
            return eval(fx)

        def err(x):
            return eval(fx) - eval(px)

        # approximation
        x = np.linspace(interval[0], interval[1], 200)
        interval_err, xn_err = [err(i) for i in x], [err(j) for j in xn]
        # whether continue
        max_interval_error = max(interval_err)
        thershold = np.array([history_error, max(xn_err), min(interval_err)])
        new = np.array([max_interval_error, max_interval_error, min(xn_err)])
        whether_continue = (new <= thershold).all()
        if plot_all:
            whether_plot = True
        else:
            if plot_example:
                whether_plot = True if name == "poly" else False
            else:
                whether_plot = False
        if whether_continue and whether_plot:
            fig, ax = plt.subplots(1, 2, figsize=(14, 4))
            ax = ax.flatten()
            ax[0].plot(x, interval_err, linewidth=1.5, label="Approximate Error")
            ax[0].scatter(xn, xn_err, linewidth=1.5, label="Approximate Points")
            # Interpolate in n + 1 points
            a, b = interval[0], interval[1]
            xn_inter = np.array(
                [
                    (a + b + (b - a) * np.cos(((2 * k + 1) * np.pi) / (2 * (n + 1)))) / 2
                    for k in range(n + 1)
                ]
            )
            y_inter = np.array([f(xi) for xi in xn_inter])
            interpolate = BarycentricInterpolator(xn_inter, y_inter)
            ax[0].plot(x, f(x) - interpolate(x), "-", label="Interpolate Error")
            ax[0].plot(
                xn_inter,
                f(xn_inter) - interpolate(xn_inter),
                "o",
                linewidth=1.5,
                label="Interpolate Points",
            )
            ax[0].set(
                xlabel="X Interval",
                ylabel="Error",
                title=f"Error On Interval {interval} Of {name} (n = {n})",
            )
            box = ax[0].get_position()
            ax[0].set_position([box.x0, box.y0, box.width, box.height * 0.8])
            ax[0].legend(loc="center", bbox_to_anchor=(0.5, 1.2), ncol=4)
            ax[1].plot(epoch_list, err_list, "d--", linewidth=1.5)
            ax[1].set(
                xlabel="Epoch",
                ylabel="Max Error",
                title=f"Error On Points In Every Epoches",
                xticks=np.linspace(1, 10, 10),
                xlim=(0.5, len(epoch_list) + 0.5),
            )
            plt.tight_layout()
            plt.show()
            print(
                f"f(x) = {fx}, n = {n}, iter epoches: {epoches}, time collapsed: {toc - tic:.2f}s\npn(x) = {px}"
            )
            print(
                f"max approximation error: {max_interval_error}\nmax interpolation error: {max(abs(f(x)-interpolate(x)))}"
            )
            print(f"points: {xn}" + "\n" + "##" * 70)
        return max_interval_error if whether_continue else None

    # set some hyper-parameter
    interval = [-1, 1]
    degree = 40
    plot_all = True  # whether want to plot all the figures
    plot_example = True  # whether plot the example
    plot_df = True  # whether want to plot the dataframe
    # first comes a polynomial example to verify the correctness and reliability of the algorithm
    # then give more examples of smooth and nonsmooth functions
    half = degree // 2
    fun_dict = {
        "poly": [("4*x**4+2*x**3-5*x**2+8*x-2.5", "16*x**3+6*x**2-10*x+8", 3)],
        r"$\sin{x}$": [("np.sin(x)", "np.cos(x)", 2 * i) for i in range(1, half)],
        r"$e^x$": [("np.exp(x)", "np.exp(x)", 2 * i + 1) for i in range(half)],
        r"$e^{1+x^2}$": [("np.exp(1+x**2)", "2*x*np.exp(1+x**2)", 2 * i + 1) for i in range(half)],
        r"$\ln{(1+x^2)}$": [("np.log(1+x**2)", "2*x/(1+x**2)", 2 * i + 1) for i in range(half)],
        r"$\sqrt{1+x^2}$": [
            ("np.sqrt(1+x**2)", "x/np.sqrt(1+x**2)", 2 * i + 1) for i in range(half)
        ],
        r"$\frac{1}{1+x^2}$": [
            ("(1+x**2)**-1", "-2*x*(1+x**2)**-2", 2 * i + 1) for i in range(half)
        ],
    }

    df = pd.DataFrame([2 * i + 1 for i in range(half)], columns=["n"])
    df.set_index(keys="n", inplace=True)
    for name, function in fun_dict.items():
        store, history_error = [], np.inf
        for fx, fx_der, n in function:
            if name == "poly" and (plot_example or plot_all):
                print(
                    "Exsample:\nFor n = 3, interval = [-1, 1] and polynomial: f(x) = 4x^4 + 2x^3 - 5x^2 + 8x - 2.5"
                )
                print(
                    "The exact pn can be caculated theoretically as: pn(x) = -3 + 8x - x^2 + 2x^3"
                )
                print(
                    "If the interval is not [-1, 1], you have to do some transformations to verify."
                )
            max_interval_error = err_analysis(
                fx, fx_der, n, interval, history_error, plot_all=plot_all
            )
            if max_interval_error is not None:
                store.append((n, max_interval_error))
                history_error = max_interval_error
            else:
                if name != r"$\sin{x}$":
                    y = pd.DataFrame(
                        [j[1] for j in store], columns=[name], index=[j[0] for j in store]
                    )
                    df = pd.merge(df, y, left_index=True, right_index=True, how="left")
                break
    df.dropna(how="all", inplace=True)
    x_length = df.index[-1]
    if plot_df or plot_all:
        df.plot(
            figsize=(12, 6),
            title=f"log(Max Error) Versus Degree Increase With Interval {interval}",
            logy=True,
            style=["rd-", "d--", "d--", "d--", "d--"],
            xlim=(0.5, x_length + 0.5),
            ylim=[5 * 10**-16, 1],
            xticks=np.arange(1, x_length + 1, 2),
            yticks=[10**-i for i in range(16)],
            fontsize=12,
        )
        plt.show()
        print(df)
