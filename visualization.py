import matplotlib.pyplot as plt
import numpy as np
import scipy

from comparison import compare_approximation_interpolation

plt.style.use("ggplot")


def visualization_px_with_fx(fx, px, xn, history, interval, n):
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

    plt.figure(figsize=(14, 5))
    plt.scatter(xn, Approx_err_xn, linewidth=1.5, label="Approximate Points")
    plt.plot(x, Approx_err_interval, linewidth=3, label="Approximate Error")
    plt.scatter(Interp_x, Interp_err_xn, linewidth=1.5, label="Interpolate Points")
    plt.plot(x, Interp_err_interval, linewidth=1.5, label="Interpolate Error")

    plt.xlabel("X Interval")
    plt.ylabel("Error")
    plt.title(f"Error on Interval {interval}, degree = {n}\nF(x)= {fx}")
    # plt.legend()
    # plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
    plt.legend(loc="center", bbox_to_anchor=(0.5, -0.25), ncol=4)
    plt.tight_layout()
    plt.show()

    print(f"f(x) = {fx}, interval = {interval}\n")
    print(f"polynomial degree = {n}")
    print(f"pn(x) = \n{px}\n")
    print(f"xn points:\n{xn}\n")
    print(f"converge iteration: {len(history['e'])}")
    print(f"MAE of approximation: {max(np.abs(Approx_err_interval))}")
    print(f"MAE of interpolation: {max(np.abs(Interp_err_interval))}")


def visualization(fx, px, xn, history, interval, n):
    x, xn, Interp_x, Approx_err_xn, Approx_err_interval, Interp_err_xn, Interp_err_interval = (
        compare_approximation_interpolation(fx, px, xn, interval, n)
    )

    plt.figure(figsize=(14, 5))
    plt.scatter(xn, Approx_err_xn, linewidth=1.5, label="Approximate Points")
    plt.plot(x, Approx_err_interval, linewidth=3, label="Approximate Error")
    plt.scatter(Interp_x, Interp_err_xn, linewidth=1.5, label="Interpolate Points")
    plt.plot(x, Interp_err_interval, linewidth=1.5, label="Interpolate Error")

    plt.xlabel("X Interval")
    plt.ylabel("Error")
    plt.title(f"Error on Interval {interval}, degree = {n}\nF(x)= {fx}")
    # plt.legend()
    # plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
    plt.legend(loc="center", bbox_to_anchor=(0.5, -0.25), ncol=4)
    plt.tight_layout()
    plt.show()

    print(f"f(x) = {fx}, interval = {interval}\n")
    print(f"polynomial degree = {n}")
    print(f"pn(x):\n{px}\n")
    print(f"xn points:\n{xn}\n")
    print(f"converge iteration: {len(history['e'])}")
    print(f"MAE of approximation: {max(np.abs(Approx_err_interval))}")
    print(f"MAE of interpolation: {max(np.abs(Interp_err_interval))}")
