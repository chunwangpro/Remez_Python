import matplotlib.pyplot as plt
import numpy as np

from comparison import compare_approximation_interpolation

plt.style.use("ggplot")


def visualization(fx, px, xn, history, interval, n, compare_interpolation=True):
    x, xn, Interp_x, Approx_err_xn, Approx_err_interval, Interp_err_xn, Interp_err_interval = (
        compare_approximation_interpolation(fx, px, xn, interval, n)
    )

    plt.figure(figsize=(14, 5))
    plt.scatter(xn, Approx_err_xn, linewidth=1.5, label="Approximate Points")
    plt.plot(x, Approx_err_interval, linewidth=3, label="Approximate Error")
    if compare_interpolation:
        plt.scatter(Interp_x, Interp_err_xn, linewidth=1.5, label="Interpolate Points")
        plt.plot(x, Interp_err_interval, linewidth=1.5, label="Interpolate Error")

    plt.xlabel("X Interval")
    plt.ylabel("Error")
    plt.title(f"Error on Interval {interval}, degree = {n}\nF(x)= {fx}")
    plt.legend(loc="center", bbox_to_anchor=(0.5, -0.25), ncol=4)
    plt.tight_layout()
    plt.savefig(f"./images/single_plot/{fx}_{interval}_{n}.png", dpi=300)
    plt.show()

    print(f"f(x) = {fx}, interval = {interval}\n")
    print(f"polynomial degree = {n}")
    print(f"pn(x):\n{px}\n")
    print(f"xn points:\n{xn}\n")
    print(f"converge iteration: {len(history['e'])}")
    print(f"MAE of approximation: {max(np.abs(Approx_err_interval))}")
    print(f"MAE of interpolation: {max(np.abs(Interp_err_interval))}")


def visualization_pipeline(fx, px, xn, history, interval, n, history_error=None):
    """
    Visualization for plenty of functions, auto-pipeline, including a 'whether continue' to decide whether to plot.
    """
    x, xn, Interp_x, Approx_err_xn, Approx_err_interval, Interp_err_xn, Interp_err_interval = (
        compare_approximation_interpolation(fx, px, xn, interval, n)
    )

    # whether continue to plot
    max_interval_error = np.max(Approx_err_interval)
    thershold = np.array([history_error, max(Approx_err_xn), min(Approx_err_interval)])
    new = np.array([max_interval_error, max_interval_error, min(Approx_err_xn)])
    whether_continue = (new <= thershold).all()

    if whether_continue:
        plt.figure(figsize=(14, 5))
        plt.scatter(xn, Approx_err_xn, linewidth=1.5, label="Approximate Points")
        plt.plot(x, Approx_err_interval, linewidth=3, label="Approximate Error")
        plt.scatter(Interp_x, Interp_err_xn, linewidth=1.5, label="Interpolate Points")
        plt.plot(x, Interp_err_interval, linewidth=1.5, label="Interpolate Error")

        plt.xlabel("X Interval")
        plt.ylabel("Error")
        plt.title(f"Error on Interval {interval}, degree = {n}\nF(x)= {fx}")
        plt.legend(loc="center", bbox_to_anchor=(0.5, -0.25), ncol=4)
        plt.tight_layout()
        plt.savefig(f"./images/pipeline_plot/{fx}_{interval}_{n}.png", dpi=300)
        plt.show()

        print(f"f(x) = {fx}, interval = {interval}\n")
        print(f"polynomial degree = {n}")
        print(f"pn(x):\n{px}\n")
        print(f"xn points:\n{xn}\n")
        print(f"converge iteration: {len(history['e'])}")
        print(f"MAE of approximation: {max(np.abs(Approx_err_interval))}")
        print(f"MAE of interpolation: {max(np.abs(Interp_err_interval))}")
    return max_interval_error if whether_continue else None


def visualization_px_with_fx(fx, px, xn, history, interval, n):
    def f(x):
        return eval(fx)

    x, xn, Interp_x, Approx_err_xn, Approx_err_interval, Interp_err_xn, Interp_err_interval = (
        compare_approximation_interpolation(fx, px, xn, interval, n)
    )

    _, ax = plt.subplots(1, 2, figsize=(14, 5))
    ax = ax.flatten()
    ax[0].scatter(xn, Approx_err_xn, linewidth=1.5, label="Approximate Points")
    ax[0].plot(x, Approx_err_interval, linewidth=1.5, label="Approximate Error")
    ax[0].scatter(Interp_x, Interp_err_xn, linewidth=1.5, label="Interpolate Points")
    ax[0].plot(x, Interp_err_interval, "-", label="Interpolate Error")

    ax[0].set(
        xlabel="X Interval",
        ylabel="Error",
        title=f"Error on Interval {interval}, degree = {n}\nF(x)= {fx}",
    )
    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width, box.height * 0.8])
    ax[0].legend(loc="center", bbox_to_anchor=(0.5, 1.2), ncol=4)

    # err_list = history["e"]
    # ax[1].plot(range(len(err_list)), err_list, "d--", linewidth=1.5)
    # ax[1].set(
    #     xlabel="Epoch",
    #     ylabel="Max Error",
    #     title=f"Error On Points In Every Epoches",
    #     xticks=np.linspace(1, 10, 10),
    #     xlim=(0.5, len(err_list) + 0.5),
    # )

    ax[1].plot(x, [f(j) for j in x], linewidth=1.5)
    ax[1].set(
        xlabel="X Interval",
        ylabel="Y",
        title=f"F(x)",
        # xticks=np.linspace(1, 10, 10),
        # xlim=(0.5, len(err_list) + 0.5),
    )
    plt.tight_layout()
    plt.savefig(f"./images/compare_plot/{fx}_{interval}_{n}.png", dpi=300)
    plt.show()

    print(f"f(x) = {fx}, interval = {interval}\n")
    print(f"polynomial degree = {n}")
    print(f"pn(x) = \n{px}\n")
    print(f"xn points:\n{xn}\n")
    print(f"converge iteration: {len(history['e'])}")
    print(f"MAE of approximation: {max(np.abs(Approx_err_interval))}")
    print(f"MAE of interpolation: {max(np.abs(Interp_err_interval))}")
