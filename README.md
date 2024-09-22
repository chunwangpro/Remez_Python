# Remez_Python

Python implementation of  [Remez algorithm](https://en.wikipedia.org/wiki/Remez_algorithm),  provides the best uniform polynomial approximation to functions based on the minimization of the maximum absolute error and the [equioscillation theorem](https://en.wikipedia.org/wiki/Equioscillation_theorem).

## Similar libarary

A similar C++ [library](https://github.com/samhocevar/lolremez?tab=readme-ov-file#docker) achieves the same accuracy as ours. MATLAB version can be found [here](https://ww2.mathworks.cn/matlabcentral/fileexchange/8094-remez-algorithm), but it often fails.

## Usage

Write the math expression of F(x) and its derivative in string like Numpy or Sympy style, then set approximation interval and polynomial degrees.

```python
fx = "np.sin(x)"
fx_der = "np.cos(x)"
interval = [-10, 10]
n = 10
px, xn, history = remez(fx, fx_der, interval, n)
visualization(fx, px, xn, history, interval, n, compare_interpolation=True)
```



## Examples

Approximate `atan(sqrt(3+x³)-exp(1+x))` over the interval `[sqrt(2),pi²]` with 5-th degree polynomial:

```python
fx = "np.arctan(np.sqrt(3 + x**3) - np.exp(1 + x))"
g = "(np.sqrt(3 + x**3) - np.exp(1 + x))"
g_prime = "((3 * x**2) / (2 * np.sqrt(3 + x**3)) - np.exp(1 + x))"
fx_der = f"({g_prime}) / (1 + {g}**2)"
interval = [np.sqrt(2), np.pi**2]
n = 5
px, xn, history = remez(fx, fx_der, interval, n)
visualization(fx, px, xn, history, interval, n)
```

Results:

```bash
f(x) = np.arctan(np.sqrt(3 + x**3) - np.exp(1 + x)), interval = [1.4142135623730951, 9.869604401089358]

polynomial degree = 5
pn(x):
- 3.955756933047265e-05 * x**5 + 0.0012947712130833584 * x**4 - 0.01654139703555944 * x**3 + 0.10351664953941357 * x**2 - 0.32051562487328494 * x - 1.1703528319321932

xn points:
[1.41421356 1.83693284 3.14845216 5.17561358 7.42752857 9.19850669
 9.8696044 ]

converge iteration: 7
MAE of approximation: 0.0012079008992569307
MAE of interpolation: 0.0021889162615582602
```







