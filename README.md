# Remez_Python

Python implementation of  [Remez algorithm](https://en.wikipedia.org/wiki/Remez_algorithm),  provides the best uniform polynomial approximation to functions based on the minimization of the maximum absolute error and the [equioscillation theorem](https://en.wikipedia.org/wiki/Equioscillation_theorem).

## Similar libarary

A similar C++ [library](https://github.com/samhocevar/lolremez?tab=readme-ov-file#docker) achieves the same accuracy as ours. MATLAB version can be found [here](https://ww2.mathworks.cn/matlabcentral/fileexchange/8094-remez-algorithm), but it often fails.

## Examples

Approximate `atan(sqrt(3+x³)-exp(1+x))` over the interval `[sqrt(2),pi²]` with 5-th degree polynomial:

```python
```







