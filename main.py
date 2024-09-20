from remez import *
from visualization import *

print("Example:\nFor polynomial: f(x) = 4x^4 + 2x^3 - 5x^2 + 8x - 2.5")
print("We want to find optimal uniform approximation of degree 3 on interval [-1, 1].\n")
print("It can be caculated theoretically as:")
print("pn(x) = 2x^3 - x**2 + 8x - 3\n")

print("Now let's calculate it using Remez algorithm:")
fx = "4 * x**4 + 2 * x**3 - 5 * x**2 + 8 * x - 2.5"
fx_der = "16 * x**3 + 6 * x**2 - 10 * x + 8"
interval = [-1, 1]
px, xn, history = remez(fx, fx_der, interval, n=3)
print(f"px =", px)

print("\nFurther, if we use degree of 4, we will get exactly the same polynomial as above.")
px, xn, history = remez(fx, fx_der, interval, n=4)
print(f"px =", px)

print("\nSame results as we increase degree s.t. it >= original.")
px, xn, history = remez(fx, fx_der, interval, n=6)
print(f"px =", px)

print("Visualization:")
visualization(fx, fx_der, interval, n=3, func_name=fx)
