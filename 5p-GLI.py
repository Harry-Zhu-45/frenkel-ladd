#!./venv/bin/python

"""
Calulate the integral of a function using 5-point Gaussian-Legendre quadrature
"""

import numpy as np


####### 5-point G-L quadrature #######

abscissa_5 = [
    -0.9061798459386640,
    -0.5384693101056831,
    0.0000000000000000,
    0.5384693101056831,
    0.9061798459386640
]

weight_5 = [
    0.2369268850561891,
    0.4786286704993665,
    0.5688888888888889,
    0.4786286704993665,
    0.2369268850561891
]

origin_lambda_min = 0
origin_lambda_max = 632.026

lambda_min = np.log(origin_lambda_min + np.e ** 3.5)
lambda_max = np.log(origin_lambda_max + np.e ** 3.5)
tmp = [(lambda_max - lambda_min) / 2 * i +
       (lambda_max + lambda_min) / 2 for i in abscissa_5]
lambda_list_5 = [np.e ** i - np.e ** 3.5 for i in tmp]
# print(lambda_list_5)

lambda_list_5 = [
    5.00421035,
    33.0590995,
    115.297688,
    299.738483,
    544.708655
]

msd_list_5 = [
    0.014400099999999999,
    0.01058859052631579,
    0.0064644305,
    0.003898474545454546,
    0.0025043180000000006
]

res = 0.0
for i in range(5):
    res += weight_5[i] * msd_list_5[i] * (lambda_list_5[i] + np.e ** 3.5)
res *= (lambda_max - lambda_min) / 2

print(-res)
