#!./venv/bin/python

"""
Calulate the integral of a function using 10-point Gaussian-Legendre quadrature
"""

import numpy as np


####### 10-point G-L quadrature #######

abscissa_10 = [
    -0.9739065285,
    -0.8650633667,
    -0.6794095683,
    -0.4333953941,
    -0.1488743390,
    0.1488743390,
    0.4333953941,
    0.6794095683,
    0.8650633667,
    0.9739065285
]

weight_10 = [
    0.0666713443,
    0.1494513492,
    0.2190863625,
    0.2692667193,
    0.2955242247,
    0.2955242247,
    0.2692667193,
    0.2190863625,
    0.1494513492,
    0.0666713443
]


origin_lambda_min = 0
origin_lambda_max = 632.026

lambda_min = np.log(origin_lambda_min + np.e ** 3.5)
lambda_max = np.log(origin_lambda_max + np.e ** 3.5)

tmp = [(lambda_max - lambda_min) / 2 * i +
       (lambda_max + lambda_min) / 2 for i in abscissa_10]
lambda_list_10 = [np.e ** i - np.e ** 3.5 for i in tmp]
# print(lambda_list_10)

lambda_list_10 = [
    1.3218454386828071,
    7.429242462300721,
    20.448998156279607,
    44.35579276774375,
    85.59512159890338,
    152.4321293389935,
    251.20242463993492,
    378.0986975526874,
    510.1482256352902,
    606.495129089744
]

msd_list_10 = [
    0.0162496,
    0.01406634,
    0.01168701,
    0.009454013,
    0.007437346,
    0.005726951,
    0.004358547,
    0.003346746,
    0.0026830330000000005,
    0.002330863
]

res = 0.0
for i in range(10):
    res += weight_10[i] * msd_list_10[i] * (lambda_list_10[i] + np.e ** 3.5)
res *= (lambda_max - lambda_min) / 2

print(-res)
