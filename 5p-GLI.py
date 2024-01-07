#!./venv/bin/python

"""
Author: HarryZ
Created: 2023-11-04
Updated: 2024-01-07
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

# 第一种取点
# lambda_min = 0
# lambda_max = 632.026
# lambda_list_5 = [(lambda_max - lambda_min) / 2 * i +
#                   (lambda_max + lambda_min) / 2 for i in abscissa_5]
# print(lambda_list_5)

lambda_list_5 = [
    29.6483883,
    145.849698,
    316.013000,
    486.176302,
    602.377612
]

msd_list_5 = [
    0.01795073,
    0.006699942,
    0.004176842,
    0.002769681,
    0.002439662
]

# 第二种取点
# lambda_min = np.log(lambda_min + np.e ** 3.5)
# lambda_max = np.log(lambda_max + np.e ** 3.5)
# tmp = [(lambda_max - lambda_min) / 2 * i +
#                   (lambda_max + lambda_min) / 2 for i in abscissa_5]
# lambda_list_5 = [np.e ** i - np.e ** 3.5 for i in tmp]
# print(lambda_list_5)

# lambda_list_5 = [
#     5.00421035,
#     33.0590995,
#     115.297688,
#     299.738483,
#     544.708655
# ]

# msd_list_5 = [
#     0.014400099999999999,
#     0.01058859052631579,
#     0.0064644305,
#     0.003898474545454546,
#     0.0025043180000000006
# ]

res = 0.0
for i in range(5):
    print(msd_list_5[i] * (lambda_list_5[i] + np.e ** 3.5))
    res += weight_5[i] * msd_list_5[i] * (lambda_list_5[i] + np.e ** 3.5)

print(-res)
