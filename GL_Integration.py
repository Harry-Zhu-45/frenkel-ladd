#!./venv/bin/python

"""
Author: HarryZ
Created: 2023-11-04
Updated: 2023-12-02
Calulate the integral of a function using Gaussian-Legendre quadrature
"""

import numpy as np


####### 5-point G-L quadrature #######

weight = [
    0.2369268850561891,
    0.4786286704993665,
    0.5688888888888889,
    0.4786286704993665,
    0.2369268850561891
]

abscissa = [
    -0.9061798459386640,
    -0.5384693101056831,
    0.0000000000000000,
    0.5384693101056831,
    0.9061798459386640
]

lambda_list = [
    5.00421035,
    33.0590995,
    115.297688,
    299.738483,
    544.708655
]

msd_list = [
    0.014400099999999999,
    0.01058859052631579,
    0.0064644305,
    0.003898474545454546,
    0.0025043180000000006
]


res = 0.0
for i in range(5):
    # print(msd_list[i] * (lambda_list[i] + np.e ** 3.5))
    res += weight[i] * msd_list[i] * (lambda_list[i] + np.e ** 3.5)


print(-res)
