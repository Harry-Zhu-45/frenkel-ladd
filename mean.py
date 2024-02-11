#!./venv/bin/python

"""
Calulate the mean (and standard deviation) of a set of data.
"""

import numpy as np


data = np.array([
    0.00233166,
    0.00232942,
    0.00233374,
    0.00233267,
    0.00233103,
    0.00233014,
    0.00233004,
    0.00232835,
    0.00233038,
    0.0023312
])

mean = np.mean(data)
# sd = np.std(data, ddof=1)

print('mean =', mean)
# print('sd =', sd)
