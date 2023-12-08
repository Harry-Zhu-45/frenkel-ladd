#!./venv/bin/python

"""
Author: HarryZ
Created: 2023-12-08
Updated: 2023-12-08
Plot the histogram of the data.
"""

import numpy as np
import matplotlib.pyplot as plt


####### histogram #######

data = [
    0.0160407,
    0.0193717,
    0.0167812,
    0.0180884,
    0.0145106,
    0.0146094,
    0.0152368,
    0.018348,
    0.0166659,
    0.0161771,
    0.0137945,
    0.018463,
    0.0141656,
    0.0233901,
    0.0149204,
    0.0165738,
    0.0204618,
    0.0159103,
    0.0156607,
    0.0214311
]

data = np.array(data)
mean = np.mean(data)
sd = np.std(data, ddof=1)

print('mean =', mean)
print('sd =', sd)

plt.hist(data,
         #  bins=10,
         bins='auto',
         label='CDF')
plt.show()
