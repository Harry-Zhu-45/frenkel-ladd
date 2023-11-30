"""
Author: HarryZ
Time: 2023-11-25
Function: Plot the mean square distance (msd) v.s. the spring constant (lambda)
"""
import matplotlib.pyplot as plt


lambda_list_1 = [
    0.001,
    0.01,
    0.1,
    1,
    10,
    100
]

msd_list_1 = [
    0.5132555,
    0.223856,
    0.143629,
    0.0268432,
    0.00703634,
    0.00279321
]

# lambda_list_2 = [0.001, 0.01, 0.1, 1, 10, 100]
# msd_list_2 = [0.0866683, 0.0667581, 0.0249357, 0.00866683, 0.000667581, 0.0000249357]

# plot
plt.figure(figsize=(8, 6))
plt.semilogx(lambda_list_1, msd_list_1, 'o-', label='msd')

plt.xlabel('lambda')
plt.ylabel('msd')
plt.legend()
# plt.show()
plt.savefig('msd.png')

plt.figure(figsize=(8, 6))
# plt.plot(lambda_list_1, msd_list_1, 'o-', label='msd')
# plt.plot(lambda_list_2, msd_list_2, 'o-', label='msd')
