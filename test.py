#!./venv/bin/python

"""
Author: HarryZ
Created: 2023-11-25
Updated: 2023-12-02
Function: Plot the mean square distance (msd) v.s. the spring constant (lambda)
"""

import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.rcParams["font.size"] = 10
mpl.rcParams['xtick.direction'] = 'in'  # 将 x-axis 的刻度线方向设置向内
mpl.rcParams['ytick.direction'] = 'in'  # 将 y-axis 的刻度线方向设置向内

mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10

mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

mpl.rcParams["legend.markerscale"] = 1

lambda_list_1 = [
    0.001,
    0.01,
    0.1,
    1,
    10,
    100
]

msd_list_1 = [
    0.01456625,
    0.014531067,
    0.0145731,
    0.01365913,
    0.0127761625,
    0.006901662
]

fig, axs = plt.subplots(
    # ncols=2,
    figsize=(8, 6),
    # figsize=(16, 6),
    layout="constrained"
)

plt.semilogx(lambda_list_1, msd_list_1, 'o-', label='msd')
# axs[1].semilogx(lambda_list_2, msd_list_2, 'o-', label='msd')

plt.xlabel('lambda')
plt.ylabel('msd')
plt.ylim(0, 0.02)

plt.errorbar(0.001, 0.01456625, yerr=0.0009041842456048431)
plt.errorbar(0.01, 0.014531067, yerr=0.0004404669832499741)
plt.errorbar(0.1, 0.0145731, yerr=0.0012012213659438457)
plt.errorbar(1, 0.01365913, yerr=0.0014643712723668588)
plt.errorbar(10, 0.0127761625, yerr=0.0010252721393958918)
plt.errorbar(100, 0.006901662, yerr=0.0004663945922928351)


# plt.legend()
# plt.show()
plt.savefig('msd.png')

# plt.plot(lambda_list_1, msd_list_1, 'o-', label='msd')
# plt.plot(lambda_list_2, msd_list_2, 'o-', label='msd')
