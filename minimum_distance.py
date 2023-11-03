"""
Author: HarryZ
Time: 2023-11-03
Function: Calculate the minimum distance between particles in a given file.
"""

import re
import numpy as np

# 打开文件
with open('fcc-3x4x6.pos', 'r', encoding='utf-8') as file:
    lines = file.readlines()

# 初始化一个列表来存储坐标信息
coordinates = []

# 定义正则表达式模式来匹配坐标信息
PATTERN = r'S0 ffff0000 ([-0-9.]+) ([-0-9.]+) ([-0-9.]+)'


# 遍历每一行并提取坐标信息
for line in lines:
    match = re.search(PATTERN, line)
    if match:
        x = float(match.group(1))
        y = float(match.group(2))
        z = float(match.group(3))
        coordinates.append((x, y, z))


# 计算粒子之间距离的最小值
min_distance = float('inf')

for i in enumerate(coordinates):
    for j in enumerate(coordinates):
        if i[0] != j[0]:
            distance = np.sqrt(
                (i[1][0] - j[1][0]) ** 2 + (i[1][1] - j[1][1]) ** 2 + (i[1][2] - j[1][2]) ** 2)
            if distance < min_distance:
                min_distance = distance
                min_distance_coord = (i[1], j[1])

print('The minimum distance between particles is: ', min_distance)
print('The coordinates of the particles are: ', min_distance_coord)
