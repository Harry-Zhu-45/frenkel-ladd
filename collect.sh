#!/bin/bash

# 设置循环次数
num_runs=10

# 输出文件名
output_file="msd_data.txt"

# 清空或创建输出文件
>$output_file

# List of values
# values=(5.00421035 33.0590995 115.297688 299.738483 544.708655)
values=(1.3218454386828071 7.429242462300721 20.448998156279607 44.35579276774375 85.59512159890338 152.4321293389935 251.20242463993492 378.0986975526874 510.1482256352902 606.495129089744)

# 循环运行命令并收集数据
for value in "${values[@]}"; do
    echo "Running for value $value..." >>$output_file
    for ((i = 1; i <= $num_runs; i++)); do
        echo "Running iteration $i with value $value..."

        # 运行命令并将输出保存到临时文件
        result=$(./build/main.out -l "$value" -n 100000)

        # 从输出中提取 MSD 数据
        msd=$(echo "$result" | grep "MSD:" | awk '{print $2}')

        # 将数据写入输出文件
        echo "$msd," >>$output_file
    done
    echo >>$output_file
done

echo # 换行，以确保进度条后面有一行空行
echo "Data collection complete. MSD data saved to $output_file"
