#!/bin/bash

# 设置循环次数
num_runs=3

# 输出文件名
output_file="msd_data.txt"

# 清空或创建输出文件
> $output_file

# 循环运行命令并收集数据
for ((i=1; i<=$num_runs; i++)); do
    echo "Running iteration $i..."

    # 运行命令并将输出保存到临时文件
    result=$(./build/main.out -l 100 -n 100000)

    # 从输出中提取 MSD 数据
    msd=$(echo "$result" | grep "MSD:" | awk '{print $2}')

    # 将数据写入输出文件
    echo "$msd" >> $output_file
done

echo  # 换行，以确保进度条后面有一行空行
echo "Data collection complete. MSD data saved to $output_file"
