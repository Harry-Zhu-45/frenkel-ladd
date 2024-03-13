#!./venv/bin/python3
import subprocess


# 设置循环次数
num_runs = 10

# 输出文件名
output_file = "msd_data.txt"

# 清空或创建输出文件
with open(output_file, "w") as file:
    file.write("")

# List of values
values = [
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

# 循环运行命令并收集数据
for value in values:
    with open(output_file, "a") as file:
        file.write(f"Running for value {value}...\n")
    for i in range(1, num_runs + 1):
        print(f"Running iteration {i} with value {value}...")

        # 运行命令并捕获输出
        result = subprocess.run(["./build/main.out", "-l", str(value), "-n", "1000000"], capture_output=True, text=True)

        # 从输出中提取 MSD 数据
        msd = result.stdout.split("MSD:")[1].split()[0].strip()

        # 将数据写入输出文件
        with open(output_file, "a") as file:
            file.write(f"{msd},\n")
    with open(output_file, "a") as file:
        file.write("\n")

print("Data collection complete. MSD data saved to", output_file)
