import numpy as np


def analyze_func(filename: str = "slurm.out") -> None:
    # 确定 lambda 取值范围
    original_lambda_min = 0
    original_lambda_max = 632.026
    lambda_min = np.log(original_lambda_min + np.e ** 3.5)
    lambda_max = np.log(original_lambda_max + np.e ** 3.5)

    # Legendre-Gauss 积分点和权重
    abscissa_10 = [
        -0.9739065285,
        -0.8650633667,
        -0.6794095683,
        -0.4333953941,
        -0.1488743390,
        0.1488743390,
        0.4333953941,
        0.6794095683,
        0.8650633667,
        0.9739065285
    ]

    weight_10 = [
        0.0666713443,
        0.1494513492,
        0.2190863625,
        0.2692667193,
        0.2955242247,
        0.2955242247,
        0.2692667193,
        0.2190863625,
        0.1494513492,
        0.0666713443
    ]

    lambda_list_10 = [
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

    # read file
    data = []
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("0.0"):
                data.append(float(line.split()[0]))

    res_list = []
    for i in range(10):
        msd_list_10 = []
        for j in range(10):
            msd_list_10.append(data[j * 10 + i])

        res = 0.0
        for i in range(10):
            res += weight_10[i] * msd_list_10[i] * (lambda_list_10[i] + np.e ** 3.5)
        res *= (lambda_max - lambda_min) / 2
        res_list.append(res)
        # print(f"Delta F_MC = {-res}")

    # errorbar of res_list
    # print("upbound:", max(res_list))
    # print("lowbound:", min(res_list))
    print("mean:", np.mean(res_list))


if __name__ == "__main__":
    analyze_func("slurm_data/3x4x6.out")
    analyze_func("slurm.jupiter.15243.out")
