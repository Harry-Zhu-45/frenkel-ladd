import numpy as np


def print_res(res: np.ndarray):
    """打印结果

    Args:
        res (np.ndarray): G-L 积分得到的结果数组
    """
    # print("median:", np.median(res))
    # print("mean:", np.mean(res))
    # print("std:", np.std(res, ddof=1))
    print(res)
    print(-np.mean(res), "+-", np.std(res, ddof=1))


def fl_original_analyze_func(filename: str, n: int = 10):
    """原始文献的 G-L 积分方法

    Args:
        filename (str): 文件名
        n (int, optional): G-L 积分点数. Defaults to 10.
    """
    # read file
    data = []
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("0.0"):
                data.append(float(line.split()[0]))

    set_of_data = len(data) // n
    print(set_of_data)
    data = np.array(data[:set_of_data * n])  # 截断
    data.resize(set_of_data, n)  # n个一组，不足n个的不处理
    # print(data.shape)

    # 确定 lambda 取值范围
    original_lambda_min = 0
    original_lambda_max = 632.026
    # original_lambda_max = 1774.927

    lambda_min = np.log(original_lambda_min + np.exp(3.5))
    lambda_max = np.log(original_lambda_max + np.exp(3.5))

    abscissa_n, weight_n = np.polynomial.legendre.leggauss(n)  # 横坐标, 权重
    lambda_list_n = np.exp((lambda_max - lambda_min) / 2 * abscissa_n +
                           (lambda_max + lambda_min) / 2) - np.exp(3.5)  # 生成 lambda list
    # print(lambda_list_n)

    res_list = np.sum(weight_n * data * (lambda_list_n + np.exp(3.5)), axis=1) * (lambda_max - lambda_min) / 2

    print_res(np.array(res_list))


def revisiting_analyze_func(filename: str, n: int):
    """n点G-L 积分方法，分析被读取的数据文件，计算 Delta F_MC

    Args:
        filename (str): 文件名
        n (int): G-L 积分点数.
    """
    # read file
    data = []
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("0.0"):
                data.append(float(line.split()[0]))

    set_of_data = len(data) // n
    print(set_of_data)
    data = np.array(data[:set_of_data * n])  # 截断
    data.resize(set_of_data, n)  # n个一组，不足n个的不处理
    # print(data.shape)

    # 确定 lambda 取值范围
    lambda_min = 0
    # lambda_max = 632.026
    lambda_max = 1000

    abscissa_n, weight_n = np.polynomial.legendre.leggauss(n)  # 积分点, 权重
    lambda_list_n = (lambda_max - lambda_min) / 2 * abscissa_n + (lambda_max + lambda_min) / 2  # 生成 lambda list
    # print(lambda_list_n)

    res_list = np.sum(weight_n * data, axis=1) * (lambda_max - lambda_min) / 2

    print_res(np.array(res_list))


if __name__ == "__main__":
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/3x3x6.out")
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/3x3x12.out")
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/3x4x6.out")
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/4x4x6.out")
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/4x4x12.out")
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/6x6x6.out")
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/cubic 2.out")
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/cubic 3.out")
    # fl_original_analyze_func("7360_original_box_10pt_slurm_data/cubic 4.out")

    # fl_original_analyze_func("7778_original_box_10pt_slurm_data/3x4x6.out")
    # fl_original_analyze_func("7778_original_box_10pt_slurm_data/4x6x6.out")
    # fl_original_analyze_func("7778_original_box_10pt_slurm_data/6x8x12.out")
    # fl_original_analyze_func("7778_original_box_10pt_slurm_data/8x12x12.out")

    # revisiting_analyze_func("20pt_slurm_data/cubic 3.out", 20)
    # revisiting_analyze_func("20pt_slurm_data/cubic 4.out", 20)
    revisiting_analyze_func("20pt_slurm_data/cubic 7.out", 20)
    revisiting_analyze_func("20pt_slurm_data/cubic 8.out", 20)
