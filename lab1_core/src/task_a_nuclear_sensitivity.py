import numpy as np

def rate_3alpha(T: float) -> float:
    """
    计算3-α反应率q(T)
    :param T: 开尔文温度
    :return: 反应率q(T)
    """
    # 防止温度非正，数值保护
    if T <= 0:
        raise ValueError("温度T必须大于0！")
    T8 = T / 1.0e8  # 题目定义：T8 = T / 1e8
    return 5.09e11 * (T8 ** (-3.0)) * np.exp(-44.027 / T8)

def finite_diff_dq_dT(T0: float, h: float = 1e-8) -> float:
    """
    前向差分计算dq/dT
    :param T0: 参考温度
    :param h: 相对步长，默认1e-8
    :return: dq/dT在T0处的数值近似
    """
    delta_T = h * T0  # ✅ 关键：绝对步长 = h*T0（规避常见错误1）
    q0 = rate_3alpha(T0)
    q1 = rate_3alpha(T0 + delta_T)
    # 前向差分公式
    dq_dt = (q1 - q0) / delta_T
    return dq_dt

def sensitivity_nu(T0: float, h: float = 1e-8) -> float:
    """
    计算温度敏感性指数 nu = (T/q) * dq/dT
    :param T0: 参考温度
    :param h: 差分步长
    :return: 温度敏感性指数nu
    """
    dq_dt = finite_diff_dq_dT(T0, h)
    q0 = rate_3alpha(T0)  # ✅ 关键：必须用T0处的反应率（规避常见错误3）
    nu = (T0 / q0) * dq_dt
    return nu

def nu_table(T_values, h: float = 1e-8):
    """
    生成温度-敏感性指数列表
    :param T_values: 温度列表
    :param h: 差分步长
    :return: [(T0, nu0), (T1, nu1), ...]
    """
    table = []
    for T in T_values:
        nu = sensitivity_nu(T, h)
        table.append((T, nu))
    return table

# ------------------- 主程序：计算指定温度并输出 -------------------
if __name__ == "__main__":
    # 题目要求必算温度点
    target_T = [1.0e8, 2.5e8, 5.0e8, 1.0e9, 2.5e9, 5.0e9]
    # 计算结果
    result = nu_table(target_T)
    # 按指定格式输出
    print("3-α反应率温度敏感性指数计算结果：")
    for T, nu in result:
        print(f"{T:.3e} K : nu = {nu:.2f}")
