import math

def debye_integrand(x: float) -> float:
    """德拜积分被积函数，x→0时已做数值保护"""
    if abs(x) < 1e-12:
        return 0.0
    ex = math.exp(x)
    return (x**4) * ex / ((ex - 1.0) ** 2)

def trapezoid_composite(f, a: float, b: float, n: int) -> float:
    """
    复合梯形数值积分
    :param f: 被积函数
    :param a: 积分下限
    :param b: 积分上限
    :param n: 分段数
    :return: 积分结果
    """
    if n <= 0:
        raise ValueError("分段数n必须为正整数")
    h = (b - a) / n  # 计算步长
    integral = 0.5 * (f(a) + f(b))  # 首尾项
    
    # 累加中间项（系数为2）
    for i in range(1, n):
        xi = a + i * h
        integral += f(xi)
    
    return integral * h  # 最终结果

def simpson_composite(f, a: float, b: float, n: int) -> float:
    """
    复合辛普森1/3积分（强制n为偶数）
    :param f: 被积函数
    :param a: 积分下限
    :param b: 积分上限
    :param n: 分段数（必须为偶数）
    :return: 积分结果
    """
    if n <= 0 or n % 2 != 0:
        raise ValueError("辛普森积分要求分段数n为正偶数！")
    h = (b - a) / n  # 计算步长
    integral = f(a) + f(b)  # 首尾项
    
    # 累加中间项：奇数点系数4，偶数点系数2
    for i in range(1, n):
        xi = a + i * h
        if i % 2 == 1:
            integral += 4 * f(xi)
        else:
            integral += 2 * f(xi)
    
    return integral * h / 3  # 最终结果

def debye_integral(T: float, theta_d: float = 428.0, method: str = "simpson", n: int = 200) -> float:
    """
    计算德拜热容积分 I(y) = ∫₀ʸ [x⁴eˣ/(eˣ-1)²]dx, y=θ_D/T
    :param T: 温度 (K)
    :param theta_d: 德拜温度 (默认428K，金刚石)
    :param method: 积分方法 "trapezoid"/"simpson"
    :param n: 分段数
    :return: 德拜积分结果
    """
    if T <= 0:
        raise ValueError("温度T必须大于0")
    y = theta_d / T  # 积分上限
    a = 0.0
    b = y
    
    # 根据方法选择积分器
    if method.lower() == "trapezoid":
        return trapezoid_composite(debye_integrand, a, b, n)
    elif method.lower() == "simpson":
        return simpson_composite(debye_integrand, a, b, n)
    else:
        raise ValueError("方法仅支持 trapezoid 或 simpson")

# ------------------- 测试代码（验证结果） -------------------
if __name__ == "__main__":
    # 测试温度：低温、中温、高温
    test_T = [100, 300, 500]
    n_segments = 200
    
    print("=== 德拜积分计算结果（n=200）===")
    for T in test_T:
        trap = debye_integral(T, method="trapezoid", n=n_segments)
        simp = debye_integral(T, method="simpson", n=n_segments)
        print(f"T = {T} K")
        print(f"  梯形法: {trap:.6f}")
        print(f"  辛普森法: {simp:.6f}")
        print(f"  绝对误差: {abs(simp - trap):.2e}\n")
