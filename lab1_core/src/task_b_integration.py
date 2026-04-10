import math

def debye_integrand(x: float) -> float:
    if abs(x) < 1e-12:
        return 0.0
    ex = math.exp(x)
    return (x**4) * ex / ((ex - 1.0) ** 2)

def trapezoid_composite(f, a: float, b: float, n: int) -> float:
    # 实现复合梯形积分
    h = (b - a) / n
    s = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        s += f(a + i * h)
    return s * h

def simpson_composite(f, a: float, b: float, n: int) -> float:
    # 实现复合 Simpson 积分，并检查 n 为偶数
    if n % 2 != 0:
        raise ValueError("n must be even for Simpson's rule")
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        x = a + i * h
        if i % 2 == 1:
            s += 4 * f(x)
        else:
            s += 2 * f(x)
    return s * h / 3

def debye_integral(T: float, theta_d: float = 428.0, method: str = "simpson", n: int = 200) -> float:
    # 计算 Debye 积分 I(theta_d/T)
    y = theta_d / T
    if method == "trapezoid":
        return trapezoid_composite(debye_integrand, 0.0, y, n)
    elif method == "simpson":
        return simpson_composite(debye_integrand, 0.0, y, n)
    else:
        raise ValueError("method must be 'trapezoid' or 'simpson'")
