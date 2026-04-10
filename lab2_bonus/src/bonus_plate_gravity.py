import numpy as np

# 万有引力常数
G = 6.674e-11

def gauss_legendre_2d(func, ax: float, bx: float, ay: float, by: float, n: int = 40) -> float:
    """
    TODO D1 已完成：二维高斯-勒让德数值积分
    积分区间：x ∈ [ax, bx], y ∈ [ay, by]
    高精度积分，适合光滑被积函数
    """
    # 获取高斯-勒让德节点和权重（默认区间 [-1,1]）
    t, w = np.polynomial.legendre.leggauss(n)
    
    # 变量替换：将 [-1,1] 映射到实际积分区间 [a,b]
    hx = (bx - ax) / 2.0
    cx = (bx + ax) / 2.0
    hy = (by - ay) / 2.0
    cy = (by + ay) / 2.0
    
    # 生成映射后的坐标
    x = hx * t + cx
    y = hy * t + cy
    
    # 二维积分：双重求和（核心公式）
    integral = 0.0
    for i in range(n):
        for j in range(n):
            integral += w[i] * w[j] * func(x[i], y[j])
    
    # 乘以变换系数
    integral *= hx * hy
    return integral

def plate_force_z(z: float, L: float = 10.0, M_plate: float = 1.0e4, m_particle: float = 1.0, n: int = 40) -> float:
    """
    TODO D2 已完成：计算方板中心正上方z处的引力z分量
    """
    # 数值保护：防止z=0除零错误
    z_safe = max(z, 1e-8)
    # 面密度
    sigma = M_plate / (L ** 2)
    # 积分区间：x,y ∈ [-L/2, L/2]
    x_min, x_max = -L/2, L/2
    y_min, y_max = -L/2, L/2
    
    # 定义被积函数
    def integrand(x, y):
        return 1.0 / (x**2 + y**2 + z_safe**2) ** 1.5
    
    # 计算二重积分
    integral = gauss_legendre_2d(integrand, x_min, x_max, y_min, y_max, n)
    
    # 核心公式计算Fz
    Fz = G * sigma * m_particle * z_safe * integral
    return Fz

def force_curve(z_values, L: float = 10.0, M_plate: float = 1.0e4, m_particle: float = 1.0, n: int = 40):
    """
    TODO D3 已完成：批量计算一组z值对应的Fz，返回数组
    """
    Fz_values = [plate_force_z(z, L, M_plate, m_particle, n) for z in z_values]
    return np.array(Fz_values)

# ------------------- 测试代码（可选，验证结果） -------------------
if __name__ == "__main__":
    # 测试：z=1m 处的引力
    z_test = 1.0
    Fz = plate_force_z(z_test)
    print(f"z = {z_test} m 时，Fz = {Fz:.6e} N")
    
    # 测试批量计算
    z_list = np.linspace(0.2, 10, 5)
    F_list = force_curve(z_list)
    print("\nz值(m) | 引力Fz(N)")
    for z, f in zip(z_list, F_list):
        print(f"{z:5.1f} | {f:.6e}")
