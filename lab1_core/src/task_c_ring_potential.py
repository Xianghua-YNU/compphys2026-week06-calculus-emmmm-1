import numpy as np
import matplotlib.pyplot as plt

def ring_potential_point(x: float, y: float, z: float, a: float = 1.0, q: float = 1.0, n_phi: int = 720) -> float:
    """
    TODO C1 已完成：用复合梯形积分计算空间单点(x,y,z)的带电圆环电势
    :param x,y,z: 空间点坐标
    :param a: 圆环半径
    :param q: 电荷参数
    :param n_phi: 角度积分分段数
    :return: 该点电势 V
    """
    # 积分变量：角度 φ 从 0 到 2π
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    d_phi = 2 * np.pi / n_phi  # 积分步长

    # 被积函数分母（距离平方+极小值防止除零）
    r2 = (x - a * np.cos(phi)) ** 2 + (y - a * np.sin(phi)) ** 2 + z ** 2
    integrand = 1.0 / np.sqrt(r2 + 1e-12)

    # 复合梯形法则积分
    integral = d_phi * (0.5 * integrand[0] + 0.5 * integrand[-1] + np.sum(integrand[1:-1]))

    # 电势公式：V = q/(2π) * 积分结果
    return q / (2 * np.pi) * integral

def ring_potential_grid(y_grid, z_grid, x0: float = 0.0, a: float = 1.0, q: float = 1.0, n_phi: int = 720):
    """
    TODO C2 已完成：在yz平面网格上计算全平面电势分布
    :param y_grid: y方向二维网格
    :param z_grid: z方向二维网格
    :param x0: 固定x坐标（默认对称面x=0）
    :return: 电势矩阵（和网格同形状）
    """
    # 初始化电势矩阵
    V_grid = np.zeros_like(y_grid, dtype=float)
    
    # 遍历所有网格点计算电势
    for i in range(y_grid.shape[0]):
        for j in range(y_grid.shape[1]):
            y = y_grid[i, j]
            z = z_grid[i, j]
            V_grid[i, j] = ring_potential_point(x0, y, z, a, q, n_phi)
    
    return V_grid

def axis_potential_analytic(z: float, a: float = 1.0, q: float = 1.0) -> float:
    """圆环轴线(z轴)解析解，用于验证数值结果"""
    return q / np.sqrt(a * a + z * z)

# ------------------- 可视化代码（实验要求：等势线+电场矢量） -------------------
def plot_ring_field(a=1.0, q=1.0, size=3.0, num=50):
    """
    绘制yz平面带电圆环的等势线+电场矢量图
    """
    # 生成网格
    y = np.linspace(-size, size, num)
    z = np.linspace(-size, size, num)
    Y, Z = np.meshgrid(y, z)
    
    # 计算电势
    V = ring_potential_grid(Y, Z, x0=0.0, a=a, q=q)
    
    # 计算电场（电场 = -电势梯度）
    Ey, Ez = np.gradient(-V, y[1]-y[0], z[1]-z[0])
    
    # 绘图
    plt.figure(figsize=(10, 8))
    # 等势线
    contour = plt.contour(Y, Z, V, levels=20, cmap='coolwarm')
    plt.clabel(contour, inline=True, fontsize=8)
    # 电场矢量
    plt.quiver(Y, Z, Ey, Ez, color='k', alpha=0.6)
    # 绘制圆环位置（x=0平面，半径a的圆）
    circle = plt.Circle((0, 0), a, color='red', fill=False, linewidth=2, label='带电圆环')
    plt.gca().add_patch(circle)
    
    plt.xlabel('y (m)')
    plt.ylabel('z (m)')
    plt.title('均匀带电圆环 yz 平面电势等势线 + 电场分布')
    plt.axis('equal')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()

# 主程序：运行可视化
if __name__ == "__main__":
    plot_ring_field()
