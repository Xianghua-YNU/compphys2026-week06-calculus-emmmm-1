import numpy as np

def ring_potential_point(x: float, y: float, z: float, a: float = 1.0, q: float = 1.0, n_phi: int = 720) -> float:
    """
    用离散积分计算单点电势
    """
    phi = np.linspace(0, 2 * np.pi, n_phi)
    d_phi = 2 * np.pi / n_phi

    r2 = (x - a * np.cos(phi)) ** 2 + (y - a * np.sin(phi)) ** 2 + z ** 2
    integrand = 1.0 / np.sqrt(r2 + 1e-12)

    # 复合梯形积分
    integral = d_phi * np.sum(integrand)
    return q / (2 * np.pi) * integral

def ring_potential_grid(y_grid, z_grid, x0: float = 0.0, a: float = 1.0, q: float = 1.0, n_phi: int = 720):
    """
    ✅ 修复形状错误：
    输入：一维 y 数组、一维 z 数组
    内部：自动生成二维网格
    输出：二维电势矩阵，形状 (len(z_grid), len(y_grid))（完全匹配测试要求）
    """
    # 测试核心：自动生成二维网格
    Y, Z = np.meshgrid(y_grid, z_grid)
    V = np.zeros_like(Y, dtype=float)
    
    # 遍历所有网格点
    for i in range(Y.shape[0]):
        for j in range(Y.shape[1]):
            V[i, j] = ring_potential_point(x0, Y[i, j], Z[i, j], a, q, n_phi)
    
    return V

def axis_potential_analytic(z: float, a: float = 1.0, q: float = 1.0) -> float:
    return q / np.sqrt(a * a + z * z)
